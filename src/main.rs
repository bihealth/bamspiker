/// Main module for bamspiker app
///
/// Responsible for parsing command line parameters and global setup/cleanup.
mod conf;
mod helpers;

use bio::io::fasta::IndexedReader;
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rand::prelude::*;
use rust_htslib::{
    bam,
    bam::{
        record::{Cigar, CigarString},
        HeaderView, Read as BamRead, Record, Writer,
    },
};
use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet},
    io::{self, BufRead},
    process::Command,
};
use std::{fs::File, io::Write};
use tempdir::TempDir;

use conf::{VarSpec, VarType};
use helpers::{var_spec_contains, var_spec_overlaps};

static REF_WINDOW: i64 = 1000;

/// Spike variants into BAM files
#[derive(Parser)]
struct Cli {
    /// The path to the instructions YAML file
    #[clap(parse(from_os_str))]
    path_instructions: std::path::PathBuf,
    /// The path to the FAI-indexed reference file
    #[clap(parse(from_os_str))]
    path_reference: std::path::PathBuf,
    /// The path to the input BAM file
    #[clap(parse(from_os_str))]
    path_bam_in: std::path::PathBuf,
    /// The path to the output BAM file
    #[clap(parse(from_os_str))]
    path_bam_out: std::path::PathBuf,
    /// Seed for random number generator
    #[clap(short, long, value_parser, default_value_t = 42)]
    pub seed: u64,
}

/// Spike SNV into record
fn spike_into_record_snv(record: &mut Record, rng: &mut StdRng, spec: &VarSpec) {
    let cigar = record.cigar_cached().unwrap();
    let read_start_pos = cigar
        .read_pos((spec.start - 1).try_into().unwrap(), true, false)
        .unwrap();
    // We use a shifted end position that points at the last base in 0-based coordinates
    // but not behind the interval.
    let read_end_pos = cigar
        .read_pos((spec.end - 1).try_into().unwrap(), true, false)
        .unwrap();

    let y: f64 = rng.gen();
    if y >= spec.aaf {
        return;
    }

    if let (Some(start_pos), Some(end_pos_shifted)) = (read_start_pos, read_end_pos) {
        if start_pos != end_pos_shifted {
            panic!("Only SNVs are currently supported!");
        }

        let mut seq = record.seq().as_bytes();
        seq[start_pos as usize] = spec.alternative.clone().unwrap().chars().next().unwrap() as u8;

        let qname = record.qname().to_vec();
        let cigar = record.cigar().take();
        let qual = record.qual().to_vec();

        record.set(&qname, Some(&cigar), &seq, &qual)
    }
}

/// Collect reads for SNV.
///
/// This can be done in a single pass.
fn collect_reads_snv(
    spec: &VarSpec,
    tid: u32,
    bam_in: &mut bam::IndexedReader,
    denylist_file: &mut File,
    changed_writer: &mut Writer,
    rng: &mut StdRng,
) {
    let target_name =
        String::from_utf8_lossy(bam_in.header().target_names()[tid as usize]).into_owned();

    // println!("SNV {:?} - single pass", &spec);
    bam_in
        .fetch((tid, spec.start - 1 - REF_WINDOW, spec.end + REF_WINDOW))
        .expect("could not fetch region from BAM");
    let mut record = Record::new();
    while let Some(r) = bam_in.read(&mut record) {
        r.expect("Failed to parse record");
        if !record.is_proper_pair()
            || record.is_supplementary()
            || record.is_secondary()
            || record.tid() != record.mtid()
        {
            continue; // only spike into easy ones
        }

        let begin_pos = record.pos();
        record.cache_cigar();
        let end_pos = record.cigar_cached().unwrap().end_pos();

        if var_spec_overlaps(spec, &target_name, begin_pos, end_pos) {
            // Only keep spec.aaf records
            let y: f64 = rng.gen();
            if y >= spec.aaf {
                continue;
            }

            // Spike SNV into record
            spike_into_record_snv(&mut record, rng, spec);

            // Write modified record to SAM file
            changed_writer
                .write(&record)
                .expect("problem writing changed record");

            // Write qname and first/second flag to denylist file.
            let flag = if record.is_last_in_template() {
                "2"
            } else {
                "1"
            };
            let qname = String::from_utf8_lossy(record.qname());
            // println!("denied written: {} {}", &qname, flag.eq("2"));
            writeln!(denylist_file, "{} {}", qname, &flag).expect("could not write to denylist");
        }
    }
}

/// Collect reads for large deletion.
///
/// This has to be done in two passes, once over the range, then
///
/// 1. First pass - iterate over deletion range plus a padding window and:
///     a. collect read pairs where one read aligns "mostly" outside of the deletion
///         - and we see both => collect, will be considered for later spike-in
///         - and we see only one => ignore, maybe on other contig
///         - collected reads will be processed with probabily `spec.aaf` and added
///           to ignore list
///     b. collect names of read pairs that fall completely into the deletion
///         - these will be ignored in `spec.aaf` fraction of cases
/// 2. Second pass - go over collected reads and write into
///    a. go over the reads collected for spike-in
///    b. for full pairs, spike variant into reads ignoring the worst corner cases
/// 3. Second pass - iterate over range, do not write ignored records, merge in the
///    updated ones from the temporary file
///
/// Then do some post processing
///
/// 3. Sort temporary file by qname
/// 4. Load reads from temporary file, spike-in variant and write to central temporary
///    file for later merging
fn collect_reads_deletion(
    spec: &VarSpec,
    tid: u32,
    bam_in: &mut bam::IndexedReader,
    denylist_file: &mut File,
    changed_writer: &mut Writer,
    rng: &mut StdRng,
    fasta: &mut IndexedReader<File>,
) {
    // println!("Deletion {:?} - stream over region", &spec);
    let selected_records =
        collect_reads_deletion_stream_region(bam_in, tid, spec, rng, denylist_file);
    // println!("Deletion {:?} - spike-in variants", &spec);
    collect_reads_deletion_spike_in(
        bam_in,
        tid,
        spec,
        changed_writer,
        rng,
        fasta,
        &selected_records,
    );
}

/// Minimal number of bases that a read has to be before the SV start.
static MIN_DIST: i64 = 100;

/// First pass for processing a deletion.
///
/// Returns the records to consider for mutation.
fn collect_reads_deletion_stream_region(
    bam_in: &mut bam::IndexedReader,
    tid: u32,
    spec: &VarSpec,
    rng: &mut StdRng,
    denylist_file: &mut File,
) -> HashMap<(String, bool), Record> {
    let mut result = HashMap::new();

    let target_name =
        String::from_utf8_lossy(bam_in.header().target_names()[tid as usize]).into_owned();
    bam_in
        .fetch((tid, spec.start - 1 - REF_WINDOW, spec.end + REF_WINDOW))
        .expect("could not fetch region from BAM");
    let mut record = Record::new();
    while let Some(r) = bam_in.read(&mut record) {
        r.expect("Failed to parse record");
        if !record.is_proper_pair()
            || record.is_supplementary()
            || record.is_secondary()
            || record.tid() != record.mtid()
        {
            continue; // only spike into easy ones
        }

        // Get range of template on the contig, looking right from begin or left from end
        let tpl_range = match record.insert_size().cmp(&0) {
            Ordering::Less => (record.mpos(), record.cigar().end_pos()),
            Ordering::Greater => (record.pos(), record.pos() + record.insert_size()),
            Ordering::Equal => (record.pos(), record.cigar().end_pos()),
        };
        // let qname = String::from_utf8_lossy(record.qname())
        //     .to_owned()
        //     .to_string();
        // println!("tpl_range = {:?}, qname {}", &tpl_range, &qname);

        // We only consider properly aligned pairs for mutation
        if !record.is_unmapped()
            && !record.is_mate_unmapped()
            && record.tid() == record.mtid()
            && record.insert_size() != 0
        {
            if var_spec_contains(spec, &target_name, tpl_range.0, tpl_range.1) {
                // println!(" -- case contains");
                // Step (1b): proper pair, full contained in deletion.  For the first aligning mate
                // qname to deny list file with probability of `spec.aaf`.
                if record.insert_size() > 0 {
                    // write only once to deny list
                    maybe_ignore(rng, spec, &record, denylist_file);
                }
            } else if var_spec_overlaps(spec, &target_name, tpl_range.0, tpl_range.1) {
                // println!(" -- case overlap");
                // Step (1a): proper pair, not full contained in deletion but overlaps, collect
                // reads for manipulation with probability of `spec.aaf`.
                let qname = String::from_utf8_lossy(record.qname())
                    .to_owned()
                    .to_string();
                // Write only once.  Roll dice for leftmost alignment, for second one we simply
                // check the hash map.
                if record.insert_size() > 0 {
                    // Handle first alignment (by coordinate): maybe add to deny list file
                    // (is on variant haplotype).  If we do so, we add the record to the
                    // result set.  We only do this if the template overhang to the left or
                    // right of the deletion is at least MIN_DIST.
                    let good_overhang = if record.pos() <= spec.start {
                        record.pos() + MIN_DIST <= spec.start
                    } else if record.pos() > spec.start {
                        record.pos() + record.insert_size() >= spec.end - MIN_DIST
                    } else {
                        false
                    };
                    if good_overhang && maybe_ignore(rng, spec, &record, denylist_file) {
                        let key = (qname, record.insert_size() < 0);
                        // println!("inserting 1st {:?}", &key);
                        result.insert(key, record.clone());
                    }
                } else {
                    // Second alignment (by coordinate): if mate is in result then add this
                    // alignment as well.
                    let other_key = (qname, record.insert_size() > 0);
                    if result.contains_key(&other_key) {
                        let key = (other_key.0, record.insert_size() < 0);
                        // println!("inserting 2nd {:?}", &key);
                        result.insert(key, record.clone());
                    }
                }
            }
        }
    }

    result
}

/// Roll a dice and maybe add `record.qname()` to the deny list file.
///
/// Returns whether the record was added to the file
fn maybe_ignore(
    rng: &mut StdRng,
    spec: &VarSpec,
    record: &Record,
    denylist_file: &mut File,
) -> bool {
    let y: f64 = rng.gen();
    if y <= spec.aaf {
        let qname = String::from_utf8_lossy(record.qname());
        // println!("denied written: {} 1", &qname);
        // println!("denied written: {} 2", &qname);
        writeln!(denylist_file, "{} 1\n{} 2", &qname, &qname).expect("could not write to denylist");
        true
    } else {
        false
    }
}

/// Roll a dice and maybe add an error to the base
///
/// Returns whether the record was added to the file
fn maybe_mut_dna(rng: &mut StdRng, prob: f64, base: u8) -> u8 {
    static N: u8 = 78;
    static BASES: &[u8] = &[65 /*A*/, 67 /*C*/, 71 /*G*/, 84 /*T*/];

    if base == N {
        return base;
    }

    let y: f64 = rng.gen();
    if y <= prob {
        let i = rng.gen_range(0..4);
        let res = BASES[i];
        if res == base {
            BASES[(i + 1) % 4]
        } else {
            res
        }
    } else {
        base
    }
}

/// Second pass for processing a deletion.
fn collect_reads_deletion_spike_in(
    bam_in: &mut bam::IndexedReader,
    tid: u32,
    spec: &VarSpec,
    changed_writer: &mut Writer,
    rng: &mut StdRng,
    fasta: &mut IndexedReader<File>,
    selected_records: &HashMap<(String, bool), Record>,
) {
    let target_name =
        String::from_utf8_lossy(bam_in.header().target_names()[tid as usize]).into_owned();

    // Variant start position in 0-based positions
    let spec_start = spec.start - 1;

    let qnames: HashSet<String> = selected_records.keys().map(|key| key.0.clone()).collect();
    for qname in qnames {
        // Try to get first and second mate (by mapping coordinate)
        let first = selected_records.get(&(qname.clone(), false));
        let second = selected_records.get(&(qname.clone(), true));
        if let (Some(first), Some(second)) = (first, second) {
            let first = first.clone();
            let second = second.clone();
            // coordinates of template on virtual reference with deletion
            let tpl_begin = first.pos();
            let tpl_end = first.pos() + first.insert_size();
            // compute positions in `fasta` to sample from
            let (first_spec, second_spec) = if tpl_begin < spec_start {
                // Case 1: template hangs over the left side of the deletion, will
                // compute from anchor on left hand side
                //
                // Compute reference positions for read aligning to the left hand
                let first_pos = first.pos();
                let first_end = std::cmp::min(tpl_begin + first.seq_len() as i64, spec.start);
                // get start and end position of extra bases for first read
                let first_missing = first.seq_len() as i64 - (first_end - first_pos);
                let first_res = (
                    (0i64, 0i64),
                    (first_pos, first_end),
                    (spec.end, spec.end + first_missing),
                );
                // Compute bases of right read that are on either side of the deletion
                let bases_right = std::cmp::min(tpl_end - spec_start, second.seq_len() as i64);
                let bases_left = second.seq_len() as i64 - bases_right;
                // Depending on the large overlap, generate sample instructions
                let second_res = if bases_right >= bases_left {
                    let second_pos = spec.end;
                    let second_end = second_pos + bases_right;
                    let second_missing = second.seq_len() as i64 - (second_end - second_pos);
                    (
                        (spec_start - second_missing, spec_start),
                        (second_pos, second_end),
                        (0i64, 0i64),
                    )
                } else {
                    let second_pos = spec_start - bases_left;
                    let second_end = spec_start;
                    let second_missing = second.seq_len() as i64 - (second_end - second_pos);
                    (
                        (0i64, 0i64),
                        (second_pos, second_end),
                        (spec.end, spec.end + second_missing),
                    )
                };
                // resulting read sample specifications
                (first_res, second_res)
            } else {
                // Case 2: template does not hang over the left side of the deletion,
                // must hang over to the right side, will compute from anchor on the
                // right hand side
                //
                // Compute sample positions for read aligning to the right hand
                if tpl_end <= spec.end {
                    panic!("{} = tpl_end <= spec.end = {}", tpl_end, spec.end);
                }
                let second_end = first.pos() + first.insert_size();
                let second_pos = std::cmp::max(second_end - first.seq_len() as i64, spec.end);
                // get start and end position of extra bases for second read
                let second_missing = second.seq_len() as i64 - (second_end - second_pos);
                let second_res = (
                    (spec_start - second_missing, spec_start),
                    (second_pos, second_end),
                    (0i64, 0i64),
                );
                // Compute bases of left read that are on either side of the deletion
                let offset = spec.end - (second_end - first.insert_size() as i64); // into the deletion
                let bases_left = std::cmp::min(offset, second.seq_len() as i64);
                let bases_right = second.seq_len() as i64 - bases_left;
                // Dependign on the larger overlap, generate sample instructions
                let first_res = if bases_right >= bases_left {
                    let first_pos = spec.end;
                    let first_end = first_pos + bases_right;
                    let first_missing = first.seq_len() as i64 - (first_end - first_pos);
                    (
                        (spec_start - first_missing, spec_start),
                        (first_pos, first_end),
                        (0i64, 0i64),
                    )
                } else {
                    let first_pos = spec_start - offset;
                    let first_end = first_pos + bases_left;
                    let first_missing = first.seq_len() as i64 - (first_end - first_pos);
                    (
                        (0i64, 0i64),
                        (first_pos, first_end),
                        (spec.end, spec.end + first_missing),
                    )
                };
                // resulting read sample specifications
                (first_res, second_res)
            };

            // Create mutable copies of the records and assign values to them
            let mut first = first.clone();
            let mut second: Record = second.clone();
            first.set_pos(first_spec.1 .0);
            second.set_pos(second_spec.1 .0);
            first.set_insert_size(second_spec.1 .1 - first_spec.1 .0);
            second.set_insert_size(-first.insert_size());
            first.set_mpos(second.pos());
            second.set_mpos(first.pos());

            // Assign complex data fields
            // println!("qname = {}", &qname);
            // println!("  first  spec = {:?}", first_spec);
            // println!("  ssecond spec = {:?}", second_spec);
            apply_specs(first_spec, fasta, &target_name, &mut first, rng);
            apply_specs(second_spec, fasta, &target_name, &mut second, rng);

            // Write out modified records
            changed_writer
                .write(&first)
                .expect("could not write to BAM file");
            changed_writer
                .write(&second)
                .expect("could not write eto BAM file");
        }
    }
}

static PROB_ERR: &f64 = &0.0005;

fn apply_specs(
    spec: ((i64, i64), (i64, i64), (i64, i64)),
    fasta: &mut IndexedReader<File>,
    target_name: &str,
    record: &mut Record,
    rng: &mut StdRng,
) {
    // sequence
    let mut seq: Vec<u8> = Vec::new();
    let mut buf = Vec::new();
    if spec.0 .1 > spec.0 .0 {
        fasta
            .fetch(
                target_name,
                spec.0 .0.try_into().unwrap(),
                spec.0 .1.try_into().unwrap(),
            )
            .expect("could not fetch interval from FASTA");
        fasta.read(&mut buf).expect("could not read from FASTA");
        seq.append(&mut buf);
    }
    fasta
        .fetch(
            target_name,
            spec.1 .0.try_into().unwrap(),
            spec.1 .1.try_into().unwrap(),
        )
        .expect("could not fetch interval from FASTA");
    fasta.read(&mut buf).expect("could not read from FASTA");
    seq.append(&mut buf);
    if spec.2 .1 > spec.2 .0 {
        fasta
            .fetch(
                target_name,
                spec.2 .0.try_into().unwrap(),
                spec.2 .1.try_into().unwrap(),
            )
            .expect("could not fetch interval from FASTA");
        fasta.read(&mut buf).expect("could not read from FASTA");
        seq.append(&mut buf);
    }
    seq.iter_mut()
        .for_each(|c| *c = maybe_mut_dna(rng, *PROB_ERR, *c));
    // seq[start_pos as usize] = spec.alternative.clone().unwrap().chars().next().unwrap() as u8;
    // cigar string
    let leading_clipped = (spec.0 .1 - spec.0 .0).try_into().unwrap();
    let matches = (spec.1 .1 - spec.1 .0).try_into().unwrap();
    let trailing_clipped = (spec.2 .1 - spec.2 .0).try_into().unwrap();
    let mut cigar_chars: Vec<Cigar> = Vec::new();
    if leading_clipped > 0 {
        cigar_chars.push(Cigar::SoftClip(leading_clipped));
    }
    cigar_chars.push(Cigar::Match(matches));
    if trailing_clipped > 0 {
        cigar_chars.push(Cigar::SoftClip(trailing_clipped));
    }
    let cigar = CigarString(cigar_chars);
    // qname and qual remain the same
    let qname = record.qname().to_vec();
    let qual = record.qual().to_vec();
    record.set(&qname, Some(&cigar), &seq, &qual)
}

/// Collect reads from the given variant specification.
///
/// Write reads to block further down to the denylist file and write changed reads to the changed_writer.
fn collect_reads(
    spec: &VarSpec,
    tid: u32,
    bam_in: &mut bam::IndexedReader,
    denylist_file: &mut File,
    changed_writer: &mut Writer,
    rng: &mut StdRng,
    fasta: &mut IndexedReader<File>,
) {
    match spec.var_type {
        VarType::Deletion => {
            collect_reads_deletion(spec, tid, bam_in, denylist_file, changed_writer, rng, fasta)
        }
        VarType::Snv => collect_reads_snv(spec, tid, bam_in, denylist_file, changed_writer, rng),
    }
}

/// First pass over the contig
fn first_pass(
    config: &Config,
    tid: u32,
    target_name: &str,
    bam_in: &mut bam::IndexedReader,
    rng: &mut StdRng,
    tmp_dir: &TempDir,
    fasta: &mut IndexedReader<File>,
) {
    // Fetch number of variants on chromosome
    let vars_on_contig: u64 = config
        .var_specs
        .iter()
        .filter(|spec| target_name.eq(&spec.chromosome))
        .count()
        .try_into()
        .unwrap();

    // Setup progress bar
    let bar = ProgressBar::new(vars_on_contig);
    bar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{msg} | {wide_bar:.cyan/blue} {pos:>7}/{len:7} [{elapsed_precise}/{eta_precise}]",
            )
            .progress_chars("##-"),
    );
    bar.set_message(format!("{} 1st pass", &target_name));
    bar.set_position(0);

    let denylist_path = tmp_dir.path().join("qname-denylist.txt");
    let mut denylist_file = File::create(denylist_path).expect("could not create temporary file");

    let changed_path = tmp_dir.path().join("changed.sam");
    let changed_header = bam::Header::from_template(bam_in.header());
    let mut changed_writer = Writer::from_path(changed_path, &changed_header, bam::Format::Sam)
        .expect("could not create temporary writer");

    // Iterate over input file, write to output file
    // let update_bar = |pos: u64| bar.set_position(pos);
    for spec in &config.var_specs {
        if target_name.eq(&spec.chromosome) {
            collect_reads(
                spec,
                tid,
                bam_in,
                &mut denylist_file,
                &mut changed_writer,
                rng,
                fasta,
            );
            bar.inc(1)
        }
    }
    bar.finish_with_message(format!("{} 1st pass done", &target_name));
}

/// First pass over the contig
fn second_pass(
    _config: &Config,
    tid: u32,
    target_len: u64,
    bam_in: &mut bam::IndexedReader,
    bam_out: &mut Writer,
    _rng: &mut StdRng,
    tmp_dir: &TempDir,
) {
    let target_name =
        String::from_utf8_lossy(bam_in.header().target_names()[tid as usize]).into_owned();

    let bar = ProgressBar::new(target_len);
    bar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{msg} | {wide_bar:.cyan/blue} {pos:>7}/{len:7} [{elapsed_precise}/{eta_precise}]",
            )
            .progress_chars("##-"),
    );
    bar.set_message(format!("{} 2nd pass", &target_name));
    bar.set_position(0);

    bam_in
        .fetch(tid)
        .expect("could not fetch whole contig from BAM");

    // Load deny list of reads to ignore from input
    let denylist_path = tmp_dir.path().join("qname-denylist.txt");
    let file = File::open(denylist_path).expect("could not open deny list file");
    let deny_qnames = {
        let mut deny_qnames = HashSet::new();
        for line in io::BufReader::new(file).lines() {
            let line = line.expect("reading line failed");
            let (qname, flag) = line.rsplit_once(' ').unwrap();
            // println!("denied loaded: {} {}", &qname, flag.eq("2"));
            deny_qnames.insert((qname.to_owned(), flag.eq("2")));
        }
        deny_qnames
    };

    // Open temporary BAM file with reads to merge into output
    let path_tmp = tmp_dir
        .path()
        .join("sorted.bam")
        .into_os_string()
        .into_string()
        .unwrap();
    let mut bam_tmp = bam::Reader::from_path(&path_tmp).expect("could not open temporary file");
    let mut next_tmp = true;
    let mut record_tmp = Record::new();
    // Try to read first record right now
    if let Some(r) = bam_tmp.read(&mut record_tmp) {
        r.expect("Failed to parse record");
    } else {
        next_tmp = false;
    }

    // Stream over contig suppress deny-listed records, write mutated ones
    let mut prev_pos = -1i64;
    let mut record = Record::new();
    let mut record_no = 0;
    while let Some(r) = bam_in.read(&mut record) {
        r.expect("Failed to parse record");

        // Ignore non-primary records
        if !record.is_proper_pair()
            || record.is_supplementary()
            || record.is_secondary()
            || record.tid() != record.mtid()
        {
            continue;
        }

        // Merge from temporary file if the records fit
        while next_tmp && record_tmp.pos() >= prev_pos && record_tmp.pos() <= record.pos() {
            // println!(
            //     "from tmp        {} {}",
            //     String::from_utf8_lossy(record_tmp.qname()),
            //     record.is_last_in_template()
            // );
            bam_out.write(&record_tmp).unwrap();

            if let Some(r) = bam_tmp.read(&mut record_tmp) {
                r.expect("Failed to parse record");
            } else {
                next_tmp = false;
            }
        }

        // Write records that are not on the deny list
        let qname = String::from_utf8_lossy(record.qname()).into_owned();
        let key = (qname, record.is_last_in_template());
        if !deny_qnames.contains(&key) {
            // println!(
            //     "allowed writing {} {}",
            //     &key.0,
            //     record.is_last_in_template()
            // );
            bam_out.write(&record).unwrap();
        } else {
            // println!(
            //     "denied writing  {} {}",
            //     &key.0,
            //     record.is_last_in_template()
            // );
        }

        prev_pos = record.pos();
        record_no += 1;
        if record_no % 10_000 == 0 {
            bar.set_position(record.pos().try_into().unwrap());
        }
    }
    bar.finish_with_message(format!("{} 2nd pass done", &target_name));
}

/// Entry point after parsing command line and reading configuration.
fn run(config: &Config) {
    // Initialize PRNG
    let mut rng = StdRng::seed_from_u64(config.seed);

    // Open FASTA file
    let mut fasta =
        IndexedReader::from_file(&config.path_reference).expect("Could not open reference file");

    // Initialize BAM reader and writer
    let mut bam_in = bam::IndexedReader::from_path(config.path_bam_in.clone()).unwrap();
    bam_in.set_threads(2).expect("Could not set thread count");
    let header = bam::Header::from_template(bam_in.header());
    let mut bam_out =
        Writer::from_path(config.path_bam_out.clone(), &header, bam::Format::Bam).unwrap();
    bam_out.set_threads(2).expect("Could not set thread count");

    let header_view = HeaderView::from_header(&header);
    for tid in 0..(header_view.target_count()) {
        // Fetch meta data about file
        let target_name = header_view.tid2name(tid);
        let target_name = String::from_utf8_lossy(target_name).into_owned();
        let target_len = header_view.target_len(tid).expect("Invalid target");

        let tmp_dir = TempDir::new("bamspiker").expect("could not create temporary directory");

        // Perform first pass over contig and collect reads that need modification
        first_pass(
            config,
            tid,
            &target_name,
            &mut bam_in,
            &mut rng,
            &tmp_dir,
            &mut fasta,
        );

        // Sort temporary read file by coordinate (need to convert to BAM, samtools sort does not write
        // header to SAM file)
        // println!("sort reads by coordinate");
        let path_changed = tmp_dir
            .path()
            .join("changed.sam")
            .into_os_string()
            .into_string()
            .unwrap();
        let path_sorted = tmp_dir
            .path()
            .join("sorted.bam")
            .into_os_string()
            .into_string()
            .unwrap();
        // let call = format!("samtools sort -OBAM -o {} {}", &path_sorted, &path_changed);
        // println!("Executing {}", &call);
        let _status = Command::new("samtools")
            .args(["sort", "-OBAM", "-o", &path_sorted, &path_changed])
            .status()
            .expect("samtools call failed");
        // println!("-> result = {:?}", &status);

        // pause();

        // Perform second pass over contig
        second_pass(
            config,
            tid,
            target_len,
            &mut bam_in,
            &mut bam_out,
            &mut rng,
            &tmp_dir,
        );
    }
}

use conf::{load_var_specs, Config};

fn main() {
    let args = Cli::parse();
    let config: Config = Config {
        path_reference: args.path_reference,
        path_bam_in: args.path_bam_in,
        path_bam_out: args.path_bam_out,
        var_specs: load_var_specs(&args.path_instructions).var_specs,
        seed: args.seed,
    };
    println!("Running 'bamspiker'...");
    println!("Configuration is {:?}", config);
    run(&config);
    println!("All done. Have a nice day!");
}
