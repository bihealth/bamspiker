/// Main module for espike app
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
    bam::{HeaderView, Read as BamRead, Record, Writer},
};
use std::{
    collections::HashSet,
    io::{self, BufRead},
    process::Command,
};
use std::{fs::File, io::Write};
use tempdir::TempDir;

use conf::{VarSpec, VarType};
use helpers::{pause, var_spec_overlaps};

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

// /// Fetch overlapping specification from configuration for the given read alignment.
// fn fetch_spec(config: &Config, target_name: &str, begin_pos: i64, end_pos: i64) -> Option<VarSpec> {
//     // println!("{:?} {} {}", &target_name, begin_pos, end_pos);
//     let small_var_spec = config
//         .instructions
//         .small_vars
//         .iter()
//         .find(|&var_spec| small_var_spec_overlaps(var_spec, &target_name, begin_pos, end_pos));
//     if let Some(inner) = small_var_spec {
//         return Some(VarSpec::SmallVar(inner.clone()));
//     }

//     let struc_var_spec = config
//         .instructions
//         .svs
//         .iter()
//         .find(|&var_spec| struc_var_spec_overlaps(var_spec, &target_name, begin_pos, end_pos));
//     if let Some(inner) = struc_var_spec {
//         Some(VarSpec::StrucVar(inner.clone()))
//     } else {
//         None
//     }
// }

// /// Apply small variant to read record.
// fn apply_small_var_spec(record: &mut Record, rng: &mut StdRng, spec: &SmallVarSpec) {
//     let cigar = record.cigar_cached().unwrap();
//     let read_start_pos = cigar
//         .read_pos((spec.start - 1).try_into().unwrap(), true, false)
//         .unwrap();
//     // We use a shifted end position that points at the last base in 0-based coordinates
//     // but not behind the interval.
//     let read_end_pos = cigar
//         .read_pos((spec.end - 1).try_into().unwrap(), true, false)
//         .unwrap();

//     let y: f64 = rng.gen();
//     if y >= spec.aaf {
//         return;
//     }

//     if let (Some(start_pos), Some(end_pos_shifted)) = (read_start_pos, read_end_pos) {
//         if start_pos != end_pos_shifted {
//             panic!("Only SNVs are currently supported!");
//         }
//         // let seq = record.seq().as_bytes();
//         // println!("pos = {:?}, from = {:?}, to = {:?}", start_pos, seq[start_pos as usize] as char, spec.alternative.chars().nth(0).unwrap());
//         // println!("SEQ BEFORE: {:?}", std::str::from_utf8(&seq));

//         let mut seq = record.seq().as_bytes();
//         seq[start_pos as usize] = spec.alternative.chars().nth(0).unwrap() as u8;

//         let qname = record.qname().to_vec();
//         let cigar = record.cigar().clone().take();
//         let qual = record.qual().to_vec();

//         record.set(&qname, Some(&cigar), &seq, &qual)
//     }

//     // let seq = record.seq().as_bytes();
//     // println!("SEQ AFTER:  {:?}", std::str::from_utf8(&seq));
// }

/// Apply structural variant to read record.
// fn apply_struc_var_spec(
//     record: &mut Record,
//     _rng: &mut StdRng,
//     spec: &StrucVarSpec,
//     _fasta: &IndexedReader<File>,
// ) {
//     let cigar = record.cigar_cached().unwrap();
//     let read_start_pos = cigar
//         .read_pos((spec.start - 1).try_into().unwrap(), true, false)
//         .unwrap();
//     // We use a shifted end position that points at the last base in 0-based coordinates
//     // but not behind the interval.
//     let read_end_pos = cigar
//         .read_pos((spec.end - 1).try_into().unwrap(), true, false)
//         .unwrap();
// }

// /// Apply variant into read record.
// fn apply_variant(
//     record: &mut Record,
//     rng: &mut StdRng,
//     var_spec: &VarSpec,
//     fasta: &IndexedReader<File>,
// ) {
//     // println!("applying changes of {:?} to {:?}", &var_spec, &record);
//     match var_spec {
//         VarSpec::SmallVar(inner_spec) => apply_small_var_spec(record, rng, inner_spec),
//         VarSpec::StrucVar(inner_spec) => apply_struc_var_spec(record, rng, inner_spec, fasta),
//     }
// }

// /// Process one contig (the call to `fetch` occurs before this function is called).
// fn process_contig<UpdateProgressBar>(
//     config: &Config,
//     target_name: &str,
//     bam_in: &mut bam::IndexedReader,
//     bam_out: &mut Writer,
//     rng: &mut StdRng,
//     fasta: &IndexedReader<File>,
//     update_bar: UpdateProgressBar,
// ) where
//     UpdateProgressBar: Fn(u64) -> (),
// {
//     let mut record = Record::new();
//     let mut record_no = 0u64;
//     let mut curr_spec: Option<VarSpec> = None;
//     while let Some(r) = bam_in.read(&mut record) {
//         r.expect("Failed to parse record");

//         let begin_pos = record.pos();
//         record.cache_cigar();
//         let end_pos = record.cigar_cached().unwrap().end_pos();

//         curr_spec = match curr_spec {
//             None => fetch_spec(&config, &target_name, begin_pos, end_pos),
//             Some(spec) => {
//                 if spec_overlaps(&spec, &target_name, begin_pos, end_pos) {
//                     Some(spec)
//                 } else {
//                     fetch_spec(&config, &target_name, begin_pos, end_pos)
//                 }
//             }
//         };

//         // Apply variant specification to read, if any. Then, write out.
//         if let Some(inner_spec) = &curr_spec {
//             apply_variant(&mut record, rng, inner_spec, fasta);
//         }
//         bam_out.write(&record).unwrap();

//         record_no += 1;
//         if record_no % 10_000 == 0 {
//             let record_pos = record.pos();
//             update_bar(record_pos.try_into().unwrap())
//         }
//     }
// }

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
        seq[start_pos as usize] = spec.alternative.clone().unwrap().chars().nth(0).unwrap() as u8;

        let qname = record.qname().to_vec();
        let cigar = record.cigar().clone().take();
        let qual = record.qual().to_vec();

        record.set(&qname, Some(&cigar), &seq, &qual)
    }
}

/// Collect reads for SNV.
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

    println!("SNV {:?} - single pass", &spec);
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

        if var_spec_overlaps(&spec, &target_name, begin_pos, end_pos) {
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
            println!("denied written: {} {}", &qname, flag.eq("2"));
            writeln!(denylist_file, "{} {}", qname, &flag).expect("could not write to denylist");
        }
    }
}

/// Collect reads for large deletion.
fn collect_reads_deletion(
    _spec: &VarSpec,
    _tid: u32,
    _bam_in: &mut bam::IndexedReader,
    _denylist_file: &mut File,
    _changed_writer: &mut Writer,
    _rng: &mut StdRng,
) {
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
) {
    match spec.var_type {
        VarType::Deletion => {
            collect_reads_deletion(spec, tid, bam_in, denylist_file, changed_writer, rng)
        }
        VarType::Snv => collect_reads_snv(spec, tid, bam_in, denylist_file, changed_writer, rng),
    }
}

/// First pass over the contig
fn first_pass(
    config: &Config,
    tid: u32,
    _target_len: u64,
    target_name: &str,
    bam_in: &mut bam::IndexedReader,
    rng: &mut StdRng,
    tmp_dir: &TempDir,
    _fasta: &IndexedReader<File>,
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
    bar.set_message(format!("{} 1st pass", target_name.clone()));
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
            );
            bar.inc(1)
        }
    }
    bar.finish_with_message(format!("{} 1st pass done", target_name.clone()));
}

/// First pass over the contig
fn second_pass(
    _config: &Config,
    tid: u32,
    _target_len: u64,
    _target_name: &str,
    bam_in: &mut bam::IndexedReader,
    bam_out: &mut Writer,
    _rng: &mut StdRng,
    tmp_dir: &TempDir,
) {
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
            let (qname, flag) = line.rsplit_once(" ").unwrap();
            println!("denied loaded: {} {}", &qname, flag.eq("2"));
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
            println!(
                "from tmp        {} {}",
                String::from_utf8_lossy(record_tmp.qname()),
                record.is_last_in_template()
            );
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
            println!(
                "allowed writing {} {}",
                &key.0,
                record.is_last_in_template()
            );
            bam_out.write(&record).unwrap();
        } else {
            println!(
                "denied writing  {} {}",
                &key.0,
                record.is_last_in_template()
            );
        }

        prev_pos = record.pos();
    }
}

/// Entry point after parsing command line and reading configuration.
fn run(config: &Config) {
    // Initialize PRNG
    let mut rng = StdRng::seed_from_u64(config.seed);

    // Open FASTA file
    let fasta =
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
        let target_name = String::from_utf8_lossy(&target_name).into_owned();
        let target_len = header_view.target_len(tid).expect("Invalid target");

        let tmp_dir = TempDir::new("espike").expect("could not create temporary directory");

        // Perform first pass over contig and collect reads that need modification
        first_pass(
            config,
            tid,
            target_len,
            &target_name,
            &mut bam_in,
            &mut rng,
            &tmp_dir,
            &fasta,
        );

        // Sort temporary read file by coordinate (need to convert to BAM, samtools sort does not write
        // header to SAM file)
        println!("sort reads by coordinate");
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
        let call = format!("samtools view -OBAM -o {} {}", &path_sorted, &path_changed);
        println!("Executing {}", &call);
        let status = Command::new("samtools")
            .args(["view", "-OBAM", "-o", &path_sorted, &path_changed])
            .status()
            .expect("samtools call failed");
        println!("-> result = {:?}", &status);

        pause();

        // Perform second pass over contig
        second_pass(
            config,
            tid,
            target_len,
            &target_name,
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
    println!("Running 'espike' - your friendly BAM spiker.");
    println!("Configuration is {:?}", config);
    run(&config);
    println!("All done. Have a nice day!");
}
