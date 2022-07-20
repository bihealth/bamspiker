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
    bam::{HeaderView, Read as BamRead, Record},
};
use std::fs::File;

use conf::{SmallVarSpec, StrucVarSpec, VarSpec};
use helpers::{small_var_spec_overlaps, spec_overlaps, struc_var_spec_overlaps};

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

/// Fetch overlapping specification from configuration for the given read alignment.
fn fetch_spec(config: &Config, target_name: &str, begin_pos: i64, end_pos: i64) -> Option<VarSpec> {
    // println!("{:?} {} {}", &target_name, begin_pos, end_pos);
    let small_var_spec = config
        .instructions
        .small_vars
        .iter()
        .find(|&var_spec| small_var_spec_overlaps(var_spec, &target_name, begin_pos, end_pos));
    if let Some(inner) = small_var_spec {
        return Some(VarSpec::SmallVar(inner.clone()));
    }

    let struc_var_spec = config
        .instructions
        .svs
        .iter()
        .find(|&var_spec| struc_var_spec_overlaps(var_spec, &target_name, begin_pos, end_pos));
    if let Some(inner) = struc_var_spec {
        Some(VarSpec::StrucVar(inner.clone()))
    } else {
        None
    }
}

/// Apply small variant to read record.
fn apply_small_var_spec(record: &mut Record, rng: &mut StdRng, spec: &SmallVarSpec) {
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
        // let seq = record.seq().as_bytes();
        // println!("pos = {:?}, from = {:?}, to = {:?}", start_pos, seq[start_pos as usize] as char, spec.alternative.chars().nth(0).unwrap());
        // println!("SEQ BEFORE: {:?}", std::str::from_utf8(&seq));

        let mut seq = record.seq().as_bytes();
        seq[start_pos as usize] = spec.alternative.chars().nth(0).unwrap() as u8;

        let qname = record.qname().to_vec();
        let cigar = record.cigar().clone().take();
        let qual = record.qual().to_vec();

        record.set(&qname, Some(&cigar), &seq, &qual)
    }

    // let seq = record.seq().as_bytes();
    // println!("SEQ AFTER:  {:?}", std::str::from_utf8(&seq));
}

/// Apply structural variant to read record.
fn apply_struc_var_spec(
    record: &mut Record,
    _rng: &mut StdRng,
    spec: &StrucVarSpec,
    _fasta: &IndexedReader<File>,
) {
    let cigar = record.cigar_cached().unwrap();
    let read_start_pos = cigar
        .read_pos((spec.start - 1).try_into().unwrap(), true, false)
        .unwrap();
    // We use a shifted end position that points at the last base in 0-based coordinates
    // but not behind the interval.
    let read_end_pos = cigar
        .read_pos((spec.end - 1).try_into().unwrap(), true, false)
        .unwrap();


}

/// Apply variant into read record.
fn apply_variant(
    record: &mut Record,
    rng: &mut StdRng,
    var_spec: &VarSpec,
    fasta: &IndexedReader<File>,
) {
    // println!("applying changes of {:?} to {:?}", &var_spec, &record);
    match var_spec {
        VarSpec::SmallVar(inner_spec) => apply_small_var_spec(record, rng, inner_spec),
        VarSpec::StrucVar(inner_spec) => apply_struc_var_spec(record, rng, inner_spec, fasta),
    }
}

/// Process one contig (the call to `fetch` occurs before this function is called).
fn process_contig<UpdateProgressBar>(
    config: &Config,
    target_name: &str,
    bam_in: &mut bam::IndexedReader,
    bam_out: &mut bam::Writer,
    rng: &mut StdRng,
    fasta: &IndexedReader<File>,
    update_bar: UpdateProgressBar,
) where
    UpdateProgressBar: Fn(u64) -> (),
{
    let mut record = Record::new();
    let mut record_no = 0u64;
    let mut curr_spec: Option<VarSpec> = None;
    while let Some(r) = bam_in.read(&mut record) {
        r.expect("Failed to parse record");

        let begin_pos = record.pos();
        record.cache_cigar();
        let end_pos = record.cigar_cached().unwrap().end_pos();

        curr_spec = match curr_spec {
            None => fetch_spec(&config, &target_name, begin_pos, end_pos),
            Some(spec) => {
                if spec_overlaps(&spec, &target_name, begin_pos, end_pos) {
                    Some(spec)
                } else {
                    fetch_spec(&config, &target_name, begin_pos, end_pos)
                }
            }
        };

        // Apply variant specification to read, if any. Then, write out.
        if let Some(inner_spec) = &curr_spec {
            apply_variant(&mut record, rng, inner_spec, fasta);
        }
        bam_out.write(&record).unwrap();

        record_no += 1;
        if record_no % 10_000 == 0 {
            let record_pos = record.pos();
            update_bar(record_pos.try_into().unwrap())
        }
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
        bam::Writer::from_path(config.path_bam_out.clone(), &header, bam::Format::Bam).unwrap();
    bam_out.set_threads(2).expect("Could not set thread count");

    let header_view = HeaderView::from_header(&header);
    for tid in 0..(header_view.target_count()) {
        // Fetch meta data about file
        let target_name = header_view.tid2name(tid);
        let target_name = String::from_utf8_lossy(&target_name).into_owned();
        let target_len = header_view.target_len(tid).expect("Invalid target");

        // Jump to position in BAM file
        bam_in.fetch(tid).expect("Could not fetch contig");

        // Setup progress bar
        let bar = ProgressBar::new(target_len);
        bar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{msg} | {wide_bar:.cyan/blue} {pos:>7}/{len:7} [{elapsed_precise}/{eta_precise}]",
                )
                .progress_chars("##-"));
        bar.set_message(format!("{} done", target_name.clone()));
        bar.set_position(0);

        // Iterate over input file, write to output file
        let update_bar = |pos: u64| bar.set_position(pos);
        process_contig(
            &config,
            &target_name,
            &mut bam_in,
            &mut bam_out,
            &mut rng,
            &fasta,
            update_bar,
        );
        bar.finish_with_message(format!("{} done", target_name.clone()));
    }
}

use conf::{load_instructions, Config};

fn main() {
    let args = Cli::parse();
    let config: Config = Config {
        path_reference: args.path_reference,
        path_bam_in: args.path_bam_in,
        path_bam_out: args.path_bam_out,
        instructions: load_instructions(&args.path_instructions),
        seed: args.seed,
    };
    println!("Running 'espik' - your friendly BAM spiker.");
    println!("Configuration is {:?}", config);
    run(&config);
    println!("All done. Have a nice day!");
}
