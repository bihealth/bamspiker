#![allow(unused)]

use std::str;
use std::{fs::File, io::Read};

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{
    bam,
    bam::{HeaderView, Read as BamRead},
};
use serde::{Deserialize, Serialize};

/// Specify spike-in of a small variant
#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct SmallVarSpec {
    /// Chromosome as in BAM file
    chromosome: String,
    /// 1-based start position
    start: i32,
    /// 1-based end position
    end: i32,
    /// Reference bases
    reference: String,
    /// Alternative bases
    alternative: String,
    /// Alternate allele fraction
    aaf: f64,
}

/// Define type for a structural variant
#[derive(Debug, PartialEq, Serialize, Deserialize)]
enum SvType {
    /// A deletion (drop in coverage, discordant reads, split reads)
    Deletion,
}

/// Specify spike-in of a structural variant
#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct StrucVarSpec {
    /// The type of the structural variant
    sv_type: SvType,
    /// Chromosome as in BAM file
    chromosome: String,
    /// 1-based start position
    start: i32,
    /// 1-based end position
    end: i32,
    /// Specify alternate allele fraction
    aaf: f64,
}

/// Instructions for spike-ins
#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct Instructions {
    /// Small variants to spike in
    small_vars: Vec<SmallVarSpec>,
    /// Structural variants to spike in
    svs: Vec<StrucVarSpec>,
}

/// Overall configuration of espike
#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct Config {
    /// Path to the FAI indexed FASTA reference file
    path_reference: std::path::PathBuf,
    /// Path to the input BAM file
    path_bam_in: std::path::PathBuf,
    /// Path to the output BAM file
    path_bam_out: std::path::PathBuf,
    /// Simulation instructions
    instructions: Instructions,
}

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
}

/// Load instructions from the YAML `file` at the given path.
fn load_instructions(file: &std::path::PathBuf) -> Instructions {
    let mut file = File::open(file).expect("Unable to open file");
    let mut contents = String::new();

    file.read_to_string(&mut contents)
        .expect("Unable to read file");

    let result: Instructions =
        serde_yaml::from_str(&contents).expect("Unable to parse file contents");
    result
}

fn run(config: &Config) {
    // Initialize BAM reader and writer
    let mut bam_in = bam::IndexedReader::from_path(config.path_bam_in.clone()).unwrap();
    bam_in.set_threads(2);
    let header = bam::Header::from_template(bam_in.header());
    let mut bam_out =
        bam::Writer::from_path(config.path_bam_out.clone(), &header, bam::Format::Bam).unwrap();
    bam_out.set_threads(2);

    let header_view = HeaderView::from_header(&header);
    for tid in (0..(header_view.target_count() - 1)) {
        if (tid != 4) {
            continue; // TODO: for testing only
        }
        let target_name = header_view.tid2name(tid);
        let target_len = header_view.target_len(tid).expect("Invalid target");
        bam_in.fetch(tid).expect("Could not fetch contig");

        println!(
            "Processing contig {:?} ({})",
            str::from_utf8(&target_name).unwrap(),
            tid
        );

        let bar = ProgressBar::new(target_len);
        bar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "[{elapsed_precise}/{eta_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
                )
                .progress_chars("##-"),
        );

        // Iterate over input file, write to output file
        for r in bam_in.records() {
            let record = r.unwrap();
            let record_pos = record.pos();
            bar.set_position(record_pos.try_into().unwrap());
            if record.is_reverse() {
                bam_out.write(&record).unwrap();
            }
        }
        bar.set_position(target_len);
    }
}

fn main() {
    let args = Cli::parse();
    let config: Config = Config {
        path_reference: args.path_reference,
        path_bam_in: args.path_bam_in,
        path_bam_out: args.path_bam_out,
        instructions: load_instructions(&args.path_instructions),
    };
    println!("Running 'espik' - your friendly BAM spiker.");
    println!("Configuration is {:?}", config);
    run(&config);
    println!("All done. Have a nice day!");
}
