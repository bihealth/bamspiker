#![allow(unused)]

use std::{fs::File, io::Read};

use clap::Parser;
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
    return result;
}

fn main() {
    let args = Cli::parse();
    let config: Config = Config {
        path_bam_in: args.path_bam_in,
        path_bam_out: args.path_bam_out,
        instructions: load_instructions(&args.path_instructions),
    };
    println!("Running 'espik' - your BAM spiker.");
    println!("Configuration is {:?}", config);
    println!("All done. Have a nice day!");
}
