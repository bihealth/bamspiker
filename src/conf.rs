/// Configuration for spiking in variants
use serde::{Deserialize, Serialize};
use std::{fs::File, io::Read};

/// Specify spike-in of a small variant
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct SmallVarSpec {
    /// Chromosome as in BAM file
    pub chromosome: String,
    /// 1-based start position
    pub start: i64,
    /// 1-based end position
    pub end: i64,
    /// Reference bases
    pub reference: String,
    /// Alternative bases
    pub alternative: String,
    /// Alternate allele fraction
    pub aaf: f64,
}

/// Define type for a structural variant
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub enum SvType {
    /// A deletion (drop in coverage, discordant reads, split reads)
    Deletion,
}

/// Specify spike-in of a structural variant
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct StrucVarSpec {
    /// The type of the structural variant
    pub sv_type: SvType,
    /// Chromosome as in BAM file
    pub chromosome: String,
    /// 1-based start position
    pub start: i64,
    /// 1-based end position
    pub end: i64,
    /// Specify alternate allele fraction
    pub aaf: f64,
}

/// Specify either a small or a structural variant
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub enum VarSpec {
    /// Small variant with specification
    SmallVar(SmallVarSpec),
    /// Structural variant with specification
    StrucVar(StrucVarSpec),
}

/// Instructions for spike-ins
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct Instructions {
    /// Small variants to spike in
    pub small_vars: Vec<SmallVarSpec>,
    /// Structural variants to spike in
    pub svs: Vec<StrucVarSpec>,
}

/// Overall configuration of espike
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct Config {
    /// Path to the FAI indexed FASTA reference file
    pub path_reference: std::path::PathBuf,
    /// Path to the input BAM file
    pub path_bam_in: std::path::PathBuf,
    /// Path to the output BAM file
    pub path_bam_out: std::path::PathBuf,
    /// Simulation instructions
    pub instructions: Instructions,
    /// Seed for random number generator
    pub seed: u64,
}

/// Load instructions from the YAML `file` at the given path.
pub fn load_instructions(file: &std::path::PathBuf) -> Instructions {
    let mut file = File::open(file).expect("Unable to open file");
    let mut contents = String::new();

    file.read_to_string(&mut contents)
        .expect("Unable to read file");

    let result: Instructions =
        serde_yaml::from_str(&contents).expect("Unable to parse file contents");
    result
}
