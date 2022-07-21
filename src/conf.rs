/// Configuration for spiking in variants
use serde::{Deserialize, Serialize};
use std::{fs::File, io::Read};

/// Define type for a structural variant
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub enum VarType {
    /// A single nucleotide variant
    Snv,
    /// A deletion (drop in coverage, discordant reads, split reads)
    Deletion,
}

/// Specify spike-in of a variant
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct VarSpec {
    /// The variant type
    pub var_type: VarType,
    /// Chromosome as in BAM file
    pub chromosome: String,
    /// 1-based start position
    pub start: i64,
    /// 1-based end position
    pub end: i64,
    /// Reference bases
    pub reference: Option<String>,
    /// Alternative bases
    pub alternative: Option<String>,
    /// Alternate allele fraction
    pub aaf: f64,
}

/// A collection of variant specifications
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct VarSpecs {
    /// The list of variant specifications
    pub var_specs: Vec<VarSpec>,
}

/// Overall configuration of bamspiker
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct Config {
    /// Path to the FAI indexed FASTA reference file
    pub path_reference: std::path::PathBuf,
    /// Path to the input BAM file
    pub path_bam_in: std::path::PathBuf,
    /// Path to the output BAM file
    pub path_bam_out: std::path::PathBuf,
    /// Simulation instructions
    pub var_specs: Vec<VarSpec>,
    /// Seed for random number generator
    pub seed: u64,
}

/// Load variant specifications from the YAML `file` at the given path.
pub fn load_var_specs(file: &std::path::PathBuf) -> VarSpecs {
    let mut file = File::open(file).expect("Unable to open file");
    let mut contents = String::new();

    file.read_to_string(&mut contents)
        .expect("Unable to read file");

    let result: VarSpecs = serde_yaml::from_str(&contents).expect("Unable to parse file contents");
    result
}
