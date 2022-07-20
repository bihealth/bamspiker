/// Main module for espike app
///
/// Responsible for parsing command line parameters and global setup/cleanup.
mod conf;

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{
    bam,
    bam::{HeaderView, Read as BamRead, Record},
};

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

/// Entry point after parsing command line and reading configuration.
fn run(config: &Config) {
    // Initialize BAM reader and writer
    let mut bam_in = bam::IndexedReader::from_path(config.path_bam_in.clone()).unwrap();
    bam_in.set_threads(2).expect("Could not set thread count");
    let header = bam::Header::from_template(bam_in.header());
    let mut bam_out =
        bam::Writer::from_path(config.path_bam_out.clone(), &header, bam::Format::Bam).unwrap();
    bam_out.set_threads(2).expect("Could not set thread count");

    let header_view = HeaderView::from_header(&header);
    for tid in 0..(header_view.target_count() - 1) {
        if tid != 4 {
            continue; // TODO: for testing only
        }
        let target_name = header_view.tid2name(tid);
        let target_name = String::from_utf8_lossy(&target_name).into_owned();
        let target_len = header_view.target_len(tid).expect("Invalid target");
        bam_in.fetch(tid).expect("Could not fetch contig");

        let bar = ProgressBar::new(target_len);
        bar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "contig {msg} | {wide_bar:.cyan/blue} {pos:>7}/{len:7} [{elapsed_precise}/{eta_precise}]",
                )
                .progress_chars("##-"),
        );
        bar.set_message(target_name.clone());

        // Iterate over input file, write to output file
        let mut record = Record::new();
        let mut record_no = 0u64;
        while let Some(r) = bam_in.read(&mut record) {
            r.expect("Failed to parse record");
            let record_pos = record.pos();
            if record.is_reverse() {
                bam_out.write(&record).unwrap();
            }
            record_no += 1;
            if record_no % 10_000 == 0 {
                bar.set_position(record_pos.try_into().unwrap());
            }
        }
        bar.finish()
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
    };
    println!("Running 'espik' - your friendly BAM spiker.");
    println!("Configuration is {:?}", config);
    run(&config);
    println!("All done. Have a nice day!");
}
