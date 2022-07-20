/// Helper functions
use super::conf::VarSpec;

use std::io;
use std::io::prelude::*;

/// Return whether two intervals overlap.
fn coords_overlap(x1: i64, x2: i64, y1: i64, y2: i64) -> bool {
    x1 < y2 && y1 < x2
}

pub fn var_spec_overlaps(
    var_spec: &VarSpec,
    itv_contig: &str,
    itv_begin: i64,
    itv_end: i64,
) -> bool {
    itv_contig.eq(&var_spec.chromosome)
        && coords_overlap(itv_begin, itv_end, var_spec.start - 1, var_spec.end)
}

#[cfg(test)]
mod tests {
    use super::super::conf::{VarSpec, VarType};
    use super::{coords_overlap, var_spec_overlaps};

    #[test]
    fn test_coords_overlap() {
        assert!(!coords_overlap(0, 100, 100, 200));
        assert!(!coords_overlap(100, 200, 0, 100));
        assert!(coords_overlap(0, 101, 100, 200));
        assert!(coords_overlap(0, 200, 150, 160));
    }

    #[test]
    fn test_var_spec_overlaps_snv() {
        let var = VarSpec {
            var_type: VarType::Snv,
            chromosome: "chr5".to_string(),
            start: 100,
            end: 100,
            reference: Some("A".to_string()),
            alternative: Some("T".to_string()),
            aaf: 0.5,
        };

        assert!(!var_spec_overlaps(&var, "chr5", 98, 99));
        assert!(var_spec_overlaps(&var, "chr5", 99, 100));
        assert!(!var_spec_overlaps(&var, "chr5", 100, 101));

        assert!(!var_spec_overlaps(&var, "chr4", 99, 100));
    }

    #[test]
    fn test_var_spec_overlaps_deletion() {
        let var = VarSpec {
            var_type: VarType::Deletion,
            chromosome: "chr5".to_string(),
            start: 1000,
            end: 2000,
            reference: None,
            alternative: None,
            aaf: 0.5,
        };

        assert!(!var_spec_overlaps(&var, "chr5", 900, 999));
        assert!(var_spec_overlaps(&var, "chr5", 900, 1000));
        assert!(var_spec_overlaps(&var, "chr5", 1999, 2100));
        assert!(!var_spec_overlaps(&var, "chr5", 2000, 2100));

        assert!(!var_spec_overlaps(&var, "chr4", 999, 2000));
    }
}

pub fn pause() {
    let mut stdin = io::stdin();
    let mut stdout = io::stdout();

    // We want the cursor to stay at the end of the line, so we print without a newline and flush manually.
    write!(stdout, "Press any key to continue...").unwrap();
    stdout.flush().unwrap();

    // Read a single byte and discard
    let _ = stdin.read(&mut [0u8]).unwrap();
}
