/// Helper functions
use super::conf::{SmallVarSpec, StrucVarSpec, SvType, VarSpec};

/// Return whether two intervals overlap.
fn coords_overlap(x1: i64, x2: i64, y1: i64, y2: i64) -> bool {
    x1 < y2 && y1 < x2
}

pub fn small_var_spec_overlaps(
    var_spec: &SmallVarSpec,
    itv_contig: &str,
    itv_begin: i64,
    itv_end: i64,
) -> bool {
    itv_contig.eq(&var_spec.chromosome)
        && coords_overlap(itv_begin, itv_end, var_spec.start - 1, var_spec.end)
}

pub fn struc_var_spec_overlaps(
    var_spec: &StrucVarSpec,
    itv_contig: &str,
    itv_begin: i64,
    itv_end: i64,
) -> bool {
    match var_spec.sv_type {
        SvType::Deletion => {
            itv_contig.eq(&var_spec.chromosome)
                && coords_overlap(itv_begin, itv_end, var_spec.start - 1, var_spec.end)
        }
    }
}

/// Return wether the given variant overlaps with the given chromosomal interval (0-based!)
pub fn spec_overlaps(var_spec: &VarSpec, itv_contig: &str, itv_begin: i64, itv_end: i64) -> bool {
    match var_spec {
        VarSpec::SmallVar(small_var_spec) => {
            small_var_spec_overlaps(&small_var_spec, itv_contig, itv_begin, itv_end)
        }
        VarSpec::StrucVar(struc_var_spec) => {
            struc_var_spec_overlaps(&struc_var_spec, itv_contig, itv_begin, itv_end)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::conf::{SmallVarSpec, StrucVarSpec, SvType, VarSpec};
    use super::{coords_overlap, spec_overlaps};

    #[test]
    fn test_coords_overlap() {
        assert!(!coords_overlap(0, 100, 100, 200));
        assert!(!coords_overlap(100, 200, 0, 100));
        assert!(coords_overlap(0, 101, 100, 200));
        assert!(coords_overlap(0, 200, 150, 160));
    }

    #[test]
    fn test_spec_overlaps_small_var_snv() {
        let var = VarSpec::SmallVar(SmallVarSpec {
            chromosome: "chr5".to_string(),
            start: 100,
            end: 100,
            reference: "A".to_string(),
            alternative: "T".to_string(),
            aaf: 0.5,
        });

        assert!(!spec_overlaps(&var, "chr5", 98, 99));
        assert!(spec_overlaps(&var, "chr5", 99, 100));
        assert!(!spec_overlaps(&var, "chr5", 100, 101));

        assert!(!spec_overlaps(&var, "chr4", 99, 100));
    }

    #[test]
    fn test_spec_overlaps_small_var_indel() {
        let var = VarSpec::SmallVar(SmallVarSpec {
            chromosome: "chr5".to_string(),
            start: 100,
            end: 101,
            reference: "AT".to_string(),
            alternative: "A".to_string(),
            aaf: 0.5,
        });

        assert!(!spec_overlaps(&var, "chr5", 98, 99));
        assert!(spec_overlaps(&var, "chr5", 99, 100));
        assert!(spec_overlaps(&var, "chr5", 100, 101));
        assert!(!spec_overlaps(&var, "chr5", 102, 102));

        assert!(!spec_overlaps(&var, "chr4", 99, 101));
    }

    #[test]
    fn test_spec_overlaps_struc_var_deletion() {
        let var = VarSpec::StrucVar(StrucVarSpec {
            sv_type: SvType::Deletion,
            chromosome: "chr5".to_string(),
            start: 1000,
            end: 2000,
            aaf: 0.5,
        });

        assert!(!spec_overlaps(&var, "chr5", 900, 999));
        assert!(spec_overlaps(&var, "chr5", 900, 1000));
        assert!(spec_overlaps(&var, "chr5", 1999, 2100));
        assert!(!spec_overlaps(&var, "chr5", 2000, 2100));

        assert!(!spec_overlaps(&var, "chr4", 999, 2000));
    }
}
