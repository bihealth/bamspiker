# espike
Spike variants into BAM files (without all the hassle)

## Running SNV Example

```
$ rm -f out-reads.bam*
$ cargo run -- \
    tests/minimal/config-snv.yaml \
    tests/minimal/ref.fa \
    tests/minimal/reads.bam \
    out-reads.bam
$ samtools index out-reads.bam
$ samtools tview -p chrom:1000 out-reads.bam tests/minimal/ref.fa
```

## Running Deletion Example

```
$ rm -f out-reads.bam*
$ cargo run -- \
    tests/minimal/config-del.yaml \
    tests/minimal/ref.fa \
    tests/minimal/reads.bam \
    out-reads.bam
$ samtools index out-reads.bam
$ samtools tview -p chrom:1000 out-reads.bam tests/minimal/ref.fa
```

## Limitations

- SNVs and deletion SVs only at the moment
- Variants must not overlap (no sanity check at the moment)
