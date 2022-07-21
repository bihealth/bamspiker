# bamspiker

Spike variants into BAM files (with limited functionality but without any hassle).

## Dependencies:

- `samtools` must be in your path

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

## Dependencies

```
$ mamba create -y -n bamspiker cmake openssl
```
