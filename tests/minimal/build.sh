#!/usr/bin/bash

samtools faidx ../../../hs38DH.fa chr1:10,000,000-10,010,000 > ref.fa
perl -p -i -e 's/^>.*/>chrom/g' ref.fa
samtools faidx ref.fa
wgsim -N $((10000 * 30 / 300)) -r 0.0 -1 150 -2 150 -d 400 ref.fa reads_1.fq reads_2.fq
gzip reads_1.fq
gzip reads_2.fq
bwa index ref.fa
bwa mem ref.fa reads_1.fq.gz reads_2.fq.gz \
| samtools sort -O BAM -o reads.bam
samtools index reads.bam
