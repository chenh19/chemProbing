#!/bin/bash
[ ! -d ./3.analysis/ ] && mkdir ./3.analysis/

# install minimap2 samtools
sudo apt-get update -qq && sudo apt-get install minimap2 samtools -y

# align reads
[ ! -d ./3.analysis/3.bam/ ] && mkdir ./3.analysis/3.bam/
for sample in ./2.fastq/*.fastq; do
  base=$(basename $sample .fastq)
  minimap2 -t 16 -ax map-ont ./3.analysis/1.refseq/refseq.fa $sample | \
    samtools sort -@ 16 -o ./3.analysis/3.bam/${base}.bam
  samtools index ./3.analysis/3.bam/${base}.bam
done
