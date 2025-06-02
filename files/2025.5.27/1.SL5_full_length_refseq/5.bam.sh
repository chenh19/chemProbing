#!/bin/bash
[ ! -d ./3.analysis/ ] && mkdir ./3.analysis/
[ ! -d ./3.analysis/5.bam/ ] && mkdir ./3.analysis/5.bam/

# install minimap2 samtools
sudo apt-get update -qq && sudo apt-get install minimap2 samtools -y

# align reads
for f in ./3.analysis/3.trim/*.fastq.gz; do
  base=$(basename $f .fastq.gz)
  minimap2 -t 16 -ax map-ont ./3.analysis/1.refseq/refseq.fa $f | \
    samtools sort -@ 16 -o ./3.analysis/5.bam/${base}.bam
  samtools index ./3.analysis/5.bam/${base}.bam
done
