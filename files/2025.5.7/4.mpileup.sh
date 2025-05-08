#!/bin/bash
[ ! -d ./3.analysis/ ] && mkdir ./3.analysis/

# install samtools
sudo apt-get update -qq && sudo apt-get install samtools -y

# pileup
[ ! -d ./3.analysis/4.mpileup/ ] && mkdir ./3.analysis/4.mpileup/
for bam in ./3.analysis/3.bam/*.bam; do
  base=$(basename "$bam" .bam)
  samtools mpileup -aa -A -B -Q 0 -d 1000000 -f ./3.analysis/1.refseq/refseq.fa "$bam" > ./3.analysis/4.mpileup/"$base".mpileup
done
