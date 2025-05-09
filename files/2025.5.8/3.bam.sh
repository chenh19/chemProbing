#!/bin/bash
[ ! -d ./3.analysis/ ] && mkdir ./3.analysis/

# install minimap2 samtools
sudo apt-get update -qq && sudo apt-get install bwa samtools -y

# align reads
[ ! -d ./3.analysis/3.bam/ ] && mkdir ./3.analysis/3.bam/
bwa index ./3.analysis/1.refseq/refseq.fa
for r1 in ./2.fastq/*_R1*.fastq; do
  base=$(basename $r1)
  r2=./2.fastq/${base/_R1_/_R2_}
  sample=$(echo "$base" | sed 's/_R1//g')
  sample=$(echo "$sample" | sed 's/.fastq//g')
  bwa mem -t 16 ./3.analysis/1.refseq/refseq.fa $r1 $r2 | \
    samtools sort -@ 16 -o ./3.analysis/3.bam/${sample}.bam
  samtools index ./3.analysis/3.bam/${sample}.bam
done
