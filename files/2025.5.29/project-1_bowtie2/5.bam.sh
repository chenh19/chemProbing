#!/bin/bash
[ ! -d ./3.analysis/ ] && mkdir ./3.analysis/
[ ! -d ./3.analysis/5.bam/ ] && mkdir ./3.analysis/5.bam/

# install minimap2 samtools
sudo apt-get update -qq && sudo apt-get install bowtie2 samtools -y

# align reads
bowtie2-build ./3.analysis/1.refseq/refseq.fa ./3.analysis/1.refseq/refseq
for r1 in ./3.analysis/3.trim/*_R1*trimmed.fastq.gz; do
  [ -f "$r1" ] || continue
  r2=$(echo "$r1" | sed -E 's/_R1/_R2/')
  base=$(basename "$r1" .trimmed.fastq.gz)
  sample=$(echo "$base" | sed -E 's/_R1//')
  bowtie2 -x ./3.analysis/1.refseq/refseq \
          -1 "$r1" -2 "$r2" \
          -p 16 | \
    samtools sort -@ 16 -o ./3.analysis/5.bam/${sample}.bam
  samtools index ./3.analysis/5.bam/${sample}.bam
done
