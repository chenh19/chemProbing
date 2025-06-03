#!/bin/bash
[ ! -d ./3.analysis/ ] && mkdir ./3.analysis/
[ ! -d ./3.analysis/3.trim/ ] && mkdir ./3.analysis/3.trim/

# install fastqc
sudo apt-get update -qq && sudo apt-get install cutadapt fastp -y


# trim by length and quality
for r1 in ./2.fastq/*_R1*.fastq*; do
    [ -f "$r1" ] || continue

    r2=$(echo "$r1" | sed -E 's/_R1/_R2/')

    if [[ "$r1" == *.fastq.gz ]]; then
        base=$(basename "$r1" .fastq.gz)
        sample=$(echo "$base" | sed -E 's/_R1//')
    elif [[ "$r1" == *.fastq ]]; then
        base=$(basename "$r1" .fastq)
        sample=$(echo "$base" | sed -E 's/_R1//')
    else
        continue
    fi

    cutadapt \
      -j 16 \
      -l 200 \
      -o "./3.analysis/3.trim/${sample}_R1.cut.fastq.gz" \
      -p "./3.analysis/3.trim/${sample}_R2.cut.fastq.gz" \
      $r1 $r2

    fastp \
      -i "./3.analysis/3.trim/${sample}_R1.cut.fastq.gz" \
      -I "./3.analysis/3.trim/${sample}_R2.cut.fastq.gz" \
      -o "./3.analysis/3.trim/${sample}_R1.trimmed.fastq.gz" \
      -O "./3.analysis/3.trim/${sample}_R2.trimmed.fastq.gz" \
      --trim_front1 5 \
      --trim_front2 5 \
      --cut_tail \
      --cut_tail_window_size 4 \
      --cut_tail_mean_quality 20 \
      --length_required 100 \
      --qualified_quality_phred 20 \
      --unqualified_percent_limit 10 \
      --thread 16 \
      -j /dev/null \
      -h /dev/null

      rm -f "./3.analysis/3.trim/${sample}_R1.cut.fastq.gz" "./3.analysis/3.trim/${sample}_R2.cut.fastq.gz"
done
