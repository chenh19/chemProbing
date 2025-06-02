#!/bin/bash
[ ! -d ./3.analysis/ ] && mkdir ./3.analysis/
[ ! -d ./3.analysis/3.trim/ ] && mkdir ./3.analysis/3.trim/

# install fastqc
sudo apt-get update -qq && sudo apt-get install fastp seqtk -y

# trim by length
for f in ./2.fastq/*.fastq ./2.fastq/*.fastq.gz; do
    [ -e "$f" ] || continue
    if [[ "$f" == *.gz ]]; then
        base=$(basename "$f" .fastq.gz)
        seqtk trimfq -L 400 "$f" | gzip > "./3.analysis/3.trim/${base}_cut.fastq.gz"
    else
        base=$(basename "$f" .fastq)
        seqtk trimfq -L 400 "$f" | gzip > "./3.analysis/3.trim/${base}_cut.fastq.gz"
    fi
done

# trim by quality
for f in ./3.analysis/3.trim/*.fastq.gz; do
    base=$(basename "$f" .fastq.gz)
    fastp \
      -i "$f" \
      -o "./3.analysis/3.trim/${base}_trimmed.fastq.gz" \
      --cut_tail \
      --cut_tail_window_size 4 \
      --cut_tail_mean_quality 20 \
      --length_required 200 \
      --qualified_quality_phred 20 \
      --unqualified_percent_limit 20 \
      --thread 16 \
      -j /dev/null \
      -h /dev/null
    rm -f "$f"
done
