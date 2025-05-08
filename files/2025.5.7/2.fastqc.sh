#!/bin/bash
[ ! -d ./3.analysis/ ] && mkdir ./3.analysis/

# install fastqc
sudo apt-get update -qq && sudo apt-get install fastqc -y

# run fastqc in parallel
[ ! -d ./3.analysis/2.fastqc/ ] && mkdir ./3.analysis/2.fastqc/
fastqc --threads 16 ./2.fastq/*.fastq --outdir=./3.analysis/2.fastqc/
