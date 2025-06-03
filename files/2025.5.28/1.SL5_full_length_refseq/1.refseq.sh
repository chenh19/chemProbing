#!/bin/bash
[ ! -d ./3.analysis/ ] && mkdir ./3.analysis/
[ ! -d ./3.analysis/1.refseq/ ] && mkdir ./3.analysis/1.refseq/

# combine refseq
cat ./1.ref/*.fa > ./3.analysis/1.refseq/refseq.fa
