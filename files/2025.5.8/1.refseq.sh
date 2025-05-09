#!/bin/bash
[ ! -d ./3.analysis/ ] && mkdir ./3.analysis/

# combine refseq
[ ! -d ./3.analysis/1.refseq/ ] && mkdir ./3.analysis/1.refseq/
cat ./1.ref/*.fa > ./3.analysis/1.refseq/refseq.fa
