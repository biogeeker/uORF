#!/bin/bash

cat $1 | grep hg19  > tmp1.txt

cat tmp1.txt | cut -d ' ' -f 1 | cut -d '_' -f 3- > transcript_name.txt

cat tmp1.txt | cut -d ' ' -f 2 | cut -d '=' -f 2 > tmp2.txt

cat tmp2.txt | cut -d ':' -f 1 > chr.txt
cat tmp2.txt | cut -d ':' -f 2 | cut -d '-' -f 1 > start.txt
cat tmp2.txt | cut -d ':' -f 2 | cut -d '-' -f 2 > end.txt

cat tmp1.txt | cut -d ' ' -f 5 | cut -d '=' -f 2 > strand.txt

/rsrch2/bcb/mxu3/scratch/tools/fastx/bin/fasta_formatter -i $1 -t | cut -f 2 > sequence.txt

awk -f vlookup.awk RefSeq_to_symbol.txt transcript_name.txt | cut -d ' ' -f 2 > symbol.txt

paste transcript_name.txt symbol.txt chr.txt start.txt end.txt strand.txt sequence.txt > parsed.txt

rm transcript_name.txt chr.txt tmp1.txt tmp2.txt start.txt end.txt strand.txt sequence.txt symbol.txt

