#!/bin/sh

rm exon_start.txt
rm exon_end.txt
rm exon_boundary.txt

for i in `cat $1`; 

do 

	grep -P $i['\t'_] RefSeq_5UTR_exon.fasta.parsed.txt > gene_specific_exons.txt

	a=`cat gene_specific_exons.txt | wc -l` 

	if [ "$a" -gt 0 ]; then 	
		cat gene_specific_exons.txt | cut -f 3 | tr '\n' ',' | sed 's/,$/\n/'	>> exon_start.txt
		cat gene_specific_exons.txt | cut -f 4 | tr '\n' ',' | sed 's/,$/\n/'   >> exon_end.txt
	else
		echo "NA" >> exon_start.txt
		echo "NA" >> exon_end.txt
	fi

	rm gene_specific_exons.txt
done

paste exon_start.txt exon_end.txt > exon_boundary.txt
