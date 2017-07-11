#! /bin/bash
#aligns merged, trimmed, collapsed reads to IMGT reference, 
#output: ${name}_aligned.sam, ${name}_aligned.txt (txt file contains readcount, allele, start of alignment, sequence, mutations, cigar)

dir=$1

list=$(ls $dir | grep _uniq.fasta)

alleles="/data/AbX/germline/GermAb/IGHV_alleles_IMGT.fasta"

bwa index $alleles

for i in $list; do
	name=$(echo $i | sed 's/_uniq.fasta//')
	echo Analysing sample $name
	bwa mem -L7,7 $alleles $i > ${name}_aligned.sam
	cut -f 1,3,4,10,12,13,16 ${name}_aligned.sam | sed '/^@/d' > ${name}_aligned.txt
done
