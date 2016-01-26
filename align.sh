#! /bin/bash

dir=$1

list=$(ls $dir | grep _uniq.fasta)

alleles="/data/AbX/germline/GermAb/IGHV_alleles.fasta"

bwa index $alleles

for i in $list; do
	name=$(echo $i | sed 's/_uniq.fasta//')
	echo Analysing sample $name
	bwa mem $alleles $i > ${name}_aligned.sam
	cut -f 1,3,4,10,12,13 ${name}_aligned.sam | sed '/^@/d' > ${name}_aligned_comb.txt
done
