#! /bin/bash
# merge reads with pandaseq, minimal overlap
# trim random nucleotides/primers at beginning and end of read witn primer_trim.py
# collapes unique reads, keep only if more than n members
# output is ${name}_panda.fasta, ${name}_trimmed.fasta, ${name}_uniq.fasta

dir=$1

gunzip ${dir}/*.gz

### create list of unzipped fastq files
list=$(ls $dir | grep _L001_R1_001.fastq | grep -v .gz | grep -v Undetermined)

for i in $list; do
	R1=$i
	R2=$(echo $R1 | sed 's/_R1_001.fastq/_R2_001.fastq/')
	name=$(basename $R1 | sed 's/_L001_R1_001.fastq//')

	### merge reads with pandaseq
	echo merging $R1 $R2 $name
	minimal=10
	pandaseq -o $minimal -f ${dir}/${R1} -r ${dir}/${R2} -w ${name}_panda.fasta
		
	### trim random nucleotides and primer with primer_trim.py
	fwd_primer="/data/AbX/germline/GermAb/germlineprimer_fw.fasta"
	rev_primer="/data/AbX/germline/GermAb/germlineprimer_rv.fasta"
	echo trimming $name
	/data/AbX/germline/GermAb/primer_trim.py ${name}_panda.fasta $fwd_primer $rev_primer > ${name}_trimmed.fasta
	
	### collapes unique reads, keep only if more than n members
	echo collapsing $name
	grep -v M0 ${name}_trimmed.fasta | sort | uniq -c | sort -nr | awk '$1>=10 {print ">"$1"\n"$2 }' > ${name}_uniq.fasta

done

gzip ${dir}/*.fastq