#! /bin/bash

dir=$1

list=$(ls $dir | grep _alleles.txt)

> all_results.txt

for i in $list; do
	name=$(echo $i | sed 's/_R_alleles.txt//')	
	echo $name
	
	### create list of all genes
	#grep -v " \* " $i | cut -f 2 | sed 's/*.*//' | sort | uniq -c > ${name}_genes.txt
	#genes=$(grep -v " \* " $i | cut -f 2 | sed 's/*.*//' | sort | uniq)
	genes=$(cut -f 2 $i | sed 's/*.*//' | sort | uniq | grep -v allel)
		
	> ${name}_alleles_final.txt
	> ${name}_final_result.txt
		
	for g in $genes; do
		echo >> ${name}_alleles_final.txt
		echo $g | sed 's/"//g' | sed 's/*/./g' >> ${name}_alleles_final.txt
		reads=$(grep "${g}\*" $i | cut -f 1)
		#n=length($reads)
		
		### get frequency drop   
		j=$(python /data/AbX/germline/GermAb/freq_drop.py $reads)
		
		###
		grep "${g}\*" $i | head -$j | cut -f 1,2,3 | sed 's/"//g' | sed 's/*/./g' >> ${name}_alleles_final.txt
		grep "${g}\*" $i | head -$j | cut -f 2,3 | uniq | sed 's/"//g' | sed 's/*/./g' >> ${name}_final_result.txt		

		echo $name $g $j | sed 's/"//g' >> all_results.txt

		#if [ $n > 1 ]
		#then
		echo -------------------- >> ${name}_alleles_final.txt
		((j++))
		grep "${g}\*" $i | tail -n +$j | cut -f 1,2,3 | sed 's/"//g' | sed 's/*/./g' >> ${name}_alleles_final.txt
		#fi
		
	done
done
