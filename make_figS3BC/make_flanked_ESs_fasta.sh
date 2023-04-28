#!/bin/bash

make_400_flanked_fasta() {
	cat $1 | while read sequence loc
	do
		grep "^>$sequence" $2 -A 1 -m 1 | tail -n +2 > seq.tmp
		roi=''
		start=$(($loc - 400))
		
		if [ $start -le 0 ]
		then
			roi+=$(python -c "print 'N' * $((-$start + 1))")
			start=1
		fi
		
		end=$(($loc + 400))
		roi+=$(cut -c $start-$end seq.tmp)
		
		missing_tail=$(echo $(($end - $(wc -c < seq.tmp))))
		if [ $missing_tail -ge 0 ]
		then
			roi+=$(python -c "print 'N' * $(($missing_tail + 1))")
		fi
		echo ">"$sequence":"$loc
		echo $roi
	done 
}

transcriptome="/Users/matthewbirk/Documents/Career/Bioinformatics_databases/Exomes/Octopus_bimaculoides_Albertin2015_CDS/Octopus_bimaculoides_CDS.fasta"

genome="/Users/matthewbirk/Documents/Career/Bioinformatics_databases/Genomes/Obimaculoides_genome_Albertin_etal.2015/Obimaculoides_genome_Albertin_etal.2015.fasta"







######################### COLD #########################

xlsx2csv ../Data/cold_differential_sites_sorted.xlsx -d "tab" | cut -f 1-3 | tr ":" "\t" | tail -n +2 > cold_sites.txt
cut cold_sites.txt -f 3-4 > tmp
make_400_flanked_fasta tmp $transcriptome > cold_transcriptome_sites_400.fasta

cut cold_sites.txt -f 1-2 > tmp
make_400_flanked_fasta tmp $genome > cold_genome_sites_400.fasta




######################### COLD > 10% #########################

xlsx2csv ../Data/cold_differential_sites_sorted.xlsx -d "tab" | cut -f 1-3,8-9 | tr ":" "\t" | tail -n +2 > tmp

cat tmp | awk '{ $7 = $5 - $6 } 1' | while read scaffold scaffold_loc transcript trans_loc cold warm diff
do
	if [ $(echo "$diff > 0.1" | bc) -eq 1 ]
	then
		echo -e $scaffold"\t"$scaffold_loc"\t"$transcript"\t"$trans_loc
	fi
done > cold10_sites.txt

cut cold10_sites.txt -f 3-4 > tmp
make_400_flanked_fasta tmp $transcriptome > cold10_transcriptome_sites_400.fasta

cut cold10_sites.txt -f 1-2 > tmp
make_400_flanked_fasta tmp $genome > cold10_genome_sites_400.fasta


######################### WARM #########################

xlsx2csv ../Data/warm_differential_sites_sorted.xlsx -d "tab" | cut -f 1-3 | tr ":" "\t" | tail -n +2 > warm_sites.txt
cut warm_sites.txt -f 3-4 > tmp
make_400_flanked_fasta tmp $transcriptome > warm_transcriptome_sites_400.fasta

cut warm_sites.txt -f 1-2 > tmp
make_400_flanked_fasta tmp $genome > warm_genome_sites_400.fasta


######################### INSENSITIVE #########################

cat ../Data/temp_insensitive_sites.txt | tr ":" "\t" | tail -n +2 > insen_sites.txt

cut insen_sites.txt -f 3-4 > tmp
make_400_flanked_fasta tmp $transcriptome > insen_transcriptome_sites_400.fasta

cut insen_sites.txt -f 1-2 > tmp
make_400_flanked_fasta tmp $genome > insen_genome_sites_400.fasta




rm tmp
rm seq.tmp