#!/bin/bash

sp_fasta="/Users/matthewbirk/Documents/Career/Bioinformatics_databases/uniprot_sprot.fasta"

############################## GET BIOGRID SEQS ##############################

cat sources/biogrid_Hs_adar1.txt | cut -f 12 | sort | uniq -c
cat sources/biogrid_Hs_adar2.txt | cut -f 12 | sort | uniq -c

# focusing on only in vivo evidence based on descriptions here: https://wiki.thebiogrid.org/doku.php/experimental_systems

cat $(readlink -f sources/* | grep "biogrid") | cut -f 1,8,9,12,14,15,24,27,36,37 | grep -e "Affinity Capture-MS" -e "Affinity Capture-RNA" -e "Affinity Capture-Western" -e "Co-localization" > good_biogrid_data.txt

cat good_biogrid_data.txt | cut -f 7-8 | tr '\t' '\n' | sort | uniq | grep -v -e "SWISS-PROT" -e "-" > biogrid_hits.txt

echo "After filtering for just in vivo evidence, there are $(wc -l < biogrid_hits.txt) human proteins known to interact with ADAR1 or ADAR2."

cat biogrid_hits.txt | while read i
do
	grep "^>sp|$i" $sp_fasta -A 1 -m 1
done > SOIs.fasta



############################## GET FREUND SEQS ##############################

cat sources/freund_etal_2020.txt | cut -f 1 | sort | uniq | while read i
do
	grep "GN=${i}" $sp_fasta -A 1 -m 1
done > freund.fasta

cat freund.fasta >> SOIs.fasta

grep ">" freund.fasta | cut -c 5-10 > freund_sp_IDs.txt

grep "GN=.*" -o freund.fasta | cut -f 1 -d " " | sed 's/GN=//' > tmp

paste freund_sp_IDs.txt tmp > tmp2

mv tmp2 freund_sp_IDs.txt



############################## GET STRING SEQS ##############################

#cat $(readlink -f sources/* | grep ".fa") >> SOIs.fasta	# add String sequences









tblastn -db "/Users/matthewbirk/Documents/Career/Bioinformatics_databases/Exomes/Octopus_bimaculoides_Albertin2015_CDS/Octopus_bimaculoides_CDS.db" -query SOIs.fasta -outfmt "6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs 1 -max_hsps 1 -evalue 1e-5 > tmp

cat tmp | awk '{ $14 = $5 / $4 } 1' OFS="\t" > blast_results.txt

awk -F"\t" '$14>0.5' blast_results.txt > tmp

echo "I selected the best blast hit to the octopus exome for each protein that had an e-value < 1e-5 and a hit covering > 50% of the query's length. $(wc -l < tmp) proteins passed this filter."

mv tmp blast_results.txt

cut -f 2 blast_results.txt | sort | uniq > Ob_SOIs.txt

echo "Some of those genes blasted to the same octopus sequence, so the final count of octopus genes is $(wc -l < Ob_SOIs.txt)."







#Rscript DE_analysis.R # I like to run this manually