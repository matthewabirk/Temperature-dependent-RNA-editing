#!/bin/bash

Rscript top10.R

cat top10.txt | while read i
do
	grep $i "/Users/matthewbirk/Documents/Career/Bioinformatics_databases/Exomes/Octopus_bimaculoides_Albertin2015_CDS/Octopus_bimaculoides_CDS.fasta" -A 1
done > top10.fasta

faTrans top10.fasta top10_prot.fasta

blastp -db nr -remote -query top10_prot.fasta -entrez_query "Mollusca[organism]" -max_target_seqs 50 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle' -out blast_results.txt