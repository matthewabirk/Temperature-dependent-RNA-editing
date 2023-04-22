load('/Users/matthewbirk/Documents/WDs/R/Rosenthal/Octo_temp/data_objects.Rdata')

# identify genes in results

setwd('/Users/matthewbirk/Documents/Career/Rosenthal Lab/Temperature Octopus/DAVID/cell membrane/')

tmp = read.table('../cold10_recod_geneids.txt', sep = '\t', quote = '', header = TRUE)

cell_membrane_genes = c(9101, 5144, 23167, 56098, 56132, 91010, 2028, 54510, 11188, 3631, 10893, 5536, 4326, 389072, 22883, 2040, 338567, 27253, 2678, 124930, 1464, 55576, 54769, 526, 2319, 10787, 8396, 6653, 8394, 3416, 23263, 5797, 1636, 9378, 8567, 4641, 2182, 57534, 23390, 1128, 783, 201232, 84455, 55915, 1889, 4771, 4651, 2195, 134957, 55361, 394, 1139, 58513, 55761, 1138, 1136, 1499, 145567, 64072, 5218, 91608, 8727, 23499, 7094, 80036, 287, 56105, 288, 56106, 324, 9429, 57719, 375790, 84079, 5745, 22933, 5361, 65981, 64137, 57282, 28231, 56111, 610, 5359, 10397, 123879, 4163, 8001, 91584, 8642, 9293, 200734, 56121, 26503, 29979, 80727, 80333, 2776, 55714, 747, 5924) # from list of genes in UP_keywords_only.csv KW-1003~Cell membrane

cell_membrane_genes = c(94005, 5144, 8897, 9743, 23167, 56098, 590, 55681, 80762, 2548, 23161, 284339, 201163, 11188, 3631, 10893, 3636, 3092, 389072, 65258, 92293, 22883, 57633, 2678, 55210, 1464, 55576, 5824, 10787, 2319, 4074, 25, 10075, 3416, 85439, 55227, 23263, 23385, 5832, 6262, 8443, 8567, 4641, 9099, 2218, 57534, 1128, 23390, 84455, 23274, 55236, 6259, 55915, 51439, 5066, 5189, 4771, 4651, 7248, 2232, 23001, 84918, 55361, 394, 1139, 58513, 1138, 1136, 1499, 145567, 84364, 11022, 9779, 64072, 10184, 389432, 23499, 287, 56105, 288, 56106, 120892, 11276, 9429, 5745, 22933, 65981, 23347, 23348, 28231, 1954, 56111, 610, 54971, 10397, 8001, 91584, 56121, 26503, 200576, 81544, 80333, 55275, 2137, 23230, 195814, 23352, 54867, 55714, 747, 6857, 26958, 6198, 9101, 3241, 56132, 57221, 91010, 55284, 2820, 2028, 51366, 54510, 4326, 5536, 8925, 7175, 10818, 25998, 2040, 338567, 121665, 27253, 124930, 55852, 90701, 54769, 526, 5789, 9487, 284161, 8396, 158747, 6653, 8394, 80146, 55069, 5797, 83891, 10231, 80142, 91949, 1636, 9378, 9931, 4122, 8720, 25851, 130507, 2182, 80255, 783, 201232, 4919, 85021, 9927, 1889, 25829, 2195, 134957, 51280, 55761, 56850, 10452, 5218, 91608, 8727, 84749, 7094, 80036, 10902, 324, 78999, 57719, 375790, 84079, 9706, 8618, 5361, 23505, 64137, 29966, 9043, 51061, 57282, 23061, 57169, 9955, 5359, 123879, 5250, 4163, 8642, 27315, 9293, 200734, 29979, 80727, 57732, 171586, 2776, 26100, 2531, 54901, 84892, 10564, 5924, 23190) # KW-0472~Membrane

tmp = subset(tmp, To %in% cell_membrane_genes)

library(dplyr)
cell_membrane = filter(cold_recod, `Uhuman uniprot ID` %in% tmp$From, EL_diff > 0.1) %>% select(Transcript, `Location in transcript`, Protein_full_name, `Uhuman uniprot ID`, EL_diff)

cell_membrane$residue = floor(cell_membrane$`Location in transcript` / 3)


write.table(cell_membrane, 'cell_membrane_dEL0.1_ESs.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')

system('cut -f 1 cell_membrane_dEL0.1_ESs.txt | uniq > tmp')
system('grep -f tmp /Users/matthewbirk/Documents/WDs/SnapGene/Exomes/Octopus_bimaculoides_CDS.fasta -A 1 > cell_membrane_transcripts.fasta')

# for some reason need to manually remove --\n from each entry

system('python3 conv_nt2aa.py > cell_membrane_prots.fasta')

# next step is to upload cell_membrane_prots.fasta to https://dtu.biolib.com/DeepTMHMM then save output as TMRs.gff3. Manually remove rows with just //

#TMHMM_results.txt


# used to be http://www.cbs.dtu.dk/services/TMHMM/

# now is v 1.0.13 of DeepTMHMM











tm = read.table('TMRs.gff3', col.names = c('transcript', 'region', 'start', 'end'))

tm_tm = filter(tm, region == 'TMhelix') # only choose proteins that actually have predicted TM helices

length(unique(tm_tm$transcript)) # out of 227 total transcripts, only 78 were predicted to have TM helices

cell_membrane = filter(cell_membrane, Transcript %in% unique(tm_tm$transcript))


tm_tm = filter(tm, region == 'TMhelix') # In what region do you want to look for ES?

cell_membrane$in_TM = apply(cell_membrane, 1, function(i){
		any(apply(filter(tm_tm, transcript == i['Transcript']), 1, function(j){ # among all the TM helices in this protein, does the ES fall in any of them?
			as.numeric(i['residue']) >= as.numeric(j['start']) & as.numeric(i['residue']) <= as.numeric(j['end'])
		}))
})





# This section with tm_meta pulls the metadata from TMHMM to find the total number of residues in all of the TM helix containing proteins.

tm_meta = readLines('TMRs.gff3')
tm_meta = grep('Length:', tm_meta, value = TRUE)
tm_meta = data.frame(transcript = gsub('# (\\w+) Length: (\\d+)', '\\1', tm_meta), 
					 length = as.numeric(gsub('# (\\w+) Length: (\\d+)', '\\2', tm_meta)))
tm_meta = filter(tm_meta, transcript %in% unique(tm_tm$transcript)) # how many residues are in each TM containing protein?




sum(tm_tm$end - tm_tm$start + 1) / sum(tm_meta$length) # 6.8% of residues in TM containing proteins are in TMhelices
table(cell_membrane$in_TM) / nrow(cell_membrane) # 8.4% of ESs are in TMhelices
