load('/Users/matthewbirk/Documents/WDs/R/Rosenthal/Octo_temp/data_objects.Rdata')

cold10_uniprots = subset(cold_recod, EL_diff > 0.1)$`Uhuman uniprot ID`
cold10_uniprots = cold10_uniprots[which(cold10_uniprots != 'NA')]

write.table(cold10_uniprots, '/Users/matthewbirk/Documents/Career/Rosenthal Lab/Temperature Octopus/DAVID/cold10_uniprot_ids.txt', quote = FALSE, row.names = FALSE, col.names = FALSE) 



cold5_uniprots = subset(cold_recod, EL_diff > 0.05)$`Uhuman uniprot ID`
cold5_uniprots = unique(cold5_uniprots[which(cold5_uniprots != 'NA')])

cold5_uniprots




# Get annotation for temp. insensitive sites -----------------------------

map = readxl::read_xlsx('/Users/matthewbirk/Documents/Career/Rosenthal Lab/Liscovitch-Brauer2017 data/ES, bimac genome - Liscovitch-Brauer2017, Table S5. Editing Sites Detected (Genome Alignment, Oct. Bim.).xlsx', sheet = 'Overlapping_locations_with_REDI')

all_recod_edits_annotated = left_join(all_recod_edits, unique(map[, c('Trinity name', 'Oct.bim. transcript name')]), by = c('transcript' = 'Oct.bim. transcript name'))

my_map = readxl::read_xlsx('/Users/matthewbirk/Documents/Career/Rosenthal Lab/Liscovitch-Brauer2017 data/Ob ESs - ORF vs CDS.xlsx')

all_recod_edits_annotated = left_join(all_recod_edits_annotated, unique(my_map[, c('Ob_transcriptome_id', 'Ob_exome_id')]), by = c('transcript' = 'Ob_exome_id'))

table(all_recod_edits_annotated$`Trinity name` == all_recod_edits_annotated$Ob_transcriptome_id) # my mapping agrees with Table S5 pretty well

table(is.na(all_recod_edits_annotated$`Trinity name`))
table(is.na(all_recod_edits_annotated$Ob_transcriptome_id))
table(is.na(all_recod_edits_annotated$Ob_transcriptome_id) & is.na(all_recod_edits_annotated$`Trinity name`)) # combining both annotations helps fill in a few other gaps

all_recod_edits_annotated[is.na(all_recod_edits_annotated$`Trinity name`), 'Trinity name'] = all_recod_edits_annotated[is.na(all_recod_edits_annotated$`Trinity name`), 'Ob_transcriptome_id'] # for any transcripts where Table S5 was missing the match, insert my own mapping to fill in gaps.



write.table(na.omit(unique(all_recod_edits_annotated$`Trinity name`)), '/Users/matthewbirk/Documents/Career/Rosenthal Lab/Temperature Octopus/DAVID/all_recoding_trinity_ids.txt', quote = FALSE, row.names = FALSE, col.names = FALSE) 


system('bash "/Users/matthewbirk/Documents/Career/Rosenthal Lab/Temperature Octopus/DAVID/convert_background_trinity_to_uniprot.sh"') # use these in DAVID as background
