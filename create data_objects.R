
# all sequences in exome must start with the start of the first codon
find_recode = function(df, transcript_col, es_col, exome){
	apply(df, 1, function(i){
		es = as.numeric(i[es_col])
		codon = switch(as.character(es %% 3),
									 '0' = (es-2):es,
									 '1' = es:(es+2),
									 '2' = (es-1):(es+1),
		)
		seq = tolower(strsplit(exome[i[transcript_col]][[1]], '')[[1]][codon])
		seq[es == codon] = 'a'
		orig = seqinr::translate(seq)
		seq[es == codon] = 'g'
		edit = seqinr::translate(seq)
		change = ifelse(orig == edit, 'syn', paste0(orig, edit, collapse = ''))
		return(change)
	})
}




library(readxl)
library(dplyr)


Obimac_CDS = seqinr::read.fasta('/Users/matthewbirk/Documents/WDs/SnapGene/Exomes/Octopus_bimaculoides_CDS.fasta', seqtype = 'DNA', as.string = TRUE) # data from http://octopus.unit.oist.jp/OCTDATA/BASIC/Metazome/Obimaculoides_280_cds.fa.gz


all_edits = as.data.frame(read_xlsx('/Users/matthewbirk/Documents/WDs/R/Rosenthal/Octo_temp/Data/SGTemprature_CDS_100reads.xlsx'))
all_edits$transcript = gsub('(\\w+):\\d+', '\\1', all_edits$Transcript)
all_edits$es = as.numeric(gsub('\\w+:(\\d+)', '\\1', all_edits$Transcript))

all_edits$`Codon changes` = find_recode(df = all_edits, transcript_col = 'transcript', es_col = 'es', exome = Obimac_CDS)


all_recod_edits = subset(all_edits, all_edits$`Codon changes` != 'syn')

cold = as.data.frame(read_xlsx('/Users/matthewbirk/Documents/WDs/R/Rosenthal/Octo_temp/Data/cold_differential_sites_sorted.xlsx'))
warm = as.data.frame(read_xlsx('/Users/matthewbirk/Documents/WDs/R/Rosenthal/Octo_temp/Data/warm_differential_sites_sorted.xlsx'))

cold$ES = paste(cold$Transcript, cold$`Location in transcript`, sep = '_')
warm$ES = paste(warm$Transcript, warm$`Location in transcript`, sep = '_')



all_edits$neighbors10_CDS = unlist(apply(all_edits, 1, function(i){
	inter1 = as.character(Obimac_CDS[[i['transcript']]])
	if(length(inter1) == 0) return(NA)
	inter2 = as.numeric(i['es'])
	stringr::str_sub(inter1, inter2 - 10, inter2 + 10)
}))

all_edits$neighbors50_CDS = unlist(apply(all_edits, 1, function(i){
	inter1 = as.character(Obimac_CDS[[i['transcript']]])
	if(length(inter1) == 0) return(NA)
	inter2 = as.numeric(i['es'])
	stringr::str_sub(inter1, inter2 - 50, inter2 + 50)
}))

all_edits$gc10_CDS = apply(all_edits, 1, function(i){
	inter1 = as.character(Obimac_CDS[[i['transcript']]])
	inter2 = as.numeric(i['es'])
	inter3 = stringr::str_sub(inter1, inter2 - ifelse(inter2 < 10, inter2, 10), inter2 + 10)
	seqinr::GC(seqinr::s2c(inter3))
})


cold$EL_cold = rowMeans(select(cold, starts_with('Editing cold')))
cold$EL_warm = rowMeans(select(cold, starts_with('Editing warm')))
cold$EL_inc = cold$EL_cold / cold$EL_warm
cold$EL_diff = cold$EL_cold - cold$EL_warm
cold$neighbors = apply(cold, 1, function(i){
	inter1 = as.character(Obimac_CDS[[i['Transcript']]])
	inter2 = as.numeric(i['Location in transcript'])
	stringr::str_sub(inter1, inter2 - 1, inter2 + 1)
})
cold$neighbors10 = apply(cold, 1, function(i){
	inter1 = as.character(Obimac_CDS[[i['Transcript']]])
	inter2 = as.numeric(i['Location in transcript'])
	stringr::str_sub(inter1, inter2 - 10, inter2 + 10)
})
# cold$neighbors50 = apply(cold, 1, function(i){
# 	inter1 = as.character(Obimac_CDS[[i['Transcript']]])
# 	inter2 = as.numeric(i['Location in transcript'])
# 	stringr::str_sub(inter1, inter2 - 50, inter2 + 50)
# })
cold$gc10 = apply(cold, 1, function(i){
	inter1 = as.character(Obimac_CDS[[i['Transcript']]])
	inter2 = as.numeric(i['Location in transcript'])
	inter3 = stringr::str_sub(inter1, inter2 - ifelse(inter2 < 10, inter2, 10), inter2 + 10)
	seqinr::GC(seqinr::s2c(inter3))
})

warm$EL_cold = rowMeans(select(warm, starts_with('Editing cold')))
warm$EL_warm = rowMeans(select(warm, starts_with('Editing warm')))
warm$EL_inc = warm$EL_warm / warm$EL_cold
warm$EL_diff = warm$EL_warm - warm$EL_cold
warm$neighbors = apply(warm, 1, function(i){
	inter1 = as.character(Obimac_CDS[[i['Transcript']]])
	inter2 = as.numeric(i['Location in transcript'])
	stringr::str_sub(inter1, inter2 - 1, inter2 + 1)
})
warm$neighbors10 = apply(warm, 1, function(i){
	inter1 = as.character(Obimac_CDS[[i['Transcript']]])
	inter2 = as.numeric(i['Location in transcript'])
	stringr::str_sub(inter1, inter2 - 10, inter2 + 10)
})
# warm$gc20 = apply(warm, 1, function(i){
# 	inter1 = as.character(Obimac_CDS[[i['Transcript']]])
# 	inter2 = as.numeric(i['Location in transcript'])
# 	inter3 = stringr::str_sub(inter1, inter2 - ifelse(inter2 < 20, inter2, 20), inter2 + 20)
# 	seqinr::GC(seqinr::s2c(inter3))
# })

warm$gc10 = apply(warm, 1, function(i){
	inter1 = as.character(Obimac_CDS[[i['Transcript']]])
	inter2 = as.numeric(i['Location in transcript'])
	inter3 = stringr::str_sub(inter1, inter2 - ifelse(inter2 < 10, inter2, 10), inter2 + 10)
	seqinr::GC(seqinr::s2c(inter3))
})

all_temp = rbind(cold, warm)

all_temp$dEL = (all_temp$EL_warm - all_temp$EL_cold) * 100







warm_samples = grep('Editing warm\\d', colnames(all_temp), value = TRUE)
cold_samples = grep('Editing cold\\d', colnames(all_temp), value = TRUE)

all_temp$warm_cv = apply(all_temp, 1, function(i){
	tmp = as.numeric(i[warm_samples])
	sd(tmp, na.rm = TRUE) / mean(tmp, na.rm = TRUE)
})
all_temp$cold_cv = apply(all_temp, 1, function(i){
	tmp = as.numeric(i[cold_samples])
	sd(tmp, na.rm = TRUE) / mean(tmp, na.rm = TRUE)
})

all_temp$warm_range = apply(all_temp, 1, function(i) diff(range(as.numeric(i[warm_samples]), na.rm = TRUE)))
all_temp$cold_range = apply(all_temp, 1, function(i) diff(range(as.numeric(i[cold_samples]), na.rm = TRUE)))

all_temp$warm_range_high_cov = apply(all_temp, 1, function(i) diff(range(as.numeric(i[warm_samples[1:3]]), na.rm = TRUE)))
all_temp$cold_range_high_cov = apply(all_temp, 1, function(i) diff(range(as.numeric(i[cold_samples[1:3]]), na.rm = TRUE)))

all_temp$warm_range_low_cov = apply(all_temp, 1, function(i) diff(range(as.numeric(i[warm_samples[4:6]]), na.rm = TRUE)))
all_temp$cold_range_low_cov = apply(all_temp, 1, function(i) diff(range(as.numeric(i[cold_samples[4:6]]), na.rm = TRUE)))

temp_insensitive = anti_join(all_edits, all_temp, by = c('transcript' = 'Transcript', 'es' = 'Location in transcript'))


hist(all_temp$EL_diff, main = 'Octopus whole animal', xlab = 'Cold - Warm', xlim = c(-1, 1), breaks = seq(-1, 1, by = 0.05))
abline(v = 0, col = 'red')





warm$change = find_recode(df = warm, transcript_col = 'Transcript', es_col = 'Location in transcript', exome = Obimac_CDS)
cold$change = find_recode(df = cold, transcript_col = 'Transcript', es_col = 'Location in transcript', exome = Obimac_CDS)


warm_recod = warm[warm$change != 'syn', ]
warm_syn = warm[warm$change == 'syn', ]

cold_recod = cold[cold$change != 'syn', ]
cold_syn = cold[cold$change == 'syn', ]

all_temp_recod = rbind(warm_recod, cold_recod)




# inter1 = subset(all_temp_recod, EL_diff >= 0.1 & all_temp_recod$`Uhuman uniprot ID` != 'NA')
# library(UniProt.ws)
# human = UniProt.ws()
# inter2 = select(x = human, keys = inter1$`Uhuman uniprot ID`, columns = 'GO', keytype = 'UNIPROTKB')
# inter3 = paste(inter2$GO, collapse = '; ')
# inter4 = strsplit(inter3, split = '; ')[[1]]

# subset(cold_recod, cold_recod$`Uhuman uniprot ID` != 'NA')$`Uhuman uniprot ID`
# 
# dplyr::arrange(subset(all_temp_recod, Protein_full_name != 'NA' & EL_diff > 0.25), desc(EL_diff))
#dplyr::arrange(subset(all_temp_recod, Protein_full_name != 'NA' & `Uhuman uniprot ID` %in% conserved$uniprot_name & EL_diff > 0.2), desc(EL_diff))

save.image('/Users/matthewbirk/Documents/WDs/R/Rosenthal/Octo_temp/data_objects.Rdata')
