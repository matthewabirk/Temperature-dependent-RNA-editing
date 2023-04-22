# Note that the two CSV files in the results directory will need to be gunzipped before running this code will work.

library(dplyr)

load('../data_objects.Rdata')
all_edits$dEL = all_edits$`mean cold` - all_edits$`mean warm`
colnames(all_edits) = gsub('mean ', 'mean_', colnames(all_edits))



genome = read.csv('results/genome_all_sites_with_folding_results.csv')
colnames(genome) = gsub('stracture', 'structure', colnames(genome))
genome$temp = ordered(gsub('(\\w+)_genome_sites_400.fasta', '\\1', genome$file), levels = c('cold10', 'cold', 'insen', 'warm'))
genome$dG_full = genome$energy_full_structure_temp286.15.K. - genome$energy_full_structure_temp295.15.K.
genome$dG_local = genome$dsRNA.ES._energy_temp286.15.K. - genome$dsRNA.ES._energy_temp295.15.K.

genome = dplyr::arrange(genome, temp)
genome = genome[-which(duplicated(genome$ID)), ]

genome = subset(genome, dsRNA.ES._energy_temp295.15.K. != 0 & dsRNA.ES._energy_temp286.15.K. != 0)











transcriptome = read.csv('results/transcriptome_all_sites_with_folding_results.csv')
colnames(transcriptome) = gsub('stracture', 'structure', colnames(transcriptome))
transcriptome$temp = ordered(gsub('(\\w+)_transcriptome_sites_400.fasta', '\\1', transcriptome$file), levels = c('cold10', 'cold', 'insen', 'warm'))
transcriptome$transcript = gsub('(\\w+):\\d+', '\\1', transcriptome$ID)
transcriptome$ES = as.numeric(gsub('\\w+:(\\d+)', '\\1', transcriptome$ID))

transcriptome$dG_full = transcriptome$energy_full_structure_temp286.15.K. - transcriptome$energy_full_structure_temp295.15.K.
transcriptome$dG_local = transcriptome$dsRNA.ES._energy_temp286.15.K. - transcriptome$dsRNA.ES._energy_temp295.15.K. # values are how G shifts going from warm to cold

transcriptome$local_size_13 = nchar(transcriptome$dsRNA.ES._structure_temp286.15.K.)
transcriptome$local_size_22 = nchar(transcriptome$dsRNA.ES._structure_temp295.15.K.)
transcriptome$d_local_size = transcriptome$local_size_22 - transcriptome$local_size_13

library(stringr)
transcriptome$local_paired_size_13 = str_count(transcriptome$dsRNA.ES._structure_temp286.15.K., '\\)|\\(')
transcriptome$local_paired_size_22 = str_count(transcriptome$dsRNA.ES._structure_temp295.15.K., '\\)|\\(')
transcriptome$d_local_paired_size = transcriptome$local_paired_size_22 - transcriptome$local_paired_size_13

transcriptome = dplyr::arrange(transcriptome, temp)
transcriptome = transcriptome[-which(duplicated(transcriptome$ID)), ]

transcriptome = left_join(transcriptome, select(all_edits, `Genomic location`, Transcript, mean_warm, mean_cold, dEL), by = c('ID' = 'Transcript'))

transcriptome = left_join(transcriptome, select(genome, ID, dsRNA.ES._energy_temp286.15.K., dsRNA.ES._energy_temp295.15.K., dG_local), by = c('Genomic location' = 'ID'), suffix = c('', '.genome'))

# Test if there's a trend in local structure completely falling apart in the heat --------


tmp = subset(transcriptome, dsRNA.ES._energy_temp295.15.K. == 0)
tmp = with(tmp, table(temp, dsRNA.ES._energy_temp286.15.K. != 0))
tmp / rowSums(tmp)
chisq.test(tmp)





transcriptome = subset(transcriptome, dsRNA.ES._energy_temp295.15.K. != 0 & dsRNA.ES._energy_temp286.15.K. != 0)


pos_span = function(i){
	t = gsub('\\((\\d+), \\d+, \\d+, (\\d+)\\)', '\\1 \\2', i)
	return(unlist(sapply(strsplit(t, split = ' '), function(j) diff(as.numeric(j)))))
}

transcriptome$pos_span_13 = pos_span(transcriptome$dsRNA.ES._coords_temp286.15.K.)
transcriptome$pos_span_22 = pos_span(transcriptome$dsRNA.ES._coords_temp295.15.K.)
transcriptome$d_pos_span = transcriptome$pos_span_22 - transcriptome$pos_span_13





with(transcriptome, tapply(dG_local, temp, median))
with(transcriptome, tapply(d_local_paired_size, temp, mean))
with(transcriptome, tapply(d_local_size, temp, mean))

library(ggplot2)
theme_set(theme_bw())

tmp = transcriptome
tmp$temp = factor(tmp$temp, levels = c('cold10', 'cold', 'insen', 'warm'), labels = c('Cold-induced >10%', 'Cold-induced', 'Insignificant temperature sensitivity', 'Warm-induced'))

gg = ggplot(tmp, aes(fill = temp)) +
	facet_wrap(~temp, ncol = 1, scales = 'free_y') +
	scale_fill_manual(values = c('navy', '#4385FF', 'grey', '#FF654B')) +
	coord_cartesian(xlim = c(-50, 50), expand = FALSE) +
	theme(legend.position = 'none') +
	labs(x = expression(paste(Delta, G[dsRNA])), y = '# of editing sites')


gg_trans = gg + 
	geom_histogram(aes(x = dG_local), binwidth = 1) +
	geom_vline(xintercept = 0, lty = 3)

gg_gen = gg + 
	geom_histogram(aes(x = dG_local.genome), binwidth = 1) +
	geom_vline(xintercept = 0, lty = 3)

cowplot::plot_grid(gg_gen, gg_trans, labels = 'AUTO')

ggsave('../FigS6.pdf', width = 6, height = 6)

# notice that RNA tends to have lower G in the cold than warm regardless of category


ggplot(transcriptome, aes(d_local_paired_size, dG_local, color = temp)) +
	geom_point(shape = 1) +
	geom_smooth(method = 'lm')

# notice that colder-sensitive sites have a slightly lower slope than insensitive-sites

ggplot(transcriptome, aes(pos_span_13, dG_local, color = temp)) +
	geom_point(shape = 1)

# most structures span >100 nt and have a slight decline in G in the cold. But among structures spanning <100 nt, the smaller the structure the larger the rise in G in the cold.

ggplot(transcriptome, aes(dG_local, dEL, color = temp)) +
	geom_point(shape = 1)

ggplot(transcriptome, aes(dsRNA.ES._energy_temp295.15.K., mean_warm, color = temp)) +
	geom_point(shape = 1)



# How similar are transcriptome and exome? --------------------------------

ggplot(transcriptome, aes(dsRNA.ES._energy_temp295.15.K., dsRNA.ES._energy_temp295.15.K..genome)) +
	geom_point(shape = 1, alpha = 0.1) +
	geom_abline(slope = 1, intercept = 0)

cor(transcriptome$dsRNA.ES._energy_temp295.15.K., transcriptome$dsRNA.ES._energy_temp295.15.K..genome)

ggplot(transcriptome, aes(dG_local, dG_local.genome)) +
	geom_point(shape = 1, alpha = 0.1) +
	geom_abline(slope = 1, intercept = 0)

cor(transcriptome$dG_local, transcriptome$dG_local.genome)










library(dplyr)

d = subset(transcriptome, transcript %in% names(which(sort(table(transcriptome$transcript), decreasing = T) >= 100)))


ggplot(d[sample(nrow(d)), ], aes(dG_local, dEL, color = temp)) +
	geom_point() +
	facet_wrap(~transcript) +
	scale_color_manual(values = c('navy', 'blue', 'grey', 'red')) +
	theme(legend.position = 'none') +
	labs(x = '∆G of dsRNA around ES', y = '∆EL')


ggplot(d[sample(nrow(d)), ], aes(dsRNA.ES._energy_temp295.15.K., mean_warm, color = temp)) +
	geom_point() +
	facet_wrap(~transcript) +
	scale_color_manual(values = c('navy', 'blue', 'grey', 'red')) +
	theme(legend.position = 'none') +
	labs(x = 'G of dsRNA around ES at 22°C', y = 'EL at 22°C')


ggplot(d[sample(nrow(d)), ], aes(local_size_22, mean_warm, color = temp)) +
	geom_point() +
	facet_wrap(~transcript) +
	scale_color_manual(values = c('navy', 'blue', 'grey', 'red')) +
	theme(legend.position = 'none') +
	labs(x = 'Size of dsRNA around ES at 22°C', y = 'EL at 22°C')











# Kinesin and synaptotagmin -----------------------------------------------

filter(transcriptome, transcript == 'Ocbimv22000619m' & ES == '845')[, grep('dG', colnames(transcriptome))] # Kin

filter(transcriptome, transcript == 'Ocbimv22021175m' & ES == '742')[, grep('dG', colnames(transcriptome))]	# Syt		 

