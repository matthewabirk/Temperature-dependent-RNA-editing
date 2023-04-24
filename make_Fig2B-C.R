library(ggplot2); library(colorspace); library(dplyr); library(Biostrings)


theme_set(theme_bw())

load('data_objects.Rdata')


# Recoding bias -----------------------------------------------------------

tmp = subset(temp_insensitive, temp_insensitive$`Codon changes` != 'syn')

inter1 = rbind(data.frame(ES = warm_recod$ES, type = 'warm', change = warm_recod$change, EL_cold = warm_recod$EL_cold, EL_warm = warm_recod$EL_warm),
							 data.frame(ES = gsub(':', '_', tmp$Transcript), type = 'insen', change = tmp$`Codon changes`, EL_cold = tmp$`mean cold`, EL_warm = tmp$`mean warm`),
							 data.frame(ES = cold_recod$ES, type = 'cold', change = cold_recod$change, EL_cold = cold_recod$EL_cold, EL_warm = cold_recod$EL_warm))

inter1$type = ordered(inter1$type, levels = c('warm', 'insen', 'cold'))

inter1$change = as.character(inter1$change)

inter2 = as.data.frame(table(inter1[, c('type', 'change')]))
inter2$prop = inter2$Freq / as.numeric(table(inter1$type))

inter2$type = ordered(inter2$type, levels = c('warm', 'insen', 'cold'))

table(inter1$type) # this is how many sites represent each category

chisq.test(table(inter1$type, inter1$change))


aa_sub = ggplot(inter2, aes(change, y = prop, fill = type)) +
	geom_col(position = 'dodge') +
	labs(x = 'AA subsitution', y = 'Proportion of all recoding events', fill = NULL) +
	scale_y_continuous(expand = c(0, 0), limits = c(0, 0.23)) +
	scale_x_discrete() +
	scale_fill_manual(values = c('#FF654B', 'gray', '#4385FF'), labels = c('Warm-induced', 'Temp. insensitive', 'Cold-induced')) +
	geom_vline(xintercept=seq(1.5, length(unique(inter2$change))-0.5, 1), lwd = 0.25, color = 'grey50', lty = 3) +
	theme(panel.grid.major.x = element_blank(),
				legend.position = 'none')

aa_sub
ggsave('figS2.pdf', width = 8, height = 3)

tmp = data.frame(change = unique(inter2$change))
tmp$p = sapply(unique(inter2$change), function(i){
	prop.trend.test(subset(inter2, change == i)$Freq, table(inter1$type))$p.value
})

tmp$p_adj = p.adjust(tmp$p, 'bonferroni')
subset(tmp, p_adj < 0.05)



# show how many sites represent each category
tmp1 = as.data.frame(table(inter1$type))

#geom_label x argument adjusts position within the slice

ggplot(tmp1, aes (x = '', y = Freq, fill = Var1)) + 
	geom_col(position = 'stack', width = 1, color = 'white') +
	geom_label(aes(label = Freq, x = c(0.9, 1, 1)), position = position_stack(vjust = 0.5), color = 'white') +
	theme_void() +
	coord_polar('y') +
	theme(legend.position = 'none') +
	scale_fill_manual(values = c('red', 'gray', 'blue'), labels = c('Warm-induced', 'Temp. insensitive', 'Cold-induced'))






# BLOSUM ------------------------------------------------------------------


library(Biostrings)
data(BLOSUM80)

inter1$blosum = sapply(inter1$change, function(i) BLOSUM80[substr(i, 1, 1), substr(i, 2, 2)])
inter1$sign = sign(inter1$blosum)
inter1[which(inter1$sign == -1), 'sign'] = '(-)'
inter1[which(inter1$sign == 1), 'sign'] = '(+)'
inter1$sign = ordered(inter1$sign, levels = c('(-)', '0', '(+)'))


full_table = table(inter1$type, inter1$sign)

chisq.test(full_table)
chisq.posthoc.test::chisq.posthoc.test(full_table, method = 'bonferroni')
chisq.test(full_table[c('warm', 'insen'), ])$p.value * 3
chisq.test(full_table[c('insen', 'cold'), ])$p.value * 3
chisq.test(full_table[c('warm', 'cold'), ])$p.value * 3

inter2 = as.data.frame(full_table)
inter2$prop = inter2$Freq / as.numeric(table(inter1$type))

aa_blosum = ggplot(inter2, aes(Var2, y = prop, fill = Var1)) +
	geom_col(position = 'dodge') +
	labs(x = 'BLOSUM80 score', y = 'Proportion of all recoding events') +
	scale_y_continuous(expand = c(0, 0), limits = c(0, 0.65)) +
	scale_fill_manual(values = c('#FF654B', 'gray', '#4385FF'), labels = c('Warm-induced', 'Temp. insensitive', 'Cold-induced')) +
	geom_vline(xintercept = seq(1.5, length(unique(inter2$Var2))-0.5, 1), lwd = 0.25, color = 'grey50', lty = 3) +
	theme(panel.grid.major.x = element_blank(),
				legend.position = 'none')

aa_blosum

tmp = data.frame(change = unique(inter2$Var2))
tmp$p = sapply(unique(inter2$Var2), function(i){
	prop.trend.test(subset(inter2, Var2 == i)$Freq, table(inter1$type))$p.value
})

tmp$p_adj = p.adjust(tmp$p, 'bonferroni')
subset(tmp, p_adj < 0.05)





t.test(x = filter(inter1, type == 'warm')$blosum,
			 y = filter(inter1, type == 'insen')$blosum,
			 var.equal = TRUE)$p.value * 3
t.test(x = filter(inter1, type == 'warm')$blosum,
			 y = filter(inter1, type == 'cold')$blosum,
			 var.equal = TRUE)$p.value * 3
t.test(x = filter(inter1, type == 'insen')$blosum,
			 y = filter(inter1, type == 'cold')$blosum,
			 var.equal = TRUE)$p.value * 3



# Polarity groups ---------------------------------------------------------


traits = readxl::read_xlsx('/Users/matthewbirk/Documents/WDs/R/amino_acid_traits.xlsx')

inter1$group = sapply(inter1$change, function(i){
	aa = strsplit(i, '')[[1]]
	g1 = traits[aa[1] == traits$amino_acid, 'group']
	g2 = traits[aa[2] == traits$amino_acid, 'group']
	change = ifelse(g1 == g2, 'Same', paste(g1, g2, sep = '>'))
})

inter1$group_switch = ifelse(inter1$group == 'Same', 'Same', 'Changed')

full_table = table(inter1$type, inter1$group_switch)

chisq.test(full_table)
chisq.posthoc.test::chisq.posthoc.test(full_table, method = 'bonferroni')
chisq.test(full_table[c('warm', 'insen'), ])$p.value * 3
chisq.test(full_table[c('insen', 'cold'), ])$p.value * 3
chisq.test(full_table[c('warm', 'cold'), ])$p.value * 3



inter2 = as.data.frame(full_table)
inter2$prop = inter2$Freq / as.numeric(table(inter1$type))

aa_group = ggplot(inter2, aes(Var2, y = prop, fill = Var1)) +
	geom_col(position = 'dodge') +
	labs(x = 'Polarity-based category', y = 'Proportion of all recoding events') +
	scale_y_continuous(expand = c(0, 0), limits = c(0, 0.95)) +
	scale_fill_manual(values = c('#FF654B', 'gray', '#4385FF'), labels = c('Warm-induced', 'Temp. insensitive', 'Cold-induced')) +
	geom_vline(xintercept=seq(1.5, length(unique(inter2$Var2))-0.5, 1), lwd = 0.25, color = 'grey50', lty = 3) +
	theme(legend.position = 'none',
				panel.grid.major.x = element_blank())

aa_group

tmp = data.frame(change = unique(inter2$Var2))
tmp$p = sapply(unique(inter2$Var2), function(i){
	prop.trend.test(subset(inter2, Var2 == i)$Freq, table(inter1$type))$p.value
})

tmp$p_adj = p.adjust(tmp$p, 'bonferroni')
subset(tmp, p_adj < 0.05)





tmp = cowplot::plot_grid(aa_group, aa_blosum, nrow = 1)
cowplot::plot_grid(tmp, aa_sub, nrow = 2)

ggsave('fig2B-C.pdf', width = 8, height = 6)



write.table(inter1, file = '~/Downloads/bias_test.txt', quote = FALSE, row.names = FALSE, sep = '\t')
