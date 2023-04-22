d = readxl::read_xlsx('Data/ADAR_Salmon_Temp2.xlsx')

d[grep('Cold', d$Sample), 'temp'] = 13
d[grep('Warm', d$Sample), 'temp'] = 22

d$temp = as.factor(d$temp)

library(ggplot2); library(ggsignif)
theme_set(theme_bw())

ggplot(d, aes(temp, TPM, fill = temp)) +
	geom_boxplot() +
	geom_point(shape = 21) +
	geom_signif(comparisons = list(1:2), test = 't.test', test.args = list(var.equal = TRUE), map_signif_level = c('****' = 0.0001, '***' = 0.001, '**' = 0.01, '*' = 0.05), color = 'black', textsize = 2.5) +
	scale_y_continuous(limits = c(0, NA)) +
	facet_wrap(~Enzyme, scales = 'free_y') +
	labs(x = expression(Temperature~(degree*C)), y = 'Transcripts per million (TPM)') +
	scale_fill_manual(values = c('#4385FF', '#FF654B')) +
	theme(legend.position = 'none')

summary(lm(TPM ~ temp, data = subset(d, Enzyme == 'ADAR1')))
summary(lm(TPM ~ temp, data = subset(d, Enzyme == 'ADAR2b')))

ggsave('FigS5.pdf', width = 3.5, height = 3)
