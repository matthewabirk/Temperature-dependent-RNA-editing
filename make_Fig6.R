d = readxl::read_xlsx('Data/Catalina_amplicon_EL.xlsx')

d[which(d$ES == 'KHC_845'), 'ES'] = 'Kinesin K282R'
d[which(d$ES == 'Syt_742'), 'ES'] = 'Synaptotagmin I248V'

loides = dplyr::filter(d, species == 'O. bimaculoides')
latus = dplyr::filter(d, species == 'O. bimaculatus')

loides[which(loides$trip == 'Sep 2019'), 'temp'] = 21
loides[which(loides$trip == 'Feb 2022'), 'temp'] = 15

latus[which(latus$trip == 'Sep 2019'), 'temp'] = 22
latus[which(latus$trip == 'Birk_2'), 'temp'] = 16

loides$temp = as.factor(loides$temp)
latus$temp = as.factor(latus$temp)

library(ggplot2); library(ggsignif)
ggplot(loides, aes(temp, EL * 100, fill = trip)) + 
	geom_boxplot() +
	geom_point(shape = 21) +
	geom_signif(comparisons = list(1:2), test = 't.test', test.args = list(var.equal = TRUE), map_signif_level = c('****' = 0.0001, '***' = 0.001, '**' = 0.01, '*' = 0.05), color = 'black', vjust = 0.5) +
	facet_grid(species~ES) +
	scale_fill_manual(values = c('#4385FF', '#FF654B')) +
	labs(x = 'Temperature (°C)', y = '% edited') +
	theme_bw() +
	theme(legend.position = 'none')

ggsave('fig6A.pdf', width = 3.5, height = 3)

ggplot(latus, aes(temp, EL * 100, fill = trip)) + 
	geom_boxplot() +
	geom_point(shape = 21) +
	geom_signif(comparisons = list(1:2), test = 't.test', test.args = list(var.equal = TRUE), map_signif_level = c('****' = 0.0001, '***' = 0.001, '**' = 0.01, '*' = 0.05), color = 'black', vjust = 0.5) +
	facet_grid(species~ES) +
	scale_fill_manual(values = c('#4385FF', '#FF654B')) +
	labs(x = 'Temperature (°C)', y = '% edited') +
	theme_bw() +
	theme(legend.position = 'none')

ggsave('fig6B.pdf', width = 3.5, height = 3)





library(lme4)
null = lmer(formula = EL ~ (1|ES), data = loides, REML = FALSE)
mod = lmer(formula = EL ~ trip + (1|ES), data = loides, REML = FALSE)
anova(null, mod)

summary(lm(EL ~ temp, data = subset(loides, ES == 'Kinesin K282R')))
summary(lm(EL ~ temp, data = subset(loides, ES == 'Synaptotagmin I248V')))


library(lme4)
null = lmer(formula = EL ~ (1|ES), data = latus, REML = FALSE)
mod = lmer(formula = EL ~ trip + (1|ES), data = latus, REML = FALSE)
anova(null, mod)

summary(lm(EL ~ temp, data = subset(latus, ES == 'Kinesin K282R')))
summary(lm(EL ~ temp, data = subset(latus, ES == 'Synaptotagmin I248V')))
