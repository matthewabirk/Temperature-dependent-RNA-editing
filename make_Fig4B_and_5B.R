library(dplyr)
library(ggplot2)
# Import data -------------------------------------------------------------

load('~/Documents/WDs/R/Rosenthal/Octo_temp/data_objects.Rdata')

ELs_temp = readxl::read_excel('Data/Amplicon_EL.xlsx', sheet = 'cold_induced')
ELs_ctrl = readxl::read_excel('Data/Amplicon_EL.xlsx', sheet = 'temp_independent')

ELs = rbind(ELs_temp, ELs_ctrl)

tissue_log = readxl::read_excel('Data/tissues_sampled.xlsx', sheet = 1)
inter1 = data.frame(sample = c('baseline', 't=0', 't=8h', 't=24h', 't=48h', 't=4d', 't=96h'),
										temp_dur = c(-20, 0, 8, 24, 48, 4*24, 96))
tissue_log = dplyr::left_join(tissue_log, inter1)
d = dplyr::left_join(ELs, tissue_log)

d = d[which(apply(d, 1, function(i) sum(c(as.numeric(i['A_height']), as.numeric(i['G_height'])))) > 100), ] # remove all ELs where the A and G peaks were combined < 100

inter1 = d %>% group_by(dT, ES, temp_dur) %>% summarize(max_EL = max(EL), min_EL = min(EL))

d = dplyr::left_join(d, inter1)

d$transcript = gsub('_\\d+', '', d$ES)

d = subset(d, ES != 'Ank_1105')





# Get 6 week study results ------------------------------------------------

inter1 = c(Ank = 'Ocbimv22034836m',
					 Dyn = 'Ocbimv22033392m',
					 HN2E = 'Ocbimv22036107m',
					 KHC = 'Ocbimv22000619m',
					 Syt = 'Ocbimv22021175m',
					 Tit = 'Ocbimv22028814m')

long = lapply(unique(d[d$apriori_temp_sensitive, ]$ES), function(i){
	gene = inter1[gsub('_\\d+', '', i)]
	es = as.numeric(gsub('.*_', '', i))
	long_el = subset(all_temp, Transcript == gene & `Location in transcript` == es, c('EL_cold', 'EL_warm')) * 100
	cbind(data.frame(ES = i, apriori_temp_sensitive = TRUE), long_el)
})
long = do.call('rbind', long)
long$transcript = gsub('_\\d+', '', long$ES)



# Normalize ---------------------------------------------------------------

inter1 = plyr::ddply(d, c('dT', 'ES'), function(i) c(baseline = mean(subset(i, sample == 'baseline')$EL), end = mean(subset(i, sample %in% c('t=4d', 't=96h'))$EL)))
inter1$min_at_ends = as.numeric(apply(inter1, 1, function(i) min(c(i['baseline'], i['end']))))

d = dplyr::left_join(d, inter1)

d$EL_normalized = (d$EL - d$min_at_ends) / abs(d$baseline - d$end)


c('KHC_845', 'Syt_742')
d_soi = subset(d, ES %in% c('KHC_845'))
long_soi = subset(long, ES %in% unique(d_soi$ES))

d_soi$dT = ifelse(d_soi$dT == 'CW', '14 -> 24°C', '24 -> 14°C')

ggplot(d_soi, aes(temp_dur, EL)) +
	geom_rect(xmin = -20, xmax = 0, ymin = 0, ymax = 100, fill = 'lightgrey', color = 'lightgrey') +
	geom_segment(data = long_soi, aes(x = -20, xend = 96, y = EL_cold, yend = EL_cold), linetype = 3, color = '#4385FF') +
	geom_segment(data = long_soi, aes(x = -20, xend = 96, y = EL_warm, yend = EL_warm), linetype = 3, color = '#FF654B') +
	geom_point(aes(color = dT)) +
	geom_linerange(aes(x = temp_dur, ymin = min_EL, ymax = max_EL, color = dT)) +
	labs(x = 'Hours since end of temperature change', y = '% editing', color = NULL) +
	theme_classic() +
	geom_smooth(aes(color = dT), method = 'loess', span = 1.5, se = FALSE) +
	scale_x_continuous(breaks = sort(unique(d[d$apriori_temp_sensitive, ]$temp_dur))) +
	scale_color_manual(values = c('#FF654B', '#4385FF')) +
	theme(legend.position = 'top')
ggsave(filename = 'fig4b.pdf', width = 3.2, height = 3)

ggsave(filename = 'fig5b.pdf', width = 3.2, height = 3)



