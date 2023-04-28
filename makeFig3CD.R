# Makes Fig3C,D


library(dplyr)
library(ggplot2)
# Import data -------------------------------------------------------------

load('data_objects.Rdata')

ELs_temp = readxl::read_excel('Data/Amplicon_EL.xlsx', sheet = 'cold_induced')
ELs_ctrl = readxl::read_excel('Data/Amplicon_EL.xlsx', sheet = 'temp_independent')

ELs = rbind(ELs_temp, ELs_ctrl)

tissue_log = readxl::read_excel('/Users/matthewbirk/Documents/WDs/R/Rosenthal/Octo_temp/tissues_sampled.xlsx', sheet = 1)
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

dt = d[d$apriori_temp_sensitive, ]


dt_all = dt




# Statistics of 4-day study vs 6-week study -------------------------------

longterm = data.frame(
	stringsAsFactors = FALSE,
	ES = c("KHC_845","Syt_741","Syt_742",
				 "Tit_3635","Tit_3636","Tit_3646","Tit_3647",
				 "Tit_3660","Tit_6350","Tit_6598","Tit_6664","Tit_8774",
				 "Tit_8784","Tit_8934","Tit_15116","Tit_15330","Tit_15432",
				 "Dyn_2384"),
	exome_ID = c("Ocbimv22000619m",
							 "Ocbimv22021175m","Ocbimv22021175m","Ocbimv22028814m",
							 "Ocbimv22028814m","Ocbimv22028814m","Ocbimv22028814m",
							 "Ocbimv22028814m","Ocbimv22028814m","Ocbimv22028814m",
							 "Ocbimv22028814m","Ocbimv22028814m","Ocbimv22028814m",
							 "Ocbimv22028814m","Ocbimv22028814m","Ocbimv22028814m",
							 "Ocbimv22028814m","Ocbimv22033392m"),
	exome_ES = c(845,741,742,3635,3636,
							 3646,3647,3660,6350,6598,6664,8774,8784,
							 8934,15116,15330,15432,2384)
)

tmp = apply(longterm, 1, function(i){
	filter(cold, Transcript == i['exome_ID'] & `Location in transcript` == as.numeric(i['exome_ES'])) %>% dplyr::select(starts_with('Editing '))
})
longterm = cbind(longterm, do.call(rbind, tmp))
longterm = reshape2::melt(longterm, id.vars = c('ES', 'exome_ID', 'exome_ES'))
longterm$variable = gsub('Editing (\\w+)\\d', '\\1', longterm$variable)
colnames(longterm)[4:5] = c('temp', 'EL')
longterm$exp = '6 week'

shortterm_W = filter(dt_all, dT == 'CW', temp_dur %in% 96) %>% dplyr::select(ES, EL, temp_dur)
shortterm_C = filter(dt_all, dT == 'WC', temp_dur %in% 96) %>% dplyr::select(ES, EL, temp_dur)
shortterm_C$temp = 'cold'
shortterm_W$temp = 'warm'
shortterm = rbind(shortterm_W, shortterm_C)
shortterm$exp = '4 day'
shortterm$EL = shortterm$EL / 100

tmp = sapply(unique(longterm$ES), function(i){											# find the ∆EL for each site
	-diff((longterm %>% filter(ES == i) %>% group_by(temp) %>% summarise(EL = mean(EL)))$EL)
})
shortterm$EL_adj = unlist(apply(shortterm, 1, function(i){ # adjust 4 day EL down by 1-2°C based on linear scaling and ∆EL with the ∆T of 9°C
	if(i['temp'] == 'warm') EL_adj = as.numeric(i['EL']) - 2 / 9 * tmp[i['ES']]
	if(i['temp'] == 'cold') EL_adj = as.numeric(i['EL']) - 1 / 9 * tmp[i['ES']]
	return(EL_adj)
}))

qplot(shortterm$EL, shortterm$EL_adj, color = shortterm$temp) + geom_abline(slope = 1, intercept = 0)

colnames(longterm)[5] = 'EL_adj'

exp_comp = rbind(shortterm[, c('ES', 'EL_adj', 'temp', 'exp')], longterm[, c('ES', 'EL_adj', 'temp', 'exp')])

exp_comp = exp_comp %>% group_by(ES, temp, exp) %>% summarise(mean_EL = mean(EL_adj))

t.test(x = filter(exp_comp, temp == 'cold', exp == '4 day')$mean_EL,
			 y = filter(exp_comp, temp == 'cold', exp == '6 week')$mean_EL,
			 paired = TRUE, var.equal = TRUE)

t.test(x = filter(exp_comp, temp == 'warm', exp == '4 day')$mean_EL,
			 y = filter(exp_comp, temp == 'warm', exp == '6 week')$mean_EL,
			 paired = TRUE, var.equal = TRUE)






# T-tests to compare editing levels ---------------------------------------

dt = subset(dt_all, dT == 'CW') # run this as CW to get Fig3D

tmp = dt %>% group_by(ES, temp_dur) %>% summarise(mean_EL = mean(EL))
tmp4 = list()
for(i in 2:length(unique(dt$temp_dur))){
	tmp2 = sort(unique(dt$temp_dur))
	tmp3 = t.test(filter(tmp, temp_dur == tmp2[i])$mean_EL - filter(tmp, temp_dur == tmp2[i - 1])$mean_EL)
	tmp4[[i - 1]] = data.frame(dT = unique(dt$dT),
														 t_start = tmp2[i - 1],
														 t_end = tmp2[i],
														 diff = as.numeric(tmp3$estimate),
														 pvalue = tmp3$p.value
														 )
}
tmp4 = do.call('rbind', tmp4)
tmp4$padj = p.adjust(tmp4$pvalue, method = 'bonferroni')

write.table(tmp4, file = paste0('tableS5', unique(dt$dT), '.txt'), row.names = FALSE, quote = FALSE, sep = '\t')

tmp5 = dt %>% group_by(temp_dur) %>% summarise(EL = mean(EL)) # average EL pooling all ESs
tmp6 = longterm %>% filter(temp == ifelse(unique(dt$dT) == 'CW', 'warm', 'cold')) %>% group_by(ES) %>% summarise(EL = mean(EL_adj) * 100) 
tmp7 = dt %>% group_by(temp_dur, ES) %>% summarise(EL = mean(EL)) # average EL pooling all ESs

ggplot(dt,aes(temp_dur, EL, col = ES)) +
	geom_rect(xmin = -20, xmax = 0, ymin = 0, ymax = 100, fill = 'lightgrey', color = 'lightgrey') +
	geom_point() +
	geom_line(data = tmp7, aes(temp_dur, EL, col = ES), inherit.aes = FALSE) +
	geom_line(data = tmp5, aes(temp_dur, EL), lwd = 2, inherit.aes = FALSE) +
	geom_point(data = tmp6, aes(x = 100, y = EL), shape = 4) +
	labs(x = 'Time since temperature change (hr)', y = '% editing') +
	scale_x_continuous(breaks = c(0, unique(dt$temp_dur))) +
	theme_classic() +
	scale_color_discrete(guide = 'none')

ggsave(paste0('fig3_', unique(dt$dT), '.pdf'), width = 4.6, height = 4.4)










# Plot --------------------------------------------------------------------

ggplot(d[d$apriori_temp_sensitive, ], aes(temp_dur, EL)) +
	geom_rect(xmin = -20, xmax = 0, ymin = 0, ymax = 100, fill = 'lightgrey', color = 'lightgrey') +
	geom_segment(data = long, aes(x = -20, xend = 96, y = EL_cold, yend = EL_cold), linetype = 3, color = scales::hue_pal()(2)[2]) +
	geom_segment(data = long, aes(x = -20, xend = 96, y = EL_warm, yend = EL_warm), linetype = 3, color = scales::hue_pal()(2)[1]) +
	geom_point(aes(color = dT)) +
	geom_linerange(aes(x = temp_dur, ymin = min_EL, ymax = max_EL, color = dT)) +
	labs(x = 'Hours since end of temperature change', y = '% editing', color = NULL) +
	theme_classic() +
	geom_smooth(aes(color = dT), method = 'loess', span = 1.5, se = FALSE) +
	facet_wrap(~ES) +
	scale_x_continuous(breaks = sort(unique(d[d$apriori_temp_sensitive, ]$temp_dur))) +
	theme(legend.position = 'top')
ggsave(filename = 'amplicon_timecourse.pdf', width = 8.6, height = 6)


