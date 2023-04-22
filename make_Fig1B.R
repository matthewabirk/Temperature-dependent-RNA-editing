library(ggplot2); library(colorspace); library(dplyr); library(Biostrings)


theme_set(theme_bw())

load('data_objects.Rdata')



# Pie chart ---------------------------------------------------------------

inter1 = all_temp$Differential
inter1[which(inter1 == 'COLD')] = 'Cold-induced'
inter1[which(inter1 == 'WARM')] = 'Warm-induced'
inter2 = c(inter1, rep('No effect', times = nrow(temp_insensitive)))
inter3 = data.frame(type = inter2)



tmp1 = as.data.frame(table(inter3))
tmp1 = tmp1 %>% 
	dplyr::arrange(desc(tmp1)) %>%
	mutate(prop = Freq / sum(tmp1$Freq) * 100) %>%
	mutate(ypos = cumsum(prop) - 0.5 * prop)
ggplot(tmp1, aes(x = '', y = prop, fill = type)) +
	geom_bar(stat = 'identity', color = 'white') +
	geom_label(aes(y = ypos, label = Freq), color = 'white') +
	coord_polar('y', start = 0) +
	scale_fill_manual(values = c('#4385FF', 'grey', '#FF654B'), name = NULL) +
	theme_void() +
	theme(legend.position = 'none')

ggsave('fig1B.pdf', width = 2.75, height = 2.75)

