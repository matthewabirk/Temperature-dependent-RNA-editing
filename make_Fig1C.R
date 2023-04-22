library(ggplot2); library(colorspace); library(dplyr); library(Biostrings)



theme_set(theme_bw())

load('data_objects.Rdata')


tmp = all_temp
tmp$dEL = ifelse(tmp$Differential == 'COLD', tmp$EL_diff * 100, tmp$EL_diff * -100)
tmp = hist(tmp$dEL, breaks = seq(-101.25, 100, by = 2.5), plot = FALSE)
tmp = data.frame(mid = tmp$mids, n = tmp$counts)
tmp[which(tmp$mid == 0), 'n'] = tmp[which(tmp$mid == 0), 'n'] + nrow(temp_insensitive)

ggplot(tmp, aes(mid, n, fill = mid)) +
	geom_col() +
	scale_fill_continuous_diverging(palette = 'Blue-Red 2', mid = 0, p1 = 0.5, p2 = 0.5, rev = TRUE) +
	labs(x = expression(paste(Delta, '% editing')), y = '# of editing sites') +
	scale_x_continuous(breaks = seq(-100, 100, by = 10), limits = range(subset(tmp, n != 0, 'mid'))) +
	coord_cartesian(ylim = c(0, 8500), expand = FALSE) +
	theme(legend.position = 'none')

ggsave('fig1C.pdf', width = 4, height = 2.5)
