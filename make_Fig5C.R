
itc = data.frame(kd = c('KD1', 'KD1', 'KD2', 'KD2'), prot = ordered(rep(c('WT', 'I248V'), 2), levels = c('WT', 'I248V')), mean = c(129, 206, 842, 757), sd = c(10.8, 34.6, 290.5, 193), n = c(4, 6, 4, 6))

ggplot(itc, aes(kd, y = mean, fill = prot)) +
	geom_col(position = 'dodge') +
	geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.9), width = 0.2) +
	coord_cartesian(expand = FALSE, ylim = c(50, 1500)) +
	scale_fill_manual(values = c('#9DD4F1', '#D8BF99')) +
	scale_y_log10(breaks = c(50, 100, seq(0, 1500, by = 200))) +
	scale_x_discrete(labels = c(expression(K[D[1]]), expression(K[D[2]]))) +
	labs(x = element_blank(), y = expression(paste(K[D[Ca^{'2+'}]], ' (', mu, 'M)')), fill = element_blank())

ggsave('fig5c.pdf', width = 4, height = 3)

t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0)
{
	# pooled standard deviation, scaled by the sample sizes
	se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
	df <- n1+n2-2
	t <- (m1-m2-m0)/se 
	dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
	names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
	return(dat) 
}

library(dplyr)

itc1 = subset(itc, kd == 'KD1')
t.test2(m1 = filter(itc1, prot == 'WT')$mean,
				m2 = filter(itc1, prot == 'I248V')$mean,
				s1 = filter(itc1, prot == 'WT')$sd,
				s2 = filter(itc1, prot == 'I248V')$sd,
				n1 = filter(itc1, prot == 'WT')$n,
				n2 = filter(itc1, prot == 'I248V')$n
)

itc2 = subset(itc, kd == 'KD2')
t.test2(m1 = filter(itc2, prot == 'WT')$mean,
				m2 = filter(itc2, prot == 'I248V')$mean,
				s1 = filter(itc2, prot == 'WT')$sd,
				s2 = filter(itc2, prot == 'I248V')$sd,
				n1 = filter(itc2, prot == 'WT')$n,
				n2 = filter(itc2, prot == 'I248V')$n
)
