cw = data.frame(
	dur = c(-35, -30, -30, -20, 0, 8, 24, 96),
	temp = c(14, 14, 14, 14, 24, 24, 24, 24),
	ltype = factor(c(2, 2, 1, 1, 1, 1, 1, 1))
)

library(ggplot2)

ggplot(cw, aes(dur, temp)) +
	geom_vline(xintercept = cw$dur[-1:-3], color = 'grey', linetype = 3) +
	geom_path(aes(lty = ltype)) +
	labs(x = 'Time since temperature change (hr)', y = 'Temperature (°C)') +
	scale_x_continuous(breaks = cw$dur[-1:-3]) +
	theme_classic() +
	theme(legend.position = 'none')
ggsave('/Users/matthewbirk/Documents/Career/Writing/Manuscripts in prep/Temp sensitive RNA editing/Figures/temp_timeline-CW.pdf', width = 4, height = 2)







wc = data.frame(
	dur = c(-35, -30, -30, -20, 0, 24, 48, 96),
	temp = c(24, 24, 24, 24, 14, 14, 14, 14),
	ltype = factor(c(2, 2, 1, 1, 1, 1, 1, 1))
)

library(ggplot2)

ggplot(wc, aes(dur, temp)) +
	geom_vline(xintercept = wc$dur[-1:-3], color = 'grey', linetype = 3) +
	geom_path(aes(lty = ltype)) +
	labs(x = 'Time since temperature change (hr)', y = 'Temperature (°C)') +
	scale_x_continuous(breaks = wc$dur[-1:-3]) +
	theme_classic() +
	theme(legend.position = 'none')
ggsave('/Users/matthewbirk/Documents/Career/Writing/Manuscripts in prep/Temp sensitive RNA editing/Figures/temp_timeline-WC.pdf', width = 4, height = 2)
