setwd('./Sep2019')

campground = read.csv('Campground.csv')
campground$site = 'Campground'

Chuck = read.csv('Chuck_logger.csv')
Chuck$site = 'Chuck'

Isthmus_reef = read.csv('Isthmus_reef.csv')
Isthmus_reef$site = 'Isthmus_reef'

FoJ = read.csv('Fourth_of_July_E9.csv')
FoJ$site = 'FoJ_E9'

Birk_1 = read.csv('../Feb2020/Birk_1.csv')
Birk_1$site = 'Birk_1'
Birk_2 = read.csv('../Feb2020/Birk_2.csv')
Birk_2$site = 'Birk_2'
Birk_3 = read.csv('../Feb2020/Birk_3.csv')
Birk_3$site = 'Birk_3'
Birk_4 = read.csv('../Feb2020/Birk_4.csv')
Birk_4$site = 'Birk_4'

LB22 = read.csv('../Feb2022/Long_Beach_Winter_2021-22.csv')
LB22$site = 'LongBeach_Feb2022'

hobo = rbind(campground, Chuck, Isthmus_reef, FoJ, Birk_1, Birk_2, Birk_3, Birk_4, LB22)

hobo = hobo[hobo$in_situ, ]

hobo$local_time = lubridate::mdy_hm(hobo$local_time)










hobo_compare = subset(hobo, site %in% c('Birk_2', 'FoJ_E9', 'Chuck', 'LongBeach_Feb2022'))

hobo_compare[which(hobo_compare$site == 'FoJ_E9'), 'species'] = 'O. bimaculatus'
hobo_compare[which(hobo_compare$site == 'Birk_2'), 'species'] = 'O. bimaculatus'
hobo_compare[which(hobo_compare$site == 'Chuck'), 'species'] = 'O. bimaculoides'
hobo_compare[which(hobo_compare$site == 'LongBeach_Feb2022'), 'species'] = 'O. bimaculoides'

hobo_compare[which(hobo_compare$site == 'FoJ_E9'), 'site'] = 'Sep 2019'
hobo_compare[which(hobo_compare$site == 'Birk_2'), 'site'] = 'Feb 2020'

hobo_compare[which(hobo_compare$site == 'Chuck'), 'site'] = 'Sep 2019'
hobo_compare[which(hobo_compare$site == 'LongBeach_Feb2022'), 'site'] = 'Feb 2022'


library(ggplot2)

ggplot(hobo_compare, aes(local_time, T, color = site)) +
	geom_line() +
	scale_x_datetime(date_breaks = '1 week') +
	theme_bw() +
	theme(axis.text.x = element_text(hjust = 0, angle = 330), legend.position = 'none') +
	facet_grid(species~site, scales = 'free_x') +
	scale_color_manual(values = c('#4385FF', '#4385FF', '#FF654B', '#FF654B')) +
	labs(x = '', y = expression(Temperature~(degree*C)))



ggsave('../../figS4.pdf', width = 7.5, height = 3.5)

with(hobo_compare, tapply(T, site, mean))
with(hobo, tapply(T, site, mean))


by(hobo_compare, hobo_compare[, c('site', 'species')], function(i) difftime(range(i$local_time)[1], range(i$local_time)[2])) / 7 # how many weeks?
