wc_a = read.csv('/Users/matthewbirk/Documents/WDs/R/Rosenthal/Octo_temp/Data/Condo A dT W->C.csv')
wc_b = read.csv('/Users/matthewbirk/Documents/WDs/R/Rosenthal/Octo_temp/Data/Condo B dT W->C.csv')
cw_c = read.csv('/Users/matthewbirk/Documents/WDs/R/Rosenthal/Octo_temp/Data/Condo C dT C->W.csv')
cw_d = read.csv('/Users/matthewbirk/Documents/WDs/R/Rosenthal/Octo_temp/Data/Condo D dT C->W.csv')

wc_a$dir = 'WC'
wc_b$dir = 'WC'
cw_c$dir = 'CW'
cw_d$dir = 'CW'

wc_a$sensor_id = 'A'
wc_b$sensor_id = 'B'
cw_c$sensor_id = 'C'
cw_d$sensor_id = 'D'

monnits = rbind(wc_a, wc_b, cw_c, cw_d)

monnits$Date = lubridate::mdy_hm(monnits$Date)

library(ggplot2)
ggplot(monnits, aes(Date, Value, color = sensor_id)) +
	geom_path() +
	facet_wrap(~dir, scales = 'free_x')
