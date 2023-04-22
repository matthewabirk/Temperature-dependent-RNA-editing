samples = list.files('make_FigS1/Run 1', full.names = TRUE)
layout(matrix(1:12, ncol = 3, byrow = TRUE))
for(i in samples){
	log = read.csv(i, skip = 1)
	log$date = lubridate::mdy_hms(log$Date.Time..GMT.04.00)
	log = log[log$date <= lubridate::ymd_hms('2015-08-03 16:00:00'), ]
	log$temp = log[, grep('Temp', colnames(log))]
	
	title = paste('Experiment_1', gsub('.csv', '', basename(i)))
	
	plot(log$date, log$temp, type = 'l', main = title, xlab = 'Date', ylab = 'Temperature (°C)', ylim = c(12, 24))
	
	tmp = range(log[log$temp < 14, 'date'])
	print(difftime(tmp[2], tmp[1]))
	print(mean(log[log$temp < 14, 'temp']))
}



samples = list.files('make_FigS1/Run 2', full.names = TRUE)
for(i in samples){
	log = read.csv(i)
	log$Date = lubridate::mdy_hm(log$Date)
	log = log[log$Date <= lubridate::ymd_hms('2017-08-31 15:02:00'), ]
	
	title = paste('Experiment_2', gsub('.csv', '', basename(i)))
	
	plot(log$Date, log$Value, type = 'l', main = title, xlab = 'Date', ylab = 'Temperature (°C)', ylim = c(12, 24))
	
	tmp = range(log[log$Value < 14, 'Date'])
	print(difftime(tmp[2], tmp[1]))
	print(mean(log[log$Value < 14, 'Value']))
}

dev.copy2pdf(file = 'figS1.pdf', height = 10, width = 12)
