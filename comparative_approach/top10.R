library(dplyr)
load('../data_objects.Rdata')

d = subset(cold_recod, !(Transcript %in% c('Ocbimv22028502m', 'Ocbimv22033634m')))

d = tail(arrange(d, EL_diff), 10)

write.table(d$Transcript, 'top10.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
