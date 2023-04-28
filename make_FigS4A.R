library(ggplot2)

# KHC multiple alignment --------------------------------------------------

khc_align = seqinr::read.fasta('Data/multiple alignment of Kinesin heavy chain - Ocbimv22000619m - top 250.fa', seqtype = 'AA', as.string = TRUE)

library(stringr)
str_split_list = function(string,size){
	str_extract_all(string, paste0('.{1,',size,'}'))
}

inter1 = sapply(khc_align, function(i) str_split_list(i, 25))
inter2 = lapply(1:unique(sapply(inter1, length)), function(i){
	as.character(sapply(inter1, function(j) j[i]))
})



ggseqlogo::ggseqlogo(inter2[[ceiling(539/25)]], method = 'prob') +
	annotate('rect', xmin = 13.5, xmax = 14.5, ymin = -0.1, ymax = 1.1, alpha = 0.1, col = 'black', fill = 'yellow') +
	ggtitle('Alignment of top 250 matches to Ocbimv22000619m (cephs to mammals)\nMotor domain site')

ggsave('figS4A.pdf', width = 8, height = 3)
