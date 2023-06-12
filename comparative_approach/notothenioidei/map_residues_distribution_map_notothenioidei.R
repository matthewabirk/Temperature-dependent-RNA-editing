library(ggplot2)

# Syt multiple alignment --------------------------------------------------

syt_align = seqinr::read.fasta('notothenioidei/Notothenioidei.fa', seqtype = 'AA', as.string = TRUE)

syt_align = syt_align[1:19] 

library(stringr)

syt_align = data.frame(species = unname(gsub('.*\\[(.*)\\]', '\\1', sapply(syt_align, function(i) attr(i, 'Annot')))),
											 msa = as.character(unlist(syt_align)))
syt_align$residue = str_sub(syt_align$msa, start = 654, 654)

syt_align = syt_align[!duplicated(syt_align$species), ]

tmp = str_sub(syt_align$msa, start = 650, 660)



ggseqlogo::ggseqlogo(tmp, method = 'prob') +
	annotate('rect', xmin = 4.5, xmax = 5.5, ymin = -0.1, ymax = 1.1, alpha = 0.1, col = 'black', fill = 'yellow') +
	ggtitle('Alignment of notothenioidei matches to Ocbimv22021175m')

ggsave('notothenoid.pdf', width = 8, height = 3)








library(rgbif)
GBIF_USER = "matthewabirk"
GBIF_PWD = "XXX"
GBIF_EMAIL = "matthewabirk@gmail.com"


tmp = sapply(syt_align$species, function(i) name_backbone(i)$usageKey)
syt_align$taxonKey = sapply(tmp, function(i) if(is.null(i)) NA else i)

occ_download(pred_in("taxonKey", na.omit(syt_align$taxonKey)), 
						 pred("hasGeospatialIssue", FALSE),
						 pred("hasCoordinate", TRUE),
						 pred("occurrenceStatus","PRESENT"), 
						 pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
						 format = "SIMPLE_CSV",
						 user = GBIF_USER, pwd = GBIF_PWD, email = GBIF_EMAIL)

d = occ_download_get('0126049-230224095556074') %>%
	occ_download_import()

d_short = select(d, 'species', 'decimalLatitude', 'decimalLongitude')

d_short = left_join(d_short, syt_align[, c('species', 'residue')])

pdf('map.pdf', height = 10, width = 20)
ggplot() + 
	geom_polygon(data = map_data('world'), aes(long, lat, group = group)) + 
	geom_point(data = d_short, aes(decimalLongitude, decimalLatitude, color = residue)) +
	facet_wrap(~species) +
	coord_fixed() + 
	metR::scale_x_longitude() + 
	metR::scale_y_latitude()
dev.off()
