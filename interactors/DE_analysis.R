setwd('interactors')

library('DESeq2'); library('tximport'); library('stringr'); library(dplyr)

quant_files = list.files('../Data/temp_sg_salmon_quant', full.names = TRUE)

sample_info = data.frame(sample_id = gsub('temp_sg_salmon_quant_(.*).txt', '\\1', basename(quant_files)))
sample_info$temperature = factor(gsub('sample\\d\\d_(\\w+)\\d', '\\1', sample_info$sample_id), levels = c('warm', 'cold'))
sample_info$batch_id = factor(c(rep(1, 6), rep(2, 6)))
rownames(sample_info) = sample_info$sample_id



names(quant_files) = gsub('temp_sg_salmon_quant_(.*)\\.txt', '\\1', basename(quant_files))





tmp = read.table(quant_files[1], header = TRUE)$Name
t2g = data.frame(tmp, tmp)

## Read Salmon abundances
txi = tximport(files = quant_files, type = 'salmon', tx2gene = t2g, countsFromAbundance = 'lengthScaledTPM')
				# the countsFromAbundance parameter changes the counts matrix in the resulting txi file


# generate DESeq data set (dds)
dds = DESeqDataSetFromTximport(txi, sample_info, design = ~temperature)


# PCA ---------------------------------------------------------------------

vst = varianceStabilizingTransformation(dds) # needed to make the variance homoskedastic, which is requirement for PCA. Also log-transforms because log-transformed read count data is roughly normally distributed.

plotPCA(vst, intgroup = 'temperature') # plots PCA for top 500 most variable genes

plotPCA(vst, intgroup = 'batch_id') # batch does not seem to be an important driver of variation. Good! Probably don't need to account for in model.





# Do hierarchical clustering ----------------------------------------------
# this method breaks samples up into a free number of clusters

vst_counts = assay(vst) # output variance stabilizing transformed count matrix
dist_matrix = dist(t(vst_counts)) # calculates Euclidean distance between samples
plot(hclust(dist_matrix))

library(pheatmap)
pheatmap(cor(vst_counts), annotation_col = sample_info['temperature']) 


# Do K-means clustering ---------------------------------------------------
# this method breaks samples up into a pre-defined number of clusters

kmeans(t(vst_counts), centers = 2)$cluster # break into 2 clusters





# -------------------------------------------------------------------------


dds = DESeq(dds) # normalize for library size, adjust dispersions, and run stats



sample_info$sizefactor = sizeFactors(dds) # QC: confirm they roughly match variation in sequencing depth (though GC bias will also skew)



plotDispEsts(dds) # QC: observe that dispersion (variance) decreases as read depth increases. Typical for RNASeq experiments. Blue are the new values that replace the black values. Blue hollow circles highlight outliers that estimateDispersions() does not attempt to change.

plotSparsity(dds) # See how evenly dispersed counts are across samples vs being mainly all from one sample.

results(dds, contrast = c('temperature', 'cold', 'warm'), alpha = 0.05)
	# QC: if LFC, pval, and padj are all NA, then the gene has zero counts in all samples and was removed before testing.
	# QC: if pval and padj are NA, then the gene has an extreme outlier and was removed before testing. However, because there are so many samples, the outlier was replaced with the average instead.
	# QC: if padj is NA, then the gene has a low mean counts and was removed before testing to improve multiple-comparisons sensitivity.
	# The LFCs here are too high, especially those genes with high dispersion and/or low counts. They need to be shrunken first.

# -------------------------------------------------------------------------

DE_results = lfcShrink(dds, coef = 'temperature_cold_vs_warm', type = 'apeglm')

plotMA(DE_results)

#DE_results = DE_results[complete.cases(DE_results), ] # get rid of genes with NAs (low counts or outliers)
#plotMA(DE_results)




# Inspect DE genes ---------------------------------------------------------

summary(DE_results, alpha = 0.05)


DE_genes = as.data.frame(DE_results) %>% filter(padj < 0.05)
DE_genes$gene_id = rownames(DE_genes)
rownames(DE_genes) = NULL
colnames(DE_genes)[2:3] = c('log2FoldChange-cold_v_warm', 'lfcSE-cold_v_warm')



# How about SOIs (ADAR interactors)? --------------------------------------


SOIs = read.table('Ob_SOIs.txt', col.names = 'seq_ID')

tmp = lapply(SOIs$seq_ID, function(i){
	if(i %in% row.names(DE_results)) as.data.frame(DE_results[i, ]) else NA
})
SOIs = do.call('rbind', tmp)
SOIs = na.omit(SOIs)

table(ifelse(SOIs$log2FoldChange < 0, '-', '+'), SOIs$padj < 0.05)



# Connect biogrid data to DE data ----------------------------------------


SOIs = cbind(Ob_ID = row.names(SOIs), SOIs)
row.names(SOIs) = NULL

blast_results = read.table('blast_results.txt')[, c(1:2)]
blast_results$V1 = gsub('sp\\|(.*)\\|.*', '\\1', blast_results$V1)
colnames(blast_results) = c('interactor_sp_ID', 'Ob_ID')

biogrid_data = read.table('good_biogrid_data.txt', sep = '\t', col.names = c('BioGRID_Interaction_ID', 'A', 'B', 'Experimental_system', 'Author', 'PubMed_ID', 'A_sp_ID', 'B_sp_ID', 'A_species', 'B_species'))

biogrid_data$interactor_sp_ID = apply(biogrid_data, 1, function(i){
	grep('P55265|P78563', i[c('A_sp_ID', 'B_sp_ID')], value = TRUE, invert = TRUE)
})

biogrid_data$interactor_name = apply(biogrid_data, 1, function(i){
	grep('ADAR|ADARB1', i[c('A', 'B')], value = TRUE, invert = TRUE)
})

biogrid_data = left_join(biogrid_data, blast_results, by = 'interactor_sp_ID')

biogrid_data = left_join(biogrid_data, SOIs, by = 'Ob_ID')

write.table(biogrid_data, 'DE_results_biogrid.txt', quote = FALSE, sep = '\t', row.names = FALSE)

# Connect Freund data to DE data ----------------------------------------

freund_data = read.table('freund_sp_IDs.txt', col.names = c('interactor_sp_ID', 'interactor_name'))

freund_data = left_join(freund_data, blast_results, by = 'interactor_sp_ID')

freund_data = left_join(freund_data, SOIs, by = 'Ob_ID')

write.table(freund_data, 'DE_results_freund.txt', quote = FALSE, sep = '\t', row.names = FALSE)






biogrid_data = cbind(dataset = 'BioGRID', biogrid_data)
freund_data = cbind(dataset = 'Freund et al. 2020', freund_data)

all_data = bind_rows(biogrid_data, freund_data)

write.table(all_data, 'DE_results.txt', quote = FALSE, sep = '\t', row.names = FALSE)

# Make figure --------------------------------------------------------

SOIs$db_origin = factor(unlist(sapply(SOIs$Ob_ID, function(i){
	in_biogrid = nrow(filter(biogrid_data, Ob_ID == i)) > 0
	in_freund = nrow(filter(freund_data, Ob_ID == i)) > 0
	if(in_biogrid & !in_freund) return('BioGRID')
	if(in_freund & !in_biogrid) return('Freund et al. 2020')
	if(in_biogrid & in_freund) return('Both')
	if(!in_biogrid & !in_freund) return('Neither')
})), levels = c('BioGRID', 'Freund et al. 2020', 'Both', 'Neither'))

SOIs = filter(SOIs, db_origin != 'Neither')

library(scales)

theme_set(theme_bw())

SOIs$alpha = SOIs$padj < 0.05

ggplot(SOIs, aes(log2FoldChange, padj, alpha = alpha, color = db_origin)) +
	geom_point(pch = 16) +
	geom_hline(yintercept = 0.05, lty = 2) +
	scale_alpha_discrete(range = c(0.25, 1), guide = 'none') +
	scale_y_continuous(trans = c("log10", "reverse")) +
	scale_color_brewer(type = 'qual', palette = 'Dark2') +
	labs(x = 'log2-fold change from warm to cold', y = 'Adjusted p-value', color = 'Source') +
	theme(legend.position = c(0.3, 0.75))

ggsave('../figS7.pdf', width = 3.5, height = 3)
