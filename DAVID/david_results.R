david_results = read.table('/Users/matthewbirk/Documents/Career/Rosenthal Lab/Temperature Octopus/DAVID/david_functional_annotation_chart_cold10_UP_keywords_only.txt', sep = '\t', header = TRUE)

signif_david_results = subset(david_results, FDR < 0.05)

View(signif_david_results)

write.csv(david_results, 'UP_keywords_only.csv')
