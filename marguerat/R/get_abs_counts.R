dir.create('../abs_counts', showWarnings = F)
url <- 'https://ars.els-cdn.com/content/image/1-s2.0-S0092867412011269-mmc1.xlsx'
download.file(url, '../abs_counts/mmc1.xlsx')
table <- openxlsx::read.xlsx('../abs_counts/mmc1.xlsx', sheet = 3, startRow = 2)

table <- dplyr::select(table, 1, 2, 11, 14, 17, 20)
table$fold_change <- log2(rowMeans(as.matrix(table[,5:6])) / rowMeans(as.matrix(table[,3:4])))

file_name <- '~/sleuth_paper_analysis/annotation/ASM294v2.pombase.all_annos.txt'
yeast_annos <- read.table(file_name, header = T, sep = "\t", quote = "",
                          stringsAsFactors = FALSE)

cov <- matrixStats::rowSds(as.matrix(table[,3:6])) / rowMeans(as.matrix(table[,3:6]))
names(cov) <- table$Systematic.name
denom_gene <- names(which(cov == min(cov[which(cov>0)], na.rm = T)))

denom <- yeast_annos$target_id[which(yeast_annos$gene_id == denom_gene)]
write.table(table, '../abs_counts/abs_counts.txt', sep = "\t", row.names = F, quote = F)
saveRDS(denom, '../abs_counts/denom.rds')
