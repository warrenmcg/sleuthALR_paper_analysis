.libPaths('~/R_library')
library("dplyr")

url <- "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1154/E-MTAB-1154.sdrf.txt"

dir.create('../metadata', showWarnings = F)
file_name <- file.path('../metadata', 'E-MTAB-1154.sdrf.txt')
download.file(url = url, destfile = file_name, method = 'auto')

table <- read.table(file_name, sep = "\t", header = T, stringsAsFactors = F)
all_samples <- dplyr::select(table, sample = Source.Name,
                          accession = Comment.ENA_RUN.,
                          condition = FactorValue.MEDIA.)
polya_samples <- all_samples[grepl("pA$", all_samples$sample), ]
write(polya_samples$accession, file = '../metadata/accessions.txt')

write.table(polya_samples, file = '../metadata/exp_info.txt', sep = "\t",
            quote = FALSE, row.names = FALSE)
