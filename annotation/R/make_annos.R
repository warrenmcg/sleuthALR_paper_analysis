system('grep "^>" ../ASM294v2.pombase.all.fa | cut -d">" -f2 > ../headers.txt')
lines <- read.table("../headers.txt", sep = "\t", quote = "", stringsAsFactors = FALSE)$V1
splits <- strsplit(lines, split = " ", fixed = T)
data <- lapply(splits, function(x) {
   target_id <- x[1]
   gene_id <- strsplit(x[4], ":", fixed = T)[[1]][2]
   gene_symbol <- strsplit(x[7], ":", fixed = T)[[1]][2]
   gene_biotype <- strsplit(x[5], ":", fixed = T)[[1]][2]
   transcript_biotype <- strsplit(x[6], ":", fixed = T)[[1]][2]
   if (length(x) < 8)
     description <- NA
   else {
     description <- paste(x[8:length(x)], collapse = " ")
     description <- gsub("description:", "", description)
   }
   c(target_id = target_id, gene_id = gene_id, gene_symbol = gene_symbol, gene_biotype = gene_biotype, transcript_biotype = transcript_biotype, full_name = description)
})

new_data <- do.call(rbind, data)
new_data <- as.data.frame(new_data)
write.table(new_data, "../ASM294v2.pombase.all_annos.txt", sep = "\t", row.names = F, quote = F)
