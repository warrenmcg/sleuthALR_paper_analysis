.libPaths(c("~/R_library", .libPaths()))
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop('Usage: Rscript yeast_results.R N_CPU')
}

n_cpu <- args[1]


source("../../geuvadis/R/benchmark_methods.R")
source("../../geuvadis/R/gene_common.R")

library(absSimSeq)
library("parallel")
library("mamabear")
options(mc.cores = n_cpu)

annos <- '~/sleuth_paper_analysis/annotation/ASM294v2.pombase.all_annos.txt'
transcript_gene_mapping <- read.table(annos, header = T, sep = "\t",
                                      quote = "", stringsAsFactors = FALSE)

###
# new refactoring
###

all_results_kal <- list()

si <- read.table('../metadata/exp_info.txt', header = T, sep = "\t", stringsAsFactors = FALSE)
si$path <- list.files('../results/kallisto', pattern = "ERR", full.names = TRUE)
si$condition <- factor(si$condition)
levels(si$condition) <- c('A', 'B')


sir <- run_sleuth(si, gene_mode = NULL, max_bootstrap = 100, num_cores = n_cpu)
filter <- sir$so$filter_bool
counts <- sleuth:::spread_abundance_by(sir$so$obs_raw, "est_counts")
tpms <- sleuth:::spread_abundance_by(sir$so$obs_raw, "tpm")
all_results_kal$sleuth.lrt <- sir$sleuth.lrt
all_results_kal$sleuth.wt <- sir$sleuth.wt
rm(sir)

message(paste('running sleuth-ALR', Sys.time()))

denom <- readRDS('../abs_counts/denom.rds')
alr_res <- run_alr(si, max_bootstrap = 100, denom = denom, delta = 0.001,
                   num_cores = n_cpu, which_var = "obs_tpm")
all_results_kal$sleuthALR.lrt <- alr_res$sleuthALR.lrt
all_results_kal$sleuthALR.wt <- alr_res$sleuthALR.wt
rm(alr_res)

## Choose the denom with the most consistent abundance across samples
cov <- matrixStats::rowSds(tpms) / rowMeans(tpms)
cons_denom <- names(which(cov == min(cov[which(cov>0)], na.rm=T)))
## This should be SPAC1142.01
message(paste('The feature with the most consistent abundance is', cons_denom))
alr_res_cons <- run_alr(si, max_bootstrap = 100, denom = cons_denom, delta = 0.001,
                             num_cores = n_cpu, which_var = "obs_tpm")
all_results_kal$sleuthALR.cons.lrt <- alr_res_cons$sleuthALR.lrt
all_results_kal$sleuthALR.cons.wt <- alr_res_cons$sleuthALR.wt
rm(alr_res_cons)

message(paste('running sleuth-ALR counts', Sys.time()))
alr_counts_res <- run_alr(si, max_bootstrap = 100, denom = denom, delta = 0.01,
                          num_cores = n_cpu, which_var = "obs_counts")
all_results_kal$sleuthALR.counts.lrt <- alr_counts_res$sleuthALR.lrt
all_results_kal$sleuthALR.counts.wt <- alr_counts_res$sleuthALR.wt
rm(alr_counts_res)

message(paste('running ALDEx2 iqlr', Sys.time()))
aldex2 <- runALDEx2(counts, si$condition, denom = 'iqlr', FALSE, 't', 'all', filter)
all_results_kal$ALDEx2.overlap <- aldex2$overlap
all_results_kal$ALDEx2.welch <- aldex2$welch
all_results_kal$ALDEx2.wilcoxon <- aldex2$wilcoxon

message(paste('running ALDEx2 denom', Sys.time()))
aldex2 <- runALDEx2(counts, si$condition, denom = denom, FALSE, 't', 'all', filter)
all_results_kal$ALDEx2.alr.overlap <- aldex2$overlap
all_results_kal$ALDEx2.alr.welch <- aldex2$welch
all_results_kal$ALDEx2.alr.wilcoxon <- aldex2$wilcoxon

message(paste('running limma', Sys.time()))
rownames(si) <- si$sample
counts <- round(counts)
mode(counts) <- 'integer'
cds <-  make_count_data_set(counts[filter, ], si)
all_results_kal$limmaVoom <- runVoom(cds, FALSE, FALSE)

message(paste('running limma with denom', Sys.time()))
all_results_kal$limmaVoom_denom <- runVoom(cds, FALSE, FALSE, denom = denom)

message(paste('running limma with iDEGES/TCC', Sys.time()))
all_results_kal$limmaVoom_iDEGES <- run_iDEGES(cds, FALSE, FALSE, TRUE, norm_method = 'edger', test_method = 'voom')$results

message(paste('running DESeq2', Sys.time()))
all_results_kal$DESeq2 <- runDESeq2(cds, FALSE, FALSE, TRUE)$results

message(paste('running DESeq2 with RUVg', Sys.time()))
all_results_kal$DESeq2_RUVg <- runDESeq2(cds, FALSE, FALSE, TRUE, denom = denom, RUVg = TRUE)$results

message(paste('running DESeq2 with denom', Sys.time()))
all_results_kal$DESeq2_denom <- runDESeq2(cds, FALSE, FALSE, TRUE, denom = denom)$results

message(paste('running DESeq2 with iDEGES/TCC', Sys.time()))
all_results_kal$DESeq2_iDEGES <- run_iDEGES(cds, FALSE, FALSE, TRUE, norm_method = 'deseq2', test_method = 'deseq2')$results

message(paste('running edgeR', Sys.time()))
design <- NULL
all_results_kal$edgeR <- runEdgeR(cds, FALSE, FALSE, TRUE, design)

message(paste('running edgeR with RUVg', Sys.time()))
design <- NULL
all_results_kal$edgeR_RUVg <- runEdgeR(cds, FALSE, FALSE, TRUE, design, denom = denom, RUVg = TRUE)

message(paste('running edgeR with denom', Sys.time()))
design <- NULL
all_results_kal$edgeR_denom <- runEdgeR(cds, FALSE, FALSE, TRUE, design, denom = denom)

message(paste('running edgeR with iDEGES/TCC', Sys.time()))
all_results_kal$edgeR_iDEGES <- run_iDEGES(cds, FALSE, FALSE, TRUE, norm_method = 'edger', test_method = 'edger')$results

all_results_kal <- absSimSeq::rename_fc_list(all_results_kal)

saveRDS(all_results_kal, file = file.path('../results', 'yeast_results.rds'))

sig_tab <- t(sapply(all_results_kal, function(x) {
  if(ncol(x) == 3) return(c(n_up = NA, n_down = NA, n_ns = NA))
  n_up <- sum(x$qval <= 0.05 & x$log_fc > 0)
  n_down <- sum(x$qval <= 0.05 & x$log_fc < 0)
  n_ns <- sum(x$qval > 0.05)
  c(n_up = n_up, n_down = n_down, n_ns = n_ns)
}))
mod_sig_tab <- sig_tab[!grepl("lrt", rownames(sig_tab)), ]
mod_sig_tab <- mod_sig_tab[!grepl("counts", rownames(mod_sig_tab)), ]
rownames(mod_sig_tab) <- c("sleuth", "sleuth-ALR", "sleuth-ALR (trend)",
  "ALDEx2 IQLR overlap", "ALDEx2 IQLR Welch", "ALDEx2 IQLR Wilcoxon",
  "ALDEx2 + CN overlap", "ALDEx2 + CN Welch", "ALDEx2 + CN Wilcoxon",
  "limma", "limma + CN", "DESeq2", "DESeq2 + RUVg",
  "DESeq2 + CN", "edgeR", "edgeR + RUVg", "edgeR + CN")

yeast_annos <- transcript_gene_mapping
tested_ids <- yeast_annos$gene_id[which(yeast_annos$target_id %in% all_results_kal[[1]]$target_id)]
## match current annotations to the old annotations where they don't match
tested_ids[c(741, 1032, 4042, 6325, 6333, 6334)] <- c("SPBC713.13", "SPAC823.02", "SPAC1556.06.1", "SPNCRNA.98", "SPSNRNA.03", "SPSNORNA.40")

abs_counts <- read.table('../abs_counts/abs_counts.txt', header = T, sep = "\t", quote = "", stringsAsFactors=F)
filt_counts <- abs_counts[which(abs_counts$Systematic.name %in% tested_ids), ]
mod_counts <- filt_counts
mod_counts[,3:6] <- apply(mod_counts[,3:6], 2, function(x) ifelse(is.na(x), 0, x))
mod_counts$mod_fc <- log2(rowMeans(mod_counts[,5:6]) / rowMeans(mod_counts[,3:4]))

n_up <- sum(mod_counts$mod_fc > 0, na.rm = T)
n_down <- sum(mod_counts$mod_fc < 0, na.rm = T)
n_na <- sum(is.na(mod_counts$mod_fc))
abs_mat <- matrix(data = c(n_up = n_up, n_down = n_down, n_ns = n_na), nrow = 1, dimnames = list(c("Absolute Counts"), c("n_up", "n_down", "n_ns")))

full_sig_tab <- rbind(mod_sig_tab, abs_mat)
full_sig_df <- data.frame(method = rownames(full_sig_tab), full_sig_tab)
write.table(full_sig_df, file = "../results/table1.txt", sep = "\t", quote = F, row.names = F)
