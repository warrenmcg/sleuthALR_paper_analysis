.libPaths(c("~/R_library", .libPaths()))
args <- commandArgs(trailing = TRUE)

if (length(args) != 1) {
  stop('Must supply only 1 argument (the number of cores).')
}

cores <- 20

cores <- args[1]

# This file creates a reference for each method on the complete data set
# at the GENE level.

source('../../geuvadis/R/benchmark_methods.R')
source('../../geuvadis/R/gene_common.R')

library('dplyr')
library('mamabear')
library('parallel')

options(mc.cores = cores)

# load all of the metadata
# training_sets <- readRDS('../metadata/training_sets.rds')
# validation_sets <- readRDS('../metadata/validation_sets.rds')
source('get_metadata.R')

transcript_gene_mapping <- get_mouse_gene_names()

# this is a temporary fix
transcript_gene_mapping <- dplyr::mutate(transcript_gene_mapping,
  target_id_backup = target_id)
transcript_gene_mapping <- dplyr::mutate(transcript_gene_mapping,
  target_id = ensembl_transcript_id)

###
# loading featureCounts data
###
training_sets <- lapply(training_sets,
  function(df) {
    df$featureCounts <- file.path('..', 'results', 'single', df$sample,
      'featureCounts.txt')
    df
  })
validation_sets <- lapply(validation_sets,
  function(df) {
    df$featureCounts <- file.path('..', 'results', 'single', df$sample,
      'featureCounts.txt')
    df
  })

# restrict to only the genes that are in the ensembl database
gene_names <- unique(transcript_gene_mapping$ens_gene)

dummy_filter <- rep(TRUE, length(gene_names))
names(dummy_filter) <- gene_names

training_counts <- lapply(training_sets,
  function(df) {
    obs <- load_union_counts_general(df$featureCounts, df$sample)
    current_filter <- intersect(rownames(obs), names(dummy_filter))
    obs[current_filter, ]
  })
validation_counts <- lapply(validation_sets,
  function(df) {
    obs <- load_union_counts_general(df$featureCounts, df$sample)
    current_filter <- intersect(rownames(obs), names(dummy_filter))
    obs[current_filter, ]
  })

# dummy_filter now contains the intersection between sleuth and featureCounts
dummy_filter <- rep(TRUE, nrow(validation_counts[[1]]))
names(dummy_filter) <- rownames(validation_counts[[1]])
dummy_filter_df <- data.frame(target_id = names(dummy_filter),
  stringsAsFactors = FALSE)

###
# run the gene methods
###

all_validation <- list()

message(paste("running sleuth", Sys.time()))
sleuth_validation <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    sir <- run_sleuth(validation, gene_mode = 'aggregate', gene_column = 'ens_gene')
    list(sleuth.lrt = sir$sleuth.lrt, sleuth.wt = sir$sleuth.wt)
  })
all_validation$sleuth.lrt <- lapply(sleuth_validation, '[[', 'sleuth.lrt')
all_validation$sleuth.wt <- lapply(sleuth_validation, '[[', 'sleuth.wt')

message(paste("running sleuth-ALR", Sys.time()))
denom <- "ENSMUSG00000038014.7"
alr_res <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    run_alr(validation, denom = denom, num_cores = 1, gene_column = 'ens_gene')
    list(sleuthALR.lrt = sir$sleuth.lrt, sleuthALR.wt = sir$sleuth.wt)
  })
all_validation$sleuthALR.lrt <- lapply(alr_res, '[[', 'sleuthALR.lrt')
all_validation$sleuthALR.wt <- lapply(alr_res, '[[', 'sleuthALR.wt')

message(paste("running DESeq2", Sys.time()))
all_validation$DESeq2 <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    obs <- validation_counts[[i]]
    DESeq2_filter_and_run_intersect(obs, validation, dummy_filter)$result
  })

message(paste("running edgeR", Sys.time()))
edgeR_results <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    obs <- validation_counts[[i]]
    edgeR_filter_and_run(obs, validation, dummy_filter)
  })
edgeR_filter_validation <- lapply(edgeR_results, '[[', 'filter')
all_validation$edgeR <- lapply(edgeR_results, '[[', 'result')

# use edgeR as the filter for limmaVoom and EBSeq
edgeR_filter_validation <- lapply(edgeR_filter_validation,
  function(x) {
    y <- dummy_filter
    y <- !y
    y[x] <- TRUE
    y
  })

message(paste("running limma", Sys.time()))
all_validation$limmaVoom <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    obs <- validation_counts[[i]]
    current_filter <- edgeR_filter_validation[[i]]
    limma_filter_and_run(obs, validation, current_filter)$result
  })

message(paste("running ALDEx2-filtered with IQLR", Sys.time()))
all_validation$ALDEx2.filt <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    obs <- validation_counts[[i]]
    current_filter <- edgeR_filter_validation[[i]]
    aldex2_filter_and_run(obs, 'iqlr', validation, current_filter, "all")
  })
all_validation$ALDEx2.filt.welch <- lapply(all_validation$ALDEx2.filt, '[[', 'ALDEx2.welch')
all_validation$ALDEx2.filt.wilcoxon <- lapply(all_validation$ALDEx2.filt, '[[', 'ALDEx2.wilcoxon')
all_validation$ALDEx2.filt.overlap <- lapply(all_validation$ALDEx2.filt, '[[', 'ALDEx2.overlap')
all_validation$ALDEx2.filt <- NULL

message(paste('Adding truth column', Sys.time()))
all_validation <- lapply(all_validation,
  function(validation) {
    lapply(validation,
      function(x) {
        dplyr::mutate(x, is_de = ifelse(qval <= 0.05, TRUE, FALSE), log_fc = NA)
      })
  })

###
# run on the training set
###

all_training <- list()

message(paste("running sleuth on training sets", Sys.time()))
sleuth_training <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    sir <- run_sleuth(training, gene_mode = 'aggregate', gene_column = 'ens_gene')
    list(sleuth.lrt = sir$sleuth.lrt, sleuth.wt = sir$sleuth.wt)
  })
all_training$sleuth.lrt <- lapply(sleuth_training, '[[', 'sleuth.lrt')
all_training$sleuth.wt <- lapply(sleuth_training, '[[', 'sleuth.wt')

message(paste("running sleuth-ALR on training sets", Sys.time()))
alr_res <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    run_alr(training, denom = denom, num_cores = 1, gene_column = 'ens_gene')
    list(sleuthALR.lrt = sir$sleuth.lrt, sleuthALR.wt = sir$sleuth.wt)
  })
all_training$sleuthALR.lrt <- lapply(alr_res, '[[', 'sleuthALR.lrt')
all_training$sleuthALR.wt <- lapply(alr_res, '[[', 'sleuthALR.wt')

message(paste("running DESeq2 on training sets", Sys.time()))
all_training$DESeq2 <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    obs <- training_counts[[i]]
    DESeq2_filter_and_run_intersect(obs, training, dummy_filter)$result
  })

message(paste("running edgeR on training sets", Sys.time()))
edgeR_training <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    obs <- training_counts[[i]]
    edgeR_filter_and_run(obs, training, dummy_filter)
  })
edgeR_filter_training <- lapply(edgeR_training, '[[', 'filter')
all_training$edgeR <- lapply(edgeR_training, '[[', 'result')

# use edgeR as the filter for limmaVoom and EBSeq
edgeR_filter_training <- lapply(edgeR_filter_training,
  function(x) {
    y <- dummy_filter
    y <- !y
    y[x] <- TRUE
    y
  })

message(paste("running limma on training sets", Sys.time()))
all_training$limmaVoom <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    obs <- training_counts[[i]]
    current_filter <- edgeR_filter_training[[i]]
    limma_filter_and_run(obs, training, current_filter)$result
  })

message(paste("running ALDEx2-filtered with IQLR on training sets", Sys.time()))
all_training$ALDEx2.filt <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    obs <- training_counts[[i]]
    current_filter <- edgeR_filter_training[[i]]
    aldex2_filter_and_run(obs, 'iqlr', training, current_filter, "all")
  })
all_training$ALDEx2.filt.welch <- lapply(all_training$ALDEx2.filt, '[[', 'ALDEx2.welch')
all_training$ALDEx2.filt.wilcoxon <- lapply(all_training$ALDEx2.filt, '[[', 'ALDEx2.wilcoxon')
all_training$ALDEx2.filt.overlap <- lapply(all_training$ALDEx2.filt, '[[', 'ALDEx2.overlap')
all_training$ALDEx2.filt <- NULL

message(paste("combining validation and training sets", Sys.time()))
self_benchmark <- lapply(seq_along(all_training),
  function(i) {
    method <- names(all_training)[[i]]
    training <- all_training[[i]]
    validation <- all_validation[[method]]
    Map(
      function(x, y) new_de_benchmark(list(x), method, y),
      training, validation)
  })
message("saving the benchmarks")
names(self_benchmark) <- names(all_training)

self_fdr <- lapply(self_benchmark, average_sensitivity_specificity)

saveRDS(self_benchmark, '../results/gene_self_benchmark.rds')
saveRDS(self_fdr, '../results/gene_self_fdr.rds')
