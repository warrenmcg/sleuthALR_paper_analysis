args <- commandArgs(trailing = TRUE)

if (length(args) != 1) {
  stop('Must supply only 1 argument (the number of cores).')
}

cores <- args[1]

# This file creates a reference for each method on the complete data set
# at the ISOFORM level.

source('../../geuvadis/R/benchmark_methods.R')
source('../../geuvadis/R/gene_common.R')

library('dplyr')
library('mamabear')
library('parallel')

options(mc.cores = cores)

# load all of the metadata
source('get_metadata.R')


###
# run each method on the validation sets
###

# get the raw data from sleuth
sleuth_validation <- mclapply(validation_sets,
  function(validation) {
    sir <- run_sleuth(validation, gene_mode = NULL)
    sleuth_filter <- sir$so$filter_bool
    obs_raw <- sleuth:::spread_abundance_by(sir$so$obs_raw, "est_counts")
    dummy_filter <- rep(TRUE, nrow(obs_raw))
    names(dummy_filter) <- rownames(obs_raw)

    list(sleuth.lrt = sir$sleuth.lrt, sleuth.wt = sir$sleuth.wt,
      obs_raw = obs_raw, sleuth_filter = sleuth_filter)
  })

obs_raw <- lapply(sleuth_validation, '[[', 'obs_raw')
# reorder to match validation_sets
obs_raw <- lapply(seq_along(validation_sets),
  function(i) {
    obs <- obs_raw[[i]]
    validation <- validation_sets[[i]]
    obs[, rownames(validation)]
  })

sleuth_filter_bool <- lapply(sleuth_validation, '[[', 'sleuth_filter')
sleuth_filter_bool <- lapply(sleuth_filter_bool, function(x) {
    names(x) <- rownames(obs_raw[[1]])
    x
  })

dummy_filter <- rep(TRUE, nrow(obs_raw[[1]]))
names(dummy_filter) <- rownames(obs_raw[[1]])

dummy_filter_df <- data.frame(target_id = names(dummy_filter),
  stringsAsFactors = FALSE)
dummy_filter_df$target_id <- sub('\\.[0-9]+', '', dummy_filter_df$target_id)

all_validation <- list()
all_validation$sleuth.lrt <- lapply(sleuth_validation, '[[', 'sleuth.lrt')
all_validation$sleuth.wt <- lapply(sleuth_validation, '[[', 'sleuth.wt')

rm(sleuth_validation)

all_validation$DESeq2 <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    obs <- obs_raw[[i]]
    DESeq2_filter_and_run_intersect(obs, validation, dummy_filter)$result
  })

all_validation$edgeR <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    obs <- obs_raw[[i]]
    edgeR_filter_and_run(obs, validation, dummy_filter)$result
  })

# use the sleuth filter since limma doesn't have a recommended filter
all_validation$limmaVoom <- mclapply(seq_along(validation_sets),
  function(i) {
    validation <- validation_sets[[i]]
    obs <- obs_raw[[i]]
    sf <- sleuth_filter_bool[[i]]
    limma_filter_and_run(obs, validation, sf)$result
  })

all_validation <- lapply(all_validation,
  function(validation) {
    lapply(validation,
      function(x) {
        dplyr::mutate(x, is_de = ifelse(qval <= 0.05, TRUE, FALSE), log_fc = NA)
      })
  })

###
# run each program as you would in practice
###

sleuth_training <- mclapply(training_sets,
  function(training) {
    sir <- run_sleuth(training, gene_mode = NULL)
    sleuth_filter <- sir$so$filter_bool
    obs_raw <- sleuth:::spread_abundance_by(sir$so$obs_raw, "est_counts")
    dummy_filter <- rep(TRUE, nrow(obs_raw))
    names(dummy_filter) <- rownames(obs_raw)

    list(sleuth.lrt = sir$sleuth.lrt, sleuth.wt = sir$sleuth.wt,
      obs_raw = obs_raw, sleuth_filter = sleuth_filter)
  })

obs_raw <- lapply(sleuth_training, '[[', 'obs_raw')
# reorder to match training_sets
obs_raw <- lapply(seq_along(training_sets),
  function(i) {
    obs <- obs_raw[[i]]
    training <- training_sets[[i]]
    obs[, rownames(training)]
  })

sleuth_filter_bool <- lapply(sleuth_training, '[[', 'sleuth_filter')
sleuth_filter_bool <- lapply(sleuth_filter_bool, function(x) {
    names(x) <- rownames(obs_raw[[1]])
    x
  })

dummy_filter <- rep(TRUE, nrow(obs_raw[[1]]))
names(dummy_filter) <- rownames(obs_raw[[1]])

all_training <- list()
all_training$sleuth.lrt <- lapply(sleuth_training, '[[', 'sleuth.lrt')
all_training$sleuth.wt <- lapply(sleuth_training, '[[', 'sleuth.wt')

all_training$DESeq2 <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    obs <- obs_raw[[i]]
    DESeq2_filter_and_run_intersect(obs, training, dummy_filter)$result
  })

all_training$edgeR <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    obs <- obs_raw[[i]]
    edgeR_filter_and_run(obs, training, dummy_filter)$result
  })

# use the sleuth filter since limma doesn't have a recommended filter
all_training$limmaVoom <- mclapply(seq_along(training_sets),
  function(i) {
    training <- training_sets[[i]]
    obs <- obs_raw[[i]]
    sf <- sleuth_filter_bool[[i]]
    limma_filter_and_run(obs, training, sf)$result
  })

self_benchmark <- lapply(seq_along(all_training),
  function(i) {
    method <- names(all_training)[[i]]
    training <- all_training[[i]]
    validation <- all_validation[[method]]
    Map(
      function(x, y) new_de_benchmark(list(x), method, y),
      training, validation)
  })
names(self_benchmark) <- names(all_training)

self_fdr <- lapply(self_benchmark, average_sensitivity_specificity)

saveRDS(self_benchmark, '../results/isoform_self_benchmark.rds')
saveRDS(self_fdr, '../results/isoform_self_fdr.rds')
