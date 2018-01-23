.libPaths(c("~/R_library", .libPaths()))
args <- commandArgs(trailing = TRUE)

if (length(args) != 1) {
  stop('Must supply only 1 argument (the number of cores).')
}

cores <- 20

cores <- args[1]

source('benchmark_methods.R')
source('gene_common.R')

library('dplyr')
library('mamabear')
library('parallel')

options(mc.cores = cores)

transcript_gene_mapping <- get_human_gene_names()


sample_info <- readRDS('../results/finn_subsamples.rds')

sleuth_training <- mclapply(sample_info,
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
# reorder to match sample_info
obs_raw <- lapply(seq_along(sample_info),
  function(i) {
    obs <- obs_raw[[i]]
    training <- sample_info[[i]]
    obs[, rownames(training)]
  })

sleuth_filter_bool <- lapply(sleuth_training, '[[', 'sleuth_filter')
sleuth_filter_bool <- lapply(sleuth_filter_bool, function(x) {
    names(x) <- rownames(obs_raw[[1]])
    x
  })

# make the oracle results
oracle <- data.frame(target_id = names(sleuth_filter_bool[[1]]), is_de = FALSE,
  log_fc = NA, stringsAsFactors = FALSE)
oracle <- lapply(1:length(sleuth_training), function(i) oracle)

dummy_filter <- rep(TRUE, nrow(obs_raw[[1]]))
names(dummy_filter) <- rownames(obs_raw[[1]])

all_training <- list()
all_training$sleuth.lrt <- lapply(sleuth_training, '[[', 'sleuth.lrt')
all_training$sleuth.wt <- lapply(sleuth_training, '[[', 'sleuth.wt')

denom <- "ENST00000373795.6" # this had lowest COV for all of the null samples
sleuth_alr_training <- mclapply(seq_along(sample_info),
  function(i) {
    training <- sample_info[[i]]
    alr <- run_alr(training, denom = denom, num_cores = 1)
    list(sleuthALR.wt = alr$sleuthALR.wt, sleuthALR.lrt = alr$sleuthALR.lrt)
  })
all_training$sleuthALR.lrt <- lapply(alr_training, '[[', 'sleuthALR.lrt')
all_training$sleuthALR.wt <- lapply(alr_training, '[[', 'sleuthALR.wt')

all_training$ALDEx2 <- mclapply(seq_along(sample_info),
  function(i) {
    training <- sample_info[[i]]
    obs <- obs_raw[[i]]
    aldex2_filter_and_run(obs, denom = 'iqlr', training, dummy_filter, "wilcoxon")$result
  })

all_training$DESeq2 <- mclapply(seq_along(sample_info),
  function(i) {
    training <- sample_info[[i]]
    obs <- obs_raw[[i]]
    DESeq2_filter_and_run_intersect(obs, training, dummy_filter)$result
  })

all_training$edgeR <- mclapply(seq_along(sample_info),
  function(i) {
    training <- sample_info[[i]]
    obs <- obs_raw[[i]]
    edgeR_filter_and_run(obs, training, dummy_filter)$result
  })

# use the sleuth filter since limma doesn't have a recommended filter
all_training$limmaVoom <- mclapply(seq_along(sample_info),
  function(i) {
    training <- sample_info[[i]]
    obs <- obs_raw[[i]]
    sf <- sleuth_filter_bool[[i]]
    limma_filter_and_run(obs, training, sf)$result
  })

# do some renaming
all_training$sleuth <- all_training$sleuth.lrt
all_training$`sleuth-ALR` <- all_training$sleuthALR.lrt
all_training$voom <- all_training$limmaVoom

all_training$sleuth.lrt <- NULL
all_training$sleuth.wt <- NULL
all_training$sleuthALR.lrt <- NULL
all_training$sleuthALR.wt <- NULL
all_training$limmaVoom <- NULL

self_benchmark <- lapply(seq_along(all_training),
  function(i) {
    method <- names(all_training)[[i]]
    training <- all_training[[i]]
    Map(
      function(x, y) new_de_benchmark(list(x), method, y),
      training, oracle)
    })

saveRDS(self_benchmark, file = '../results/null_resampling/isoform.rds')
