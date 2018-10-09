.libPaths(c("~/R_library", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop('Usage: Rscript bench_isoform.R N_CPU SIM_NAME')
}

n_cpu <- args[1]
sim_name <- args[2]

# temporary for debugging
# n_cpu <- 20
# sim_name <- 'isoform_5_5_15_387_1'

source("benchmark_methods.R")
source("gene_common.R")

library("parallel")
library("mamabear")
library("absSimSeq")
options(mc.cores = n_cpu)

transcript_gene_mapping <- get_human_gene_names()

###
# new refactoring
###

all_results_kal <- list()

N_SIM <- 15
GROUP_BREAKS <- c(5, 10)

sim <- parse_simulation(sim_name)

sample_info <- lapply(1:N_SIM,
  function(i) {
    n <- sim$a + sim$b
    num <- sprintf('%02d', i)
    sample_nums <- sprintf('%02d', 1:n)
    kal_dirs <- file.path('..', 'sims', sim_name, paste0("run", num),
      paste0("sample_", sample_nums), "kallisto")
    get_sample_to_condition(sim$a, sim$b, kal_dirs)
  })

initial_sleuth_res <- mclapply(1:N_SIM,
  function(i) {
    si <- sample_info[[i]]
    so <- sleuth_prep(si, max_bootstrap = 2, num_cores = 1, normalize = FALSE)
    counts <- sleuth::sleuth_to_matrix(so, "obs_raw", "est_counts")
    filter <- so$filter_bool
    list(counts = counts, filter = filter)
  })
obs_counts <- lapply(initial_sleuth_res, '[[', 'counts')
sleuth_filters <- lapply(initial_sleuth_res, '[[', 'filter')

message(paste('running sleuth', Sys.time()))
sleuth_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    filter <- sleuth_filters[[i]]
    filter_ids <- names(filter)[filter]
    sir <- run_sleuth(si, gene_mode = NULL, max_bootstrap = 100, num_cores = 1,
                      filter_target_id = filter_ids)
    list(sleuth.lrt = sir$sleuth.lrt, sleuth.wt = sir$sleuth.wt)
  })
all_results_kal$sleuth.lrt <- lapply(sleuth_res, '[[', 'sleuth.lrt')
all_results_kal$sleuth.wt <- lapply(sleuth_res, '[[', 'sleuth.wt')
rm(sleuth_res)

message(paste('running sleuth-ALR', Sys.time()))
data('ERCC92_data', package = "absSimSeq")
log_means <- log2(rowMeans(ERCC92_data[,c(4:5)]))
denoms <- names(log_means[log_means > 3])
denoms <- gsub('_', '-', denoms)

alr_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    filter <- sleuth_filters[[i]]
    filter_ids <- names(filter)[filter]
    alr <- run_alr(si, max_bootstrap = 100, denom = denoms,
                   num_cores = 1, which_var = "obs_tpm", filter_target_id = filter_ids)
    message(paste0("the denominator for sleuth-ALR run ", i, ": ", paste(alr$denoms, collapse = ", ")))
    list(sleuthALR.lrt = alr$sleuthALR.lrt, sleuthALR.wt = alr$sleuthALR.wt)
  })
all_results_kal$sleuthALR.lrt <- lapply(alr_res, '[[', 'sleuthALR.lrt')
all_results_kal$sleuthALR.wt <- lapply(alr_res, '[[', 'sleuthALR.wt')
rm(alr_res)

message(paste('running sleuth-ALR counts', Sys.time()))
alr_counts_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    filter <- sleuth_filters[[i]]
    filter_ids <- names(filter)[filter]
    alr <- run_alr(si, max_bootstrap = 100, denom = denoms,
                   num_cores = 1, which_var = "obs_counts", filter_target_id = filter_ids)
    message(paste0("the denominator for sleuth-ALR run ", i, ": ", paste(alr$denoms, collapse = ", ")))
    list(sleuthALR.lrt = alr$sleuthALR.lrt, sleuthALR.wt = alr$sleuthALR.wt)
  })
all_results_kal$sleuthALR.counts.lrt <- lapply(alr_counts_res, '[[', 'sleuthALR.lrt')
all_results_kal$sleuthALR.counts.wt <- lapply(alr_counts_res, '[[', 'sleuthALR.wt')
rm(alr_counts_res)

message(paste('running sleuth-ALR with delta = 0.1', Sys.time()))
alr_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    filter <- sleuth_filters[[i]]
    filter_ids <- names(filter)[filter]
    alr <- run_alr(si, max_bootstrap = 100, denom = denoms, delta = 0.1,
                   num_cores = 1, which_var = "obs_tpm", filter_target_id = filter_ids)
    message(paste0("the denominator for sleuth-ALR run ", i, ": ", paste(alr$denoms, collapse = ", ")))
    list(sleuthALR.lrt = alr$sleuthALR.lrt, sleuthALR.wt = alr$sleuthALR.wt)
  })
all_results_kal$sleuthALR.0.1.lrt <- lapply(alr_res, '[[', 'sleuthALR.lrt')
all_results_kal$sleuthALR.0.1.wt <- lapply(alr_res, '[[', 'sleuthALR.wt')
rm(alr_res)

message(paste('running sleuth-ALR counts with delta = 0.001', Sys.time()))
alr_counts_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    filter <- sleuth_filters[[i]]
    filter_ids <- names(filter)[filter]
    alr <- run_alr(si, max_bootstrap = 100, denom = denoms, delta = 0.001,
                   num_cores = 1, which_var = "obs_counts", filter_target_id = filter_ids)
    message(paste0("the denominator for sleuth-ALR run ", i, ": ", paste(alr$denoms, collapse = ", ")))
    list(sleuthALR.lrt = alr$sleuthALR.lrt, sleuthALR.wt = alr$sleuthALR.wt)
  })
all_results_kal$sleuthALR.counts.0.001.lrt <- lapply(alr_counts_res, '[[', 'sleuthALR.lrt')
all_results_kal$sleuthALR.counts.0.001.wt <- lapply(alr_counts_res, '[[', 'sleuthALR.wt')
rm(alr_counts_res)

message(paste('running sleuth-ALR TPM with delta = 0.01', Sys.time()))
alr_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    filter <- sleuth_filters[[i]]
    filter_ids <- names(filter)[filter]
    alr <- run_alr(si, max_bootstrap = 100, denom = denoms, delta = 0.01,
                   num_cores = 1, which_var = "obs_tpm", filter_target_id = filter_ids)
    message(paste0("the denominator for sleuth-ALR run ", i, ": ", paste(alr$denoms, collapse = ", ")))
    list(sleuthALR.lrt = alr$sleuthALR.lrt, sleuthALR.wt = alr$sleuthALR.wt)
  })
all_results_kal$sleuthALR.0.01.lrt <- lapply(alr_res, '[[', 'sleuthALR.lrt')
all_results_kal$sleuthALR.0.01.wt <- lapply(alr_res, '[[', 'sleuthALR.wt')
rm(alr_res)

message(paste('running sleuth-ALR counts with delta = 0.01', Sys.time()))
alr_counts_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    filter <- sleuth_filters[[i]]
    filter_ids <- names(filter)[filter]
    alr <- run_alr(si, max_bootstrap = 100, denom = denoms, delta = 0.01,
                   num_cores = 1, which_var = "obs_counts", filter_target_id = filter_ids)
    message(paste0("the denominator for sleuth-ALR run ", i, ": ", paste(alr$denoms, collapse = ", ")))
    list(sleuthALR.lrt = alr$sleuthALR.lrt, sleuthALR.wt = alr$sleuthALR.wt)
  })
all_results_kal$sleuthALR.counts.0.01.lrt <- lapply(alr_counts_res, '[[', 'sleuthALR.lrt')
all_results_kal$sleuthALR.counts.0.01.wt <- lapply(alr_counts_res, '[[', 'sleuthALR.wt')
rm(alr_counts_res)

message(paste('running sleuth-ALR TPM with delta = 0.001', Sys.time()))
alr_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    filter <- sleuth_filters[[i]]
    filter_ids <- names(filter)[filter]
    alr <- run_alr(si, max_bootstrap = 100, denom = denoms, delta = 0.001,
                   num_cores = 1, which_var = "obs_tpm", filter_target_id = filter_ids)
    message(paste0("the denominator for sleuth-ALR run ", i, ": ", paste(alr$denoms, collapse = ", ")))
    list(sleuthALR.lrt = alr$sleuthALR.lrt, sleuthALR.wt = alr$sleuthALR.wt)
  })
all_results_kal$sleuthALR.0.001.lrt <- lapply(alr_res, '[[', 'sleuthALR.lrt')
all_results_kal$sleuthALR.0.001.wt <- lapply(alr_res, '[[', 'sleuthALR.wt')
rm(alr_res)

message(paste('running sleuth-ALR counts with delta = 0.1', Sys.time()))
alr_counts_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    filter <- sleuth_filters[[i]]
    filter_ids <- names(filter)[filter]
    alr <- run_alr(si, max_bootstrap = 100, denom = denoms, delta = 0.1,
                   num_cores = 1, which_var = "obs_counts", filter_target_id = filter_ids)
    message(paste0("the denominator for sleuth-ALR run ", i, ": ", paste(alr$denoms, collapse = ", ")))
    list(sleuthALR.lrt = alr$sleuthALR.lrt, sleuthALR.wt = alr$sleuthALR.wt)
  })
all_results_kal$sleuthALR.counts.0.1.lrt <- lapply(alr_counts_res, '[[', 'sleuthALR.lrt')
all_results_kal$sleuthALR.counts.0.1.wt <- lapply(alr_counts_res, '[[', 'sleuthALR.wt')
rm(alr_counts_res)

message(paste('running sleuth-ALR TPM with delta = 1e-4', Sys.time()))
alr_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    filter <- sleuth_filters[[i]]
    filter_ids <- names(filter)[filter]
    alr <- run_alr(si, max_bootstrap = 100, denom = denoms, delta = 1E-4,
                   num_cores = 1, which_var = "obs_tpm", filter_target_id = filter_ids)
    message(paste0("the denominator for sleuth-ALR run ", i, ": ", paste(alr$denoms, collapse = ", ")))
    list(sleuthALR.lrt = alr$sleuthALR.lrt, sleuthALR.wt = alr$sleuthALR.wt)
  })
all_results_kal$sleuthALR.1E4.lrt <- lapply(alr_res, '[[', 'sleuthALR.lrt')
all_results_kal$sleuthALR.1E4.wt <- lapply(alr_res, '[[', 'sleuthALR.wt')
rm(alr_res)

message(paste('running sleuth-ALR counts with delta = 0.5', Sys.time()))
alr_counts_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    filter <- sleuth_filters[[i]]
    filter_ids <- names(filter)[filter]
    alr <- run_alr(si, max_bootstrap = 100, denom = denoms, delta = 0.5,
                   num_cores = 1, which_var = "obs_counts", filter_target_id = filter_ids)
    message(paste0("the denominator for sleuth-ALR run ", i, ": ", paste(alr$denoms, collapse = ", ")))
    list(sleuthALR.lrt = alr$sleuthALR.lrt, sleuthALR.wt = alr$sleuthALR.wt)
  })
all_results_kal$sleuthALR.counts.0.5.lrt <- lapply(alr_counts_res, '[[', 'sleuthALR.lrt')
all_results_kal$sleuthALR.counts.0.5.wt <- lapply(alr_counts_res, '[[', 'sleuthALR.wt')
rm(alr_counts_res)

message(paste('getting benchmarks', Sys.time()))
saveRDS(all_results_kal, file = file.path('../results', sim_name, 'all_kal_filt_impute_results.rds'))
all_results_kal <- absSimSeq::rename_fc_list(all_results_kal)
saveRDS(all_results_kal, file = file.path('../results', sim_name, 'all_kal_fc_filt_impute_results.rds'))
oracles <- readRDS('../sims/polyester_oracles.rds')
method_names <- names(all_results_kal)
all_kal_benchmarks <- mclapply(seq_along(all_results_kal[[1]]),
  function(i) {
    results <- lapply(all_results_kal, function(x) x[[i]])
    oracle <- oracles$alr_oracle[[i]]
    filter <- sleuth_filters[[i]]
    filter_ids <- names(filter[filter])
    filt_oracle <- oracle[which(oracle$target_id %in% filter_ids), ]
    fold_change <- FALSE
    new_de_benchmark(results, method_names, filt_oracle,
                     de_linetypes = rep(1, length(method_names)), fold_change = fold_change)
  }, mc.cores = n_cpu)

small_kal_benchmarks <- all_kal_benchmarks[1:GROUP_BREAKS[1]]
down_kal_benchmarks <- all_kal_benchmarks[(GROUP_BREAKS[1]+1):(GROUP_BREAKS[2])]
up_kal_benchmarks <- all_kal_benchmarks[(GROUP_BREAKS[2]+1):N_SIM]
rm(all_kal_benchmarks)
suppressMessages(gc())

dir.create(file.path('../results', sim_name), showWarnings = F)
saveRDS(small_kal_benchmarks, file = file.path('../results', sim_name,
  'small_kal_filt_impute_benchmarks.rds'))
saveRDS(down_kal_benchmarks, file = file.path('../results', sim_name,
  'down_kal_filt_impute_benchmarks.rds'))
saveRDS(up_kal_benchmarks, file = file.path('../results', sim_name,
  'up_kal_filt_impute_benchmarks.rds'))
