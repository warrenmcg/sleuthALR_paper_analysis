.libPaths(c("~/R_library", .libPaths()))
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop('Usage: Rscript bench_isoform.R N_CPU SIM_NAME')
}

# temporary for debugging
# n_cpu <- 20
# sim_name <- 'isoform_3_3_20_1_1'
sim_name <- 'isoform_5_5_15_387_1'

n_cpu <- args[1]
sim_name <- args[2]


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

message(paste('running sleuth', Sys.time()))
sleuth_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    so <- run_sleuth_prep(si, num_cores = 1, normalize = FALSE)
    sir <-  list(so = so)
    filter <- sir$so$filter_bool
    counts <- sleuth::sleuth_to_matrix(sir$so, "obs_raw", "est_counts")
    list(sleuth.lrt = sir$sleuth.lrt, sleuth.wt = sir$sleuth.wt, filter = filter,
         counts = counts)
  }, mc.cores = n_cpu)
all_results_kal$sleuth.lrt <- lapply(sleuth_res, '[[', 'sleuth.lrt')
all_results_kal$sleuth.wt <- lapply(sleuth_res, '[[', 'sleuth.wt')
sleuth_filters <- lapply(sleuth_res, '[[', 'filter')
obs_counts <- lapply(sleuth_res, '[[', 'counts')
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
    alr <- run_alr(si, max_bootstrap = 100, denom = denoms, delta = 0.01,
                   num_cores = 1, which_var = "obs_tpm")
    list(sleuthALR.lrt = alr$sleuthALR.lrt, sleuthALR.wt = alr$sleuthALR.wt)
  }, mc.cores = n_cpu)
all_results_kal$sleuthALR.lrt <- lapply(alr_res, '[[', 'sleuthALR.lrt')
all_results_kal$sleuthALR.wt <- lapply(alr_res, '[[', 'sleuthALR.wt')
rm(alr_res)

message(paste('running ALDEx2 iqlr', Sys.time()))
aldex2_iqlr_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    counts <- obs_counts[[i]]
    filter <- sleuth_filters[[i]]
    aldex2 <- runALDEx2(counts, si$condition, denom = 'iqlr', FALSE, 't', 'all', filter)
    list(ALDEx2.overlap = aldex2$overlap,
         ALDEx2.welch = aldex2$welch,
         ALDEx2.wilcoxon = aldex2$wilcoxon)
  }, mc.cores = 5)
all_results_kal$ALDEx2.iqlr.overlap <- lapply(aldex2_iqlr_res, '[[', 'ALDEx2.overlap')
all_results_kal$ALDEx2.iqlr.welch <- lapply(aldex2_iqlr_res, '[[', 'ALDEx2.welch')
all_results_kal$ALDEx2.iqlr.wilcoxon <- lapply(aldex2_iqlr_res, '[[', 'ALDEx2.wilcoxon')
rm(aldex2_iqlr_res)

message(paste('running ALDEx2 clr', Sys.time()))
aldex2_clr_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    counts <- obs_counts[[i]]
    filter <- sleuth_filters[[i]]
    aldex2 <- runALDEx2(counts, si$condition, denom = 'all', FALSE, 't', 'all', filter)
    list(ALDEx2.overlap = aldex2$overlap,
         ALDEx2.welch = aldex2$welch,
         ALDEx2.wilcoxon = aldex2$wilcoxon)
  }, mc.cores = 5)
all_results_kal$ALDEx2.clr.overlap <- lapply(aldex2_clr_res, '[[', 'ALDEx2.overlap')
all_results_kal$ALDEx2.clr.welch <- lapply(aldex2_clr_res, '[[', 'ALDEx2.welch')
all_results_kal$ALDEx2.clr.wilcoxon <- lapply(aldex2_clr_res, '[[', 'ALDEx2.wilcoxon')
rm(aldex2_clr_res)

message(paste('running ALDEx2 with spike-ins', Sys.time()))
aldex2_denom_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    si <- sample_info[[i]]
    counts <- obs_counts[[i]]
    filter <- sleuth_filters[[i]]
    aldex2 <- runALDEx2(counts, si$condition, denom = denoms, FALSE, 't', 'all', filter)
    list(ALDEx2.overlap = aldex2$overlap,
         ALDEx2.welch = aldex2$welch,
         ALDEx2.wilcoxon = aldex2$wilcoxon)
  }, mc.cores = 5)
all_results_kal$ALDEx2.denom.overlap <- lapply(aldex2_denom_res, '[[', 'ALDEx2.overlap')
all_results_kal$ALDEx2.denom.welch <- lapply(aldex2_denom_res, '[[', 'ALDEx2.welch')
all_results_kal$ALDEx2.denom.wilcoxon <- lapply(aldex2_denom_res, '[[', 'ALDEx2.wilcoxon')
rm(aldex2_denom_res)

message(paste('running limma', Sys.time()))
all_results_kal$limmaVoom <- mclapply(1:N_SIM,
  function(i) {
    cat('Sample: ', i, '\n')
    si <- sample_info[[i]]
    counts <- obs_counts[[i]]
    counts <- round(counts)
    mode(counts) <- 'integer'
    filter <- sleuth_filters[[i]]
    filt_counts <- counts[filter, ]
    cds <-  make_count_data_set(counts[filter, ], si)
    runVoom(cds, FALSE, FALSE)
  }, mc.cores = n_cpu)

message(paste('running limma with spike-ins', Sys.time()))
all_results_kal$limmaVoom_denom <- mclapply(1:N_SIM,
  function(i) {
    cat('Sample: ', i, '\n')
    si <- sample_info[[i]]
    counts <- obs_counts[[i]]
    counts <- round(counts)
    mode(counts) <- 'integer'
    filter <- sleuth_filters[[i]]
    filt_counts <- counts[filter, ]
    cds <-  make_count_data_set(counts[filter, ], si)
    runVoom(cds, FALSE, FALSE, denom = denoms)
  }, mc.cores = n_cpu)

message(paste('running DESeq2', Sys.time()))
all_results_kal$DESeq2 <- mclapply(1:N_SIM,
  function(i) {
    cat('Sample: ', i, '\n')
    si <- sample_info[[i]]
    counts <- obs_counts[[i]]
    counts <- round(counts)
    mode(counts) <- 'integer'
    filter <- sleuth_filters[[i]]
    cds <- make_count_data_set(counts[filter, ], si)
    runDESeq2(cds, FALSE, FALSE, TRUE)$results
  }, mc.cores = n_cpu)

message(paste('running DESeq2 with RUVg', Sys.time()))
all_results_kal$DESeq2_RUVg <- mclapply(1:N_SIM,
  function(i) {
    cat('Sample: ', i, '\n')
    si <- sample_info[[i]]
    counts <- obs_counts[[i]]
    counts <- round(counts)
    mode(counts) <- 'integer'
    filter <- sleuth_filters[[i]]
    denom <- denoms
    cds <- make_count_data_set(counts[filter, ], si)
    runDESeq2(cds, FALSE, FALSE, TRUE, denom = denom, RUVg = TRUE)$results
  }, mc.cores = n_cpu)

message(paste('running DESeq2 with spike-ins', Sys.time()))
all_results_kal$DESeq2_denom <- mclapply(1:N_SIM,
  function(i) {
    cat('Sample: ', i, '\n')
    si <- sample_info[[i]]
    counts <- obs_counts[[i]]
    counts <- round(counts)
    mode(counts) <- 'integer'
    filter <- sleuth_filters[[i]]
    denom <- denoms
    cds <- make_count_data_set(counts[filter, ], si)
    runDESeq2(cds, FALSE, FALSE, TRUE, denom = denom)$results
  }, mc.cores = n_cpu)

message(paste('running edgeR', Sys.time()))
all_results_kal$edgeR <- mclapply(1:N_SIM,
  function(i) {
    cat('Sample: ', i, '\n')
    si <- sample_info[[i]]
    counts <- obs_counts[[i]]
    counts <- round(counts)
    mode(counts) <- 'integer'
    filter <- sleuth_filters[[i]]
    cds <- make_count_data_set(counts[filter, ], si)
    design <- NULL
    runEdgeR(cds, FALSE, FALSE, TRUE, design)
  }, mc.cores = n_cpu)

message(paste('running edgeR with RUVg', Sys.time()))
all_results_kal$edgeR_RUVg <- mclapply(1:N_SIM,
  function(i) {
    cat('Sample: ', i, '\n')
    si <- sample_info[[i]]
    counts <- obs_counts[[i]]
    counts <- round(counts)
    mode(counts) <- 'integer'
    filter <- sleuth_filters[[i]]
    cds <- make_count_data_set(counts[filter, ], si)
    design <- NULL
    denom <- denoms
    runEdgeR(cds, FALSE, FALSE, TRUE, design, denom = denom, RUVg = TRUE)
  }, mc.cores = n_cpu)

message(paste('running edgeR with spike-ins', Sys.time()))
all_results_kal$edgeR_denom <- mclapply(1:N_SIM,
  function(i) {
    cat('Sample: ', i, '\n')
    si <- sample_info[[i]]
    counts <- obs_counts[[i]]
    counts <- round(counts)
    mode(counts) <- 'integer'
    filter <- sleuth_filters[[i]]
    cds <- make_count_data_set(counts[filter, ], si)
    design <- NULL
    denom <- denoms
    runEdgeR(cds, FALSE, FALSE, TRUE, design, denom = denom)
  }, mc.cores = n_cpu)

message(paste('getting benchmarks', Sys.time()))
saveRDS(all_results_kal, file = file.path('../results', sim_name, 'all_kal_results.rds'))
all_results_kal <- absSimSeq::rename_fc_list(all_results_kal)
saveRDS(all_results_kal, file = file.path('../results', sim_name, 'all_kal_fc_results.rds'))
oracles <- readRDS('../sims/polyester_oracles.rds')
all_kal_benchmarks <- lapply(seq_along(all_results_kal),
  function(i) {
    method_name <- names(all_results_kal)[i]
    result <- all_results_kal[[i]]
    mclapply(seq_along(result),
      function(j) {
        x <- result[j]
        if(grepl('sleuthALR', method_name) | grepl('ALDEx2', method_name)) {
          oracle <- oracles$alr_oracle[[j]]
        } else {
          oracle <- oracles$main_oracle[[j]]
        }
        filter <- sleuth_filters[[j]]
        filt_oracle <- oracle[which(oracle$target_id %in% names(filter[filter])), ]
        fold_change <- FALSE
        new_de_benchmark(x, method_name, filt_oracle, de_colors = full_method_colors[method_name],
                         de_linetypes = method_kal_ltys[method_name], fold_change = fold_change)
      }, mc.cores = n_cpu)
  })
names(all_kal_benchmarks) <- names(all_results_kal)

small_kal_benchmarks <- lapply(all_kal_benchmarks, function(x) x[1:GROUP_BREAKS[1]])
down_kal_benchmarks <- lapply(all_kal_benchmarks, function(x) x[(GROUP_BREAKS[1]+1):(GROUP_BREAKS[2])])
up_kal_benchmarks <- lapply(all_kal_benchmarks, function(x) x[(GROUP_BREAKS[2]+1):N_SIM])

dir.create(file.path('../results', sim_name), showWarnings = F)
saveRDS(small_kal_benchmarks, file = file.path('../results', sim_name,
  'small_kal_benchmarks.rds'))
saveRDS(down_kal_benchmarks, file = file.path('../results', sim_name,
  'down_kal_benchmarks.rds'))
saveRDS(up_kal_benchmarks, file = file.path('../results', sim_name,
  'up_kal_benchmarks.rds'))
