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
options(mc.cores = n_cpu)

transcript_gene_mapping <- get_human_gene_names()

###
# new refactoring
###

all_results_sal <- list()

N_SIM <- 15
GROUP_BREAKS <- c(5, 10)

message(paste('running sleuth', Sys.time()))
sleuth_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    sir <- load_isoform_results_intersect(sim_name, i, 'salmon', 'sleuth',
      run_sleuth)
    list(sleuth.lrt = sir$sleuth.lrt, sleuth.wt = sir$sleuth.wt)
  })
all_results_sal$sleuth.lrt <- lapply(sleuth_res, '[[', 'sleuth.lrt')
all_results_sal$sleuth.wt <- lapply(sleuth_res, '[[', 'sleuth.wt')
rm(sleuth_res)

message(paste('running sleuth-ALR', Sys.time()))
alr_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    alr <- load_isoform_results_intersect(sim_name, i, 'salmon', 'sleuthALR',
      run_alr, denom = 'ENST00000466430.5')
    list(sleuthALR.lrt = alr$sleuthALR.lrt, sleuthALR.wt = alr$sleuthALR.wt)
  })
all_results_sal$sleuthALR.lrt <- lapply(sleuth_res, '[[', 'sleuthALR.lrt')
all_results_sal$sleuthALR.wt <- lapply(sleuth_res, '[[', 'sleuthALR.wt')
rm(alr_res)

message(paste('running sleuth-ALR', Sys.time()))
aldex2_res <- mclapply(1:N_SIM,
  function(i) {
    cat('Run: ', i, '\n')
    aldex2 <- load_isoform_results_intersect(sim_name, i, 'salmon', 'ALDEx2',
      aldex2_filter_and_run, denom = 'ENST00000466430.5')
    list(ALDEx2.overlap = aldex2$ALDEx2.overlap,
         ALDEx2.welsh = aldex2$ALDEx2.welsh,
         ALDEx2.wilcoxon = aldex2$ALDEx2.wilcoxon)
  })
all_results_sal$ALDEx2.overlap <- lapply(aldex2_res, '[[', 'ALDEx2.overlap')
all_results_sal$ALDEx2.welsh <- lapply(aldex2_res, '[[', 'ALDEx2.welsh')
all_results_sal$ALDEx2.wilcoxon <- lapply(aldex2_res, '[[', 'ALDEx2.wilcoxon')

message(paste('running limma', Sys.time()))
all_results_sal$limmaVoom <- mclapply(1:N_SIM,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'salmon',
      'limmaVoom', limma_filter_and_run)$limmaVoom
  }, mc.cores = n_cpu)

message(paste('running DESeq2', Sys.time()))
all_results_sal$DESeq2 <- mclapply(1:N_SIM,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'salmon',
      'DESeq2', DESeq2_filter_and_run_intersect)$DESeq2
  }, mc.cores = n_cpu)

message(paste('running edgeR', Sys.time()))
all_results_sal$edgeR <- mclapply(1:N_SIM,
  function(i) {
    cat('Sample: ', i, '\n')
    load_isoform_results_intersect(sim_name, i, 'salmon',
       'edgeR', edgeR_filter_and_run)$edgeR
  }, mc.cores = n_cpu)

message(paste('getting benchmarks', Sys.time()))
oracles <- readRDS('../sims/polyester_oracles.rds')
all_sal_benchmarks <- lapply(seq_along(all_results_sal),
  function(i) {
    method_name <- names(all_results_sal)[i]
    result <- all_results_sal[[i]]
    mclapply(seq_along(result),
      function(j) {
        x <- result[j]
        if(grepl('sleuthALR', method_name) | grepl('ALDEx2', method_name)) {
          oracle <- oracles$alr_oracle[[j]]
        } else {
          oracle <- oracles$main_oracle[[j]]
        }
        new_de_benchmark(x, method_name, oracle)
      }, mc.cores = n_cpu)
  })
names(all_sal_benchmarks) <- names(all_results_sal)

small_sal_benchmarks <- lapply(all_sal_benchmarks, function(x) x[1:GROUP_BREAKS[1]])
down_sal_benchmarks <- lapply(all_sal_benchmarks, function(x) x[(GROUP_BREAKS[1]+1):(GROUP_BREAKS[2])])
up_sal_benchmarks <- lapply(all_sal_benchmarks, function(x) x[(GROUP_BREAKS[2]+1):N_SIM])

dir.create(file.path('../results', sim_name), showWarnings = F)
saveRDS(small_sal_benchmarks, file = file.path('../results', sim_name,
  'small_sal_benchmarks.rds'))
saveRDS(down_sal_benchmarks, file = file.path('../results', sim_name,
  'down_sal_benchmarks.rds'))
saveRDS(up_sal_benchmarks, file = file.path('../results', sim_name,
  'up_sal_benchmarks.rds'))
