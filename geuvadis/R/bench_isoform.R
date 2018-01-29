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
### TO-DO: add sleuth
### TO-DO: add sleuth-ALR
### TO-DO: add ALDEx2

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
all_sal_benchmarks <- lapply(all_results_sal,
  function(result) {
    mclapply(seq_along(result),
      function(i) {
        x <- result[[i]]
        sim_info <- get_de_info(sim_name, i, transcript_gene_mapping)
        new_de_benchmark(x, names(x), sim_info$de_info)
      }, mc.cores = n_cpu)
  })

### TO-DO: Convert this to the separate benchmarks
saveRDS(all_benchmarks, file = paste0('../results/', sim_name,
  '/isoform_benchmarks.rds'))
