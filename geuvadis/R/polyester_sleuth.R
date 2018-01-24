.libPaths(c("~/R_library", .libPaths()))
n_cpu <- 20

geuvadis_dir <- "../results/finn_samples"
source("benchmark_methods.R")
source("gene_common.R")

suppressMessages(library("parallel"))

transcript_gene_mapping <- get_human_gene_names()

load('../metadata/geu_meta.RData', verbose = TRUE)
`%>%` <- dplyr::`%>%`

sample_info <- geu_meta %>%
  dplyr::filter(population == "FIN") %>%
  dplyr::filter(sex == "female")

all_geuvadis <- dir(geuvadis_dir)

proper_names <- lapply(finn_subset$sample,
  function(x) {
    grep(x, all_geuvadis, value = TRUE)
  })

sample_info <- dplyr::mutate(sample_info,
                             path = file.path(geuvadis_base,
                                              proper_names[sample],
                                              'abundance.h5'))

res <- run_sleuth(sample_info, formula = ~1)

saveRDS(res$so, file = '../results/finn_sleuth.rds')
