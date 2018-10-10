.libPaths(c("~/R_library", .libPaths()))

sim_name <- 'isoform_5_5_15_387_1'
target_id <- 'ENST00000561452.5'
n_cpu <- 10

suppressMessages({
  source("benchmark_methods.R")
  source("gene_common.R")

  library("parallel")
  library("mamabear")
  library("cowplot")
  library("absSimSeq")
  options(mc.cores = n_cpu)
})

transcript_gene_mapping <- get_human_gene_names()

###
# new refactoring
###

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
si <- sample_info[[6]]
so <- sleuth_prep(si, max_bootstrap = 2, num_cores = n_cpu)
filter <- so$filter_bool
filter_ids <- names(filter)[filter]

data('ERCC92_data', package = "absSimSeq")
log_means <- log2(rowMeans(ERCC92_data[,c(4:5)]))
denoms <- names(log_means[log_means > 3])
denom <- gsub('_', '-', denoms)

message(paste('running sleuth-ALR', Sys.time()))
alr <- run_alr(si, max_bootstrap = 100, denom = denom,
               num_cores = n_cpu, which_var = "obs_tpm", filter_target_id = filter_ids)
alr_so <- alr$so

message(paste('running sleuth-ALR with delta = 0.1', Sys.time()))
alr01 <- run_alr(si, max_bootstrap = 100, denom = denom, delta = 0.1,
                 num_cores = n_cpu, which_var = "obs_tpm", filter_target_id = filter_ids)
alr01_so <- alr01$so

theme_hp <- function() {
  theme_cowplot(12) +
    theme(legend.key.size = unit(2, "lines"), legend.position = "none",
          axis.title.x = ggplot2::element_text(face = "bold", family = "ArialMT"),
          axis.title.y = ggplot2::element_text(face = "bold", family = "ArialMT"),
          axis.text.x = ggplot2::element_text(family = "ArialMT", angle = 50, hjust = 1),
          axis.text.y = ggplot2::element_text(family = "ArialMT"),
          panel.grid.major.y = ggplot2::element_line(size = 0.25,color = "gray"))
}

alr_graph <- sleuth::plot_bootstrap(alr_so, target_id = target_id, units = "tpm", color_by = "condition")
alr01_graph <- sleuth::plot_bootstrap(alr01_so, target_id = target_id, units = "tpm", color_by = "condition")

alr_graph <- alr_graph + theme_hp() +
  scale_fill_manual(values = c('A' = 'white', 'B' = 'grey')) +
  ggtitle('Impute < Smallest Observed') +
  ylab('TPM logratio') +
  ylim(c(-41, 0))
alr01_graph <- alr01_graph + theme_hp() +
  scale_fill_manual(values = c('A' = 'white', 'B' = 'grey')) +
  ggtitle('Impute = 0.01') +
  ylab('TPM logratio') +
  ylim(c(-41, 0))

p <- plot_grid(alr_graph, alr01_graph, labels = "AUTO", align = "vh", nrow = 1, label_size = 16, hjust = -1)

save_plot('../results/final_figures/figureS7.pdf', p, base_aspect_ratio = 2, base_height = 3)
