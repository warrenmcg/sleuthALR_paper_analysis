source('paper_fig_methods.R')

base_dir <- '../results/final_figures'
default_extension <- '.pdf'

message('loading benchmarks')
sim_name <- "isoform_5_5_15_387_1"
sim_types <- c('small', 'down', 'up')
suffix <- 'benchmarks.rds'
benchmarks <- lapply(sim_types, function(sim_type) {
  message(sim_type)
  kal_benchmarks <- readRDS(paste0('../results/', sim_name,
    '/', sim_type, paste0('_kal_', suffix)))
  kal_benchmarks
})

## FIGURE 2 and S2 ##
message('making figures 2 and S2')
to_include <- c('sleuth.lrt', 'ALDEx2.iqlr.overlap', 'DESeq2', 'edgeR',
  'limmaVoom', 'DESeq2_RUVg', 'sleuthALR.lrt', 'ALDEx2.denom.overlap',
  'DESeq2_denom', 'edgeR_denom', 'limmaVoom_denom', 'edgeR_RUVg')

new_labels <- c('sleuth', 'ALDEx2 IQLR', 'DESeq2', 'edgeR',
  'limma', 'DESeq2+RUVg', 'sleuth-ALR', 'ALDEx2+denom',
  'DESeq2+denom', 'edgeR+denom', 'limma+denom', 'edgeR+RUVg')

graphs <- lapply(1:3, function(i) {
  b <- benchmarks[[i]]
  invisible(make_benchmark_plots(b, to_include, new_labels))
})

full_graphs <- lapply(graphs, '[[', 'full')
zoom_graphs <- lapply(graphs, '[[', 'zoom')
g <- full_graphs[[1]] +
  guides(colour = guide_legend(ncol=2, title.hjust = 0.5),
         linetype = guide_legend(ncol=2, title.hjust = 0.5)) +
  theme(legend.justification="center",
        legend.box.just = "bottom")
legend <- get_legend(g)
p <- plot_grid(full_graphs[[1]] + theme(legend.position = "none"),
               full_graphs[[2]] + theme(legend.position = "none"),
               full_graphs[[3]] + theme(legend.position = "none"),
               legend, labels = c("A","B","C",""), align = 'center',
               nrow = 2, label_size = 12, hjust = -1.2)
filename <- file.path(base_dir, paste0('figureS2', default_extension))
save_plot(filename, p, base_width = 6.5)

g <- zoom_graphs[[1]] +
  guides(colour = guide_legend(ncol=2, title.hjust = 0.5),
         linetype = guide_legend(ncol=2, title.hjust = 0.5)) +
  theme(legend.justification="center",
        legend.box.just = "bottom")
legend <- get_legend(g)
p <- plot_grid(zoom_graphs[[1]] + theme(legend.position = "none"),
               zoom_graphs[[2]] + theme(legend.position = "none"),
               zoom_graphs[[3]] + theme(legend.position = "none"),
               legend, labels = c("A","B","C",""), align = 'center',
               nrow = 2, label_size = 12, hjust = -1.2)
filename <- file.path(base_dir, paste0('figure2', default_extension))
save_plot(filename, p, base_width = 6.5)

#### CODE FOR FIGURE 3 -- Bottomly Self-consistency Experiment ###
message("making Figure 3 and Figure S5")
mapping <- c('sleuth.lrt' = 'sleuth LRT', 'sleuth.wt' = 'sleuth Wald',
  'sleuthALR.lrt' = 'sleuth-ALR\nLRT', 'sleuthALR.wt' = 'sleuth-ALR\nWald',
  'limmaVoom' = 'limma +\nC.N.', 'limmaVoom_old' = 'limma',
  'DESeq2' = 'DESeq2 +\nC.N.', 'DESeq2_old' = 'DESeq2', 'DESeq2_RUVg' = 'DESeq2 +\nRUVg',
  'edgeR' = 'edgeR +\nC.N.', 'edgeR_old' = 'edgeR', 'edgeR_RUVg' = 'edgeR +\nRUVg',
  'ALDEx2.filt.overlap' = 'ALDEx2\noverlap')

iso_file <- '../../bottomly/results/isoform_self_benchmark.rds'
bench <- readRDS(iso_file)
iso_p <- make_fig3(bench, mapping)

p <- plot_grid(iso_p, labels = "AUTO", align = "vh", nrow = 1, label_size = 24, hjust = -0.5, vjust = 1.25)
filename <- file.path(base_dir, paste0('figure3', default_extension))
save_plot(filename, p, base_aspect_ratio = 1.5, base_width = 6)

iso_p <- make_fig3(bench, suppl = TRUE)

p <- plot_grid(iso_p, labels = "", align = "vh", ncol = 1, label_size = 20, hjust = -0.5, vjust = 1.25, axis = "bottom")
filename <- file.path(base_dir, paste0('figureS5', default_extension))
save_plot(filename, p, base_aspect_ratio = 1.5, base_width = 6)

### FIGURE 4 -- GEUVADIS NULL RESAMPLING EXPERIMENT ###
message("making figure 4")
self_benchmark <- readRDS('../results/null_resampling/isoform.rds')
iso_p <- make_fig4(self_benchmark)
p <- plot_grid(iso_p, labels = "AUTO", align = "vh", nrow = 1, label_size = 24, hjust = -0.5, vjust = 1.25)
filename <- file.path(base_dir, paste0('figure4', default_extension))
save_plot(filename, p, base_aspect_ratio = 1.5, base_height = 6)

## aldex parameters
message('making figure S3')
to_keep <- c(
  'ALDEx2.iqlr.overlap', 'ALDEx2.iqlr.welch', 'ALDEx2.iqlr.wilcoxon',
  'ALDEx2.clr.overlap', 'ALDEx2.clr.welch', 'ALDEx2.clr.wilcoxon',
  'ALDEx2.denom.overlap', 'ALDEx2.denom.welch', 'ALDEx2.denom.wilcoxon'
)
new_labels <- c(
  'ALDEx2 IQLR\noverlap', 'ALDEx2 IQLR\nWelch', 'ALDEx2 IQLR\nWilcoxon',
  'ALDEx2 CLR\noverlap', 'ALDEx2 CLR\nWelch', 'ALDEx2 CLR\nWilcoxon',
  'ALDEx2 denom\noverlap', 'ALDEx2 denom\nWelch', 'ALDEx2 denom\nWilcoxon'
)

aldex_graphs <- lapply(1:3, function(i) {
  b <- benchmarks[[i]]
  make_aldex_plot(b, to_keep, new_labels)
})

g <- aldex_graphs[[1]] +
  guides(colour = guide_legend(ncol=3, byrow = TRUE, title.hjust = 0.5),
         linetype = guide_legend(ncol=3, byrow = TRUE, title.hjust = 0.5)) +
  theme(legend.justification="center",
        legend.box.just = "bottom")#,
legend <- get_legend(g)
p <- plot_grid(aldex_graphs[[1]] + theme(legend.position = "none"),
               aldex_graphs[[2]] + theme(legend.position = "none"),
               aldex_graphs[[3]] + theme(legend.position = "none"),
               legend, labels = c("A","B","C",""), align = 'center',
               nrow = 2, label_size = 12, hjust = -1.2)
filename <- file.path(base_dir, paste0('figureS3', default_extension))
save_plot(filename, p, base_width = 6.5)

## sleuth parameters
message('making figure S4')
to_keep <- c('sleuth.lrt', 'sleuthALR.lrt',
  'sleuth.wt', 'sleuthALR.wt',
  'sleuthALR.counts.lrt', 'sleuthALR.counts.wt'
)
new_labels <- c('sleuth LRT', 'sleuth-ALR\nTPM LRT',
  'sleuth Wald', 'sleuth-ALR\nTPM Wald',
  'sleuth-ALR\nCounts LRT', 'sleuth-ALR\nCounts Wald')

sleuth_graphs <- lapply(1:3, function(i) {
  b <- benchmarks[[i]]
  make_sleuth_plot(b, to_keep, new_labels)
})

g <- sleuth_graphs[[1]] +
  guides(colour = guide_legend(ncol=2, title.hjust = 0.5),
         linetype = guide_legend(ncol=2, title.hjust = 0.5)) +
  theme(legend.justification="center",
        legend.box.just = "bottom")
legend <- get_legend(g)
p <- plot_grid(sleuth_graphs[[1]] + theme(legend.position = "none"),
               sleuth_graphs[[2]] + theme(legend.position = "none"),
               sleuth_graphs[[3]] + theme(legend.position = "none"),
               legend, labels = c("A","B","C",""), align = 'center',
               nrow = 2, label_size = 12, hjust = -1.2)
filename <- file.path(base_dir, paste0('figureS4', default_extension))
save_plot(filename, p, base_width = 6.5)

## impute parameters
message('making figure S6')
original_labels <- c('sleuth.wt', 'sleuthALR.wt',
  'sleuthALR.0.1.wt', 'sleuthALR.0.01.wt', 'sleuthALR.0.001.wt',
  'sleuthALR.1E4.wt', 'sleuthALR.counts.wt',
  'sleuthALR.counts.0.001.wt', 'sleuthALR.counts.0.01.wt',
  'sleuthALR.counts.0.1.wt', 'sleuthALR.counts.0.5.wt'
)
new_labels <- c('sleuth', 'sleuth-ALR TPM',
  'sleuth-ALR TPM 1e-1', 'sleuth-ALR TPM 1e-2', 'sleuth-ALR TPM 1e-3',
  'sleuth-ALR TPM 1e-4', 'sleuth-ALR counts',
  'sleuth-ALR counts 1e-3', 'sleuth-ALR counts 0.01',
  'sleuth-ALR counts 0.1', 'sleuth-ALR counts 0.5'
)

impute_graphs <- lapply(1:3, function(i) {
  b <- benchmarks[[i]]
  make_impute_plots(b, original_labels, new_labels)
  })

count_graphs <- lapply(impute_graphs, '[[', 'counts_graph')
tpm_graphs <- lapply(impute_graphs, '[[', 'tpms_graph')

prow <- plot_grid(count_graphs[[1]] + theme(legend.position="none"),
  count_graphs[[2]] + theme(legend.position="none"),
  count_graphs[[3]] + theme(legend.position="none"),
  labels = "AUTO", align = 'vh', nrow = 1, label_size = 12, hjust = -1.2)
legend <- get_legend(count_graphs[[1]])
p1 <- plot_grid(prow, legend, rel_widths = c(3, 1))

tpm_prow <- plot_grid(tpm_graphs[[1]] + theme(legend.position="none"),
  tpm_graphs[[2]] + theme(legend.position="none"),
  tpm_graphs[[3]] + theme(legend.position="none"),
  labels = c("D", "E", "F"), align = 'vh', nrow = 1, label_size = 12, hjust = -1.2)
legend <- get_legend(tpm_graphs[[1]])
p2 <- plot_grid(tpm_prow, legend, rel_widths = c(3, 1))

p <- plot_grid(p1, p2, align = 'v', ncol = 1)
filename <- file.path(base_dir, paste0('figureS6', default_extension))
save_plot(filename, p, base_aspect_ratio = 1.75, base_height = 3.76)
