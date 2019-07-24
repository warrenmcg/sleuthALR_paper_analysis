.libPaths("~/R_library", .libPaths()))
message('loading packages')
suppressMessages({
  library(absSimSeq)
  library(ggplot2)
  library(cowplot)
})

data(ERCC92_data, package = "absSimSeq")
all_spikeins <- rownames(ERCC92_data)
logmeans <- log2(rowMeans(ERCC92_data[,c(4:5)]))
denoms <- names(logmeans[logmeans > 3])

sim_name <- 'isoform_5_5_30_645175_1'
base_dir <- '../results/final_figures'
default_extension <- '.pdf'

message('processing ground truth counts')
true_counts <- lapply(1:15, function(i) {
  file_name <- file.path(paste0('run', sprintf('%02d', i)),
                         'sim_counts_matrix.rda')
  load(file.path('../sims', sim_name, file_name))
  counts_matrix
})

run_names <- c('Small', 'Down', 'Up')

spike_percent <- sapply(true_counts, function(x)  {
  colSums(x[all_spikeins, ]) / colSums(x)
})
rownames(spike_percent) <- paste0('sample_', sprintf('%02d',1:10))
colnames(spike_percent) <- paste0('run', sprintf('%02d', 1:15))
spike_df <- as.data.frame(spike_percent)
spike_df$sample <- rownames(spike_df)
spike_df$condition <- c(rep("Ctr", 5), rep("Exp", 5))
spike_count_df <- tidyr::gather(spike_df, run, percent, -sample, -condition)

message('loading ground truth expected reads per transcript')
load('../results/polyester_ground_truth.RData')
x <- final_results$results[[1]]$expected_reads

theme_hp <- function() {
  theme_cowplot(12) +
    theme(legend.key.size = unit(2, "lines"),
          axis.title.x = ggplot2::element_text(size = 10, face = "bold", family = "ArialMT"),
          axis.title.y = ggplot2::element_text(size = 10, face = "bold", family = "ArialMT"),
          axis.text.x = ggplot2::element_text(family = "ArialMT"),
          axis.text.y = ggplot2::element_text(family = "ArialMT"))
}

message('making figure 3')
graphs <- lapply(1:3, function(i) {
  runs_to_use <- paste0("run", sprintf('%02d', ((i-1)*5+1):(5*i)))
  data <- spike_count_df[which(spike_count_df$run %in% runs_to_use),]
  expected_percent <- c(colSums(x[all_spikeins,]) / colSums(x))[1]
  x_label <- paste(run_names[i], "Experiments")
  ggplot(data, aes(run, percent, group = sample, fill = condition)) +
    geom_bar(position = "dodge", stat = "identity", color = 'black') +
    geom_hline(yintercept = expected_percent, linetype = "dotted") +
    scale_fill_manual(values = c('Ctr' = 'white', 'Exp' = '#E69F00')) +
    xlab(x_label) +
    ylab("Reads Mapping\nto Spike-ins (%)") +
    theme_hp()
})

legend <- get_legend(graphs[[1]] + theme(legend.justification="center",
                                         legend.box.just = "bottom"))
fig3 <- plot_grid(graphs[[1]] + theme(legend.position = "none"),
                   graphs[[2]] + theme(legend.position = "none"),
                   graphs[[3]] + theme(legend.position = "none"),
                   legend, labels = c("A","B","C",""), align = 'center',
                   nrow = 2, label_size = 12, hjust = -1.2)
filename <- file.path(base_dir, paste0('figure3', default_extension))
save_plot(filename, fig3, base_width = 6.5)

message('making figure S5')
spike_fc <- sapply(true_counts, function(x) {
  spikeins <- x[all_spikeins, ]
  log2(rowMeans(spikeins[, 6:10]) / rowMeans(spikeins[, 1:5]))
})

colnames(spike_fc) <- paste0('run', sprintf('%02d', 1:15))
spike_fc <- as.data.frame(spike_fc)
spike_fc$spikein <- rownames(spike_fc)
spike_fc$log2Mean <- log2(rowMeans(ERCC92_data[,c(4:5)]))
spike_fc$highBool <- spike_fc$log2Mean > 3
spike_fc_df <- tidyr::gather(spike_fc, run, fc, -spikein, -log2Mean, -highBool)

spike_geo_means <- sapply(true_counts, function(x) {
  spikeins <- x[denoms, ]
  geomeans <- apply(spikeins, 2, sleuthALR::geomean)
  log2(mean(geomeans[6:10]) / mean(geomeans[1:5]))
})

spike_deseq_means <- sapply(true_counts, function(x) {
  spikeins <- x[denoms, ]
  deseq_sf <- DESeq2::estimateSizeFactorsForMatrix(spikeins)
  log2(mean(deseq_sf[6:10]) / mean(deseq_sf[1:5]))
})

cns <- sapply(final_results$results, function(x) {
  colSums(x$copy_numbers_per_cell)
})

ideal_fc <- log2(cns[1, ] / cns[2, ])
x_vals <- ifelse(1:15 %% 5 == 0, 5, 1:15 %% 5)
cns_df <- data.frame(`RNA Content Change` = ideal_fc,
                     `Geometric Means` = spike_geo_means,
                     `DESeq2 Size Factors` = spike_deseq_means,
                     x = x_vals - 0.375,
                     xend  = x_vals + 0.375)
cns_df$run <- colnames(spike_fc)[1:15]
cns_df <- tidyr::gather(cns_df, method, y, -x, -xend, -run)
cns_df$method <- gsub('.', ' ', cns_df$method, fixed = TRUE)
colors <- c(`RNA Content Change` = "#FF0000",
            `Geometric Means` = "#E69F00",
            `DESeq2 Size Factors` = "#2166AC")
cns_df$color <- colors[cns_df$method]

spike_graphs <- lapply(1:3, function(i) {
  runs_to_use <- paste0("run", sprintf('%02d', ((i-1)*5+1):(5*i)))
  data <- spike_fc_df[which(spike_fc_df$run %in% runs_to_use),]
  x_label <- paste(run_names[i], "Experiments")
  cns_to_use <- cns_df[which(cns_df$run %in% runs_to_use), ]
  p <- ggplot(data, aes(run, fc, color = highBool)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_segment(data = cns_to_use, aes(x = x, xend = xend,
                                        y = y, yend = y),
                 linetype = "dotted", color = cns_to_use$color) +
    scale_color_manual(values = c('black', '#990000'),
                       name = "Spike-in Expression",
                       labels = c("low", "high")) +
    xlab(x_label) +
    ylab("Spike-in Fold Changes\n(Log2)") +
    theme_hp() +
    theme(axis.text.x = element_text(family = "ArialMT",
                                     angle = 90, vjust = 0.5))
  p
})
legend <- get_legend(spike_graphs[[1]] + theme(legend.justification="center",
                                               legend.box.just = "bottom"))
figS5 <- plot_grid(spike_graphs[[1]] + theme(legend.position = "none"),
                  spike_graphs[[2]] + theme(legend.position = "none"),
                  spike_graphs[[3]] + theme(legend.position = "none"),
                  legend, labels = c("A","B","C",""), align = 'center',
                  nrow = 2, label_size = 12, hjust = -1.2)
filename <- file.path(base_dir, paste0('figureS5', default_extension))
save_plot(filename, figS5, base_height = 6)

