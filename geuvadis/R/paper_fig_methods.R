.libPaths(c("~/R_library", .libPaths()))

suppressPackageStartupMessages({
  library('cowplot')
  library('data.table')
  library('dplyr')
  library('mamabear')
  library('parallel')
  source('gene_common.R')
})

# theme used in Figures 2, S2-S4, and S8
theme_wam <- theme(
  axis.title.x = ggplot2::element_text(size = 10, face = "bold", family = "ArialMT"),
  axis.title.y = ggplot2::element_text(size = 10, face = "bold", family = "ArialMT"),
  axis.text.x = ggplot2::element_text(size = 8, family = "ArialMT"),
  axis.text.y = ggplot2::element_text(size = 8, family = "ArialMT"),
  legend.text = ggplot2::element_text(size = 8, family = "ArialMT"),
  legend.title = ggplot2::element_text(size = 10, family = "ArialMT")
)

# theme used in Figures 3-5, and S5-S7
theme_hp <- function() {
  theme_cowplot(12) +
    theme(legend.key.size = unit(2, "lines"), legend.position = "none",
          axis.title.x = ggplot2::element_text(size = 10, face = "bold", family = "ArialMT"),
          axis.title.y = ggplot2::element_text(size = 10, face = "bold", family = "ArialMT"),
          axis.text.x = ggplot2::element_text(size = 8, family = "ArialMT"),
          axis.text.y = ggplot2::element_text(size = 8, family = "ArialMT"))
}

### CODE FOR Figures 2 and S2 -- Main Benchmarks ###

make_benchmark_plots <- function(benchmarks, to_include, new_labels_kal) {
  benchmarks <- benchmarks[sort(names(benchmarks))]

  subset_benchmarks <- benchmarks[to_include]
  subset_names <- names(subset_benchmarks)
  new_names_kal <- new_labels_kal[match(subset_names, to_include)]
  subset_benchmarks <- lapply(seq_along(subset_benchmarks), function(i) {
    bench_list <- subset_benchmarks[[i]]
    original_label <- subset_names[i]
    new_label <- new_names_kal[i]
    lapply(bench_list, function(bench) {
      new_bench <- rename_benchmark(bench, original_label, new_label, join_mode = "intersect")
      if (grepl("RUVg", original_label)) {
        new_bench$line_mapping <- 4
        names(new_bench$line_mapping) <- new_label
      } else if (!grepl("denom", original_label)) {
        new_bench$line_mapping <- 2
        names(new_bench$line_mapping) <- new_label
      }
      if (original_label == "ALDEx2.denom.overlap") {
        new_bench$color_mapping <- c("ALDEx2+denom" = '#009E73')
      }
      new_bench
    })
  })
  color_mapping <- sapply(subset_benchmarks, function(x) x[[1]]$color_mapping)
  line_mapping <- sapply(subset_benchmarks, function(x) x[[1]]$line_mapping)
  names(subset_benchmarks) <- new_names_kal

  subset_benchmarks <- subset_benchmarks[sort(names(subset_benchmarks))]

  subset_fdr_data <- lapply(subset_benchmarks,
    function(bench) {
      fdr <- suppressMessages(get_fdr(bench))
      fdr$pvals
    })
  subset_fdr_data <- dplyr::bind_rows(subset_fdr_data)

  p <- ggplot(subset_fdr_data, aes(true_fdr, sensitivity))
  p <- p + geom_path(aes(color = method, linetype = method),
      size = 0.25, alpha = 0.8)
  p <- p + xlab("False Discovery Rate")
  p <- p + ylab("Sensitivity")
  p <- p + scale_color_manual(values = color_mapping)
  p <- p + scale_linetype_manual(values = line_mapping)
  p <- p + theme_wam

  zoom_p <- p + xlim(0, 0.25)
  zoom_p <- zoom_p + theme(axis.text.x = ggplot2::element_text(size = 6, family = "ArialMT"))

  list(full = p, zoom = zoom_p)
}

#### CODE FOR FIGURE 3 -- Bottomly Self-consistency Experiment ###

make_fig3 <- function(self_benchmark, mapping) {
  self_fdr <- mclapply(self_benchmark, average_sensitivity_specificity, mc.cores = 15)
  self_fdr <- dplyr::bind_rows(self_fdr)
  self_fdr <- dplyr::mutate(self_fdr, method = sub('qval_', '', method))
  self_fdr <- dplyr::filter(self_fdr, !grepl("(wilcoxon|welch)", method))
  self_fdr <- dplyr::mutate(self_fdr, method = mapping[method])
  
  if (!suppl) {
    self_fdr <- dplyr::filter(self_fdr, fdr_level == 0.1)
  } else {
    self_fdr <- dplyr::filter(self_fdr, fdr_level != 0.1)
  }

  self_fdr <- dplyr::mutate(self_fdr,
    fdr_level_string = paste0('eFDR = ', sprintf('%.2f', fdr_level)))
  fdr_count <- dplyr::select(self_fdr, method, fdr_level_string, true_fdr, sample) %>%
    dplyr::group_by(method, fdr_level_string) %>%
    dplyr::summarize(n = sum(!is.na(true_fdr)))
  fdr_count <- dplyr::mutate(fdr_count, y = 1)

  p <- ggplot(self_fdr, aes(method, true_fdr, color = method)) +
    geom_boxplot(outlier.shape = NA) + #avoid plotting outliers twice
    geom_jitter(size = 0.5, position = position_jitter(width = 0.375, height = 0), alpha = 0.5) +
    geom_text(aes(method, y, label = n), data = fdr_count) +
    facet_wrap(~fdr_level_string, ncol = 1) +
    geom_hline(aes(yintercept = fdr_level), linetype = 2) +
    ylim(0, 1) +
    ylab('false discovery rate') +
    theme_hp() +
    scale_color_manual(values = method_colors) # method_colors in gene_common.R
  p
}

### CODE FOR FIGURE 4 -- GEUVADIS NULL RESAMPLING EXPERIMENT ###

make_fig4 <- function(self_benchmark) {
  self_fdr <- lapply(self_benchmark, average_sensitivity_specificity,
    use_oracle = TRUE, use_fdr = TRUE)
  self_fdr <- dplyr::bind_rows(self_fdr)
  self_fdr <- dplyr::mutate(self_fdr, method = sub('pval_', '', method))
  self_fdr <- dplyr::mutate(self_fdr,
    fdr_level_string = paste0('eFDR = ', sprintf('%.2f', fdr_level)))

  library(scales)
  S_log10 <- function(x) {ifelse(x == 0, 0, sign(x) * log10(abs(x)))}
  IS_log10 <- function(x) {ifelse(x == 0, 0, 10^abs(x) * sign(x))}
  S_log10_trans <- function() trans_new("S_log10", S_log10, IS_log10)

  p <- ggplot(self_fdr, aes(method, fp, color = method)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~fdr_level_string, ncol = 1) +
    scale_y_continuous(trans = "S_log10", breaks = c(0, 10, 30, 100, 300, 1000, 3000, 10000, 30000)) +
    ylab('number of false positives') +
    theme_hp() + theme(axis.text.x = ggplot2::element_text(family = "ArialMT", angle = 45, vjust = 1, hjust=1),
                       panel.grid.major.y = element_line(colour = "grey92")) +
    geom_jitter(size = 0.5, height = 0) +
    scale_color_manual(values = method_colors) # method_colors in gene_common.R
  p
}

## CODE FOR Figure S3 -- ALDEx2 INTERNAL COMPARISON ###

make_aldex_plot <- function(benchmarks, to_keep, new_labels_kal) {
  benchmarks <- benchmarks[sort(names(benchmarks))]
  subset_benchmarks <- benchmarks[names(benchmarks) %in% to_keep]
  subset_names <- names(subset_benchmarks)
  new_names_kal <- new_labels_kal[match(subset_names, to_keep)]
  subset_benchmarks <- lapply(subset_benchmarks, function(bench) {
    bench_list <- subset_benchmarks[[i]]
    original_label <- subset_names[i]
    new_label <- new_names_kal[i]
    lapply(bench_list, function(bench) {
      new_bench <- rename_benchmark(bench, original_label, new_label, join_mode = "intersect")
      if (grepl("iqlr", original_label)) {
        new_bench$line_mapping <- 2
      } else if (grepl("clr", original_label)) {
        new_bench$line_mapping <- 4
      }
      names(new_bench$line_mapping) <- new_label
      new_bench
    })
  })
  color_mapping <- sapply(subset_benchmarks, function(x) x[[1]]$color_mapping)
  line_mapping <- sapply(subset_benchmarks, function(x) x[[1]]$line_mapping)
  names(subset_benchmarks) <- new_names_kal
  subset_benchmarks <- subset_benchmarks[sort(names(subset_benchmarks))]

  subset_fdr_data <- lapply(subset_benchmarks,
    function(bench) {
      fdr <- suppressMessages(get_fdr(bench))
      fdr$pvals
    })
  subset_fdr_data <- dplyr::bind_rows(subset_fdr_data)
  p <- ggplot(subset_fdr_data, aes(true_fdr, sensitivity))
  p <- p + geom_path(aes(color = method, linetype = method),
      size = 0.25, alpha = 0.8)
  p <- p + xlab("False Discovery Rate")
  p <- p + ylab("Sensitivity")
  p <- p + scale_color_manual(values = color_mapping)
  p <- p + scale_linetype_manual(values = line_mapping)
  p <- p + theme_wam
  p
}

## CODE FOR Figure S4 -- Internal sleuth and sleuth-ALR Comparison ###

make_sleuth_plot <- function(benchmarks, to_keep, new_labels_kal) {
  benchmarks <- benchmarks[sort(names(benchmarks))]
  subset_benchmarks <- benchmarks[names(benchmarks) %in% to_keep]
  subset_names <- names(subset_benchmarks)
  new_names_kal <- new_labels_kal[match(subset_names, to_keep)]
  subset_benchmarks <- lapply(subset_benchmarks, function(bench) {
    bench_list <- subset_benchmarks[[i]]
    original_label <- subset_names[i]
    new_label <- new_names_kal[i]
    lapply(bench_list, function(bench) {
      new_bench <- rename_benchmark(bench, original_label, new_label, join_mode = "intersect")
      if (grepl("wt", original_label)) {
        new_bench$line_mapping <- 2
      }
      names(new_bench$line_mapping) <- new_label
      new_bench
    })
  })
  color_mapping <- sapply(subset_benchmarks, function(x) x[[1]]$color_mapping)
  line_mapping <- sapply(subset_benchmarks, function(x) x[[1]]$line_mapping)
  names(subset_benchmarks) <- new_names_kal
  subset_benchmarks <- subset_benchmarks[sort(names(subset_benchmarks))]
  
  subset_fdr_data <- lapply(subset_benchmarks,
    function(bench) {
      fdr <- suppressMessages(get_fdr(bench))
      fdr$pvals
    })
  subset_fdr_data <- dplyr::bind_rows(subset_fdr_data)
  p <- ggplot(subset_fdr_data, aes(true_fdr, sensitivity))
  p <- p + geom_path(aes(color = method, linetype = method),
      size = 0.25, alpha = 0.8)
  p <- p + xlab("False Discovery Rate")
  p <- p + ylab("Sensitivity")
  p <- p + scale_color_manual(values = color_mapping)
  p <- p + scale_linetype_manual(values = line_mapping)
  p <- p + theme_wam
  p
}

## CODE for Figure S8 -- Internal sleuth-ALR imputation comparison ###

make_impute_plots <- function(benchmarks, original_labels, new_labels_kal) {
  benchmarks <- lapply(benchmarks, function(bench) {
    rename_benchmark(bench, original_labels, new_labels_kal, join_mode = "intersect") 
  })
  fdr <- suppressMessages(get_fdr(benchmarks)$pvals)
  fdr <- fdr[!grepl("lrt", fdr$method), ]
  counts_fdr <- fdr[grepl("(counts|^sleuth$)", fdr$method), ]
  tpm_fdr <- fdr[grepl("(TPM|^sleuth$)", fdr$method), ]
  p <- ggplot(counts_fdr, aes(true_fdr, sensitivity))
  p <- p + geom_path(aes(color = method, linetype = method),
      size = 0.25, alpha = 0.8)
  p <- p + xlab("False Discovery Rate")
  p <- p + ylab("Sensitivity")
  p <- p + theme_wam
  c <- p
  
  p <- ggplot(tpm_fdr, aes(true_fdr, sensitivity))
  p <- p + geom_path(aes(color = method, linetype = method),
      size = 0.25, alpha = 0.8)
  p <- p + xlab("False Discovery Rate")
  p <- p + ylab("Sensitivity")
  p <- p + theme_wam

  list(counts_graph = c, tpms_graph = p) 
}
