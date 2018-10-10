library(ggplot2)
library(cowplot)

yeast_annos <- read.table('~/sleuth_paper_analysis/annotation/ASM294v2.pombase.all_annos.txt', header = T, sep = "\t", quote = "", stringsAsFactors = F)
res <- readRDS('../results/yeast_results.rds')
sleuth_table <- res$sleuth.wt
sleuth_table$gene_id <- yeast_annos$gene_id[match(sleuth_table$target_id, yeast_annos$target_id)]
alr_table <- res$sleuthALR.wt
alr_table$gene_id <- yeast_annos$gene_id[match(alr_table$target_id, yeast_annos$target_id)]

abs_counts <- read.table('../abs_counts/abs_counts.txt', header = T, sep = "\t", stringsAsFactors=F)
abs_counts$meanMM <- apply(abs_counts[, 3:4], 1, mean)
abs_counts$meanMN <- apply(abs_counts[, 5:6], 1, mean)
filtered_counts <- abs_counts[which(!is.na(abs_counts$fold_change)), ]
filtered_counts$fold_change[!is.finite(filtered_counts$fold_change)] <- -10

plot_data_final <- reshape2::melt(filtered_counts, id.vars = "Systematic.name",
                                  measure.vars = c("meanMM", "meanMN"),
                                  variable.name = "group",
                                  value.name = "Absolute Counts")
levels(plot_data_final$group) <- c("Control", "24-hr Nitrogen\nStarved")
fc_data <- merge(unique(filtered_counts[,c("Systematic.name", "fold_change")]),
                 sleuth_table[,c("gene_id", "log_fc")],
                 by.x = "Systematic.name",
                 by.y = "gene_id")
colnames(fc_data)[3] <- "sleuth_log_fc"
fc_data <- merge(fc_data, alr_table[,c("gene_id", "log_fc")],
                 by.x = "Systematic.name",
                 by.y = "gene_id")
colnames(fc_data)[4] <- "alr_log_fc"
fc_data_final <- reshape2::melt(fc_data, id.vars = "Systematic.name",
                                measure.vars = colnames(fc_data)[-1],
                                variable.name = "group",
                                value.name = "fold_change")
fc_data_final$group <- factor(fc_data_final$group,
                              levels = c("sleuth_log_fc", "alr_log_fc", "fold_change"))
levels(fc_data_final$group) <- c("sleuth", "sleuth-ALR", "Absolute\nCounts")

#levels(fc_data_final$group) <- c("Absolute\nCounts", "DESeq2")
#fc_data_final$group <- relevel(fc_data_final$group, "DESeq2")

pct_up <- sum(fc_data$fold_change > 0) / nrow(fc_data) * 100
pct_down <- sum(fc_data$fold_change < 0) / nrow(fc_data) * 100

theme_hp <- function() {
   cowplot::theme_cowplot(12) +
     ggplot2::theme(legend.key.size = unit(2, "lines"), legend.position = "none",
           axis.title.x = ggplot2::element_text(face = "bold", family = "ArialMT"),
           axis.title.y = ggplot2::element_text(face = "bold", family = "ArialMT"),
           axis.text.x = ggplot2::element_text(family = "ArialMT"),
           axis.text.y = ggplot2::element_text(family = "ArialMT"))
}

fc_violin_plot <- ggplot2::ggplot(fc_data_final,
                                    ggplot2::aes(x = group, y = fold_change, fill = group)) +
  ggplot2::geom_violin(alpha = 0.4) +
  ggplot2::scale_fill_manual(values = c("#CCCCCC", "#3f17b5", "#56B4E9")) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  theme_hp() +
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 panel.grid.major.x = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_line(size = 0.25,color = "gray"),
                 legend.position = "none") +
  ggplot2::annotate("text", hjust = 0, lineheight = 0.75,
           label = paste0(format(signif(pct_down, 3), nsmall = 1),
                          "%\nDown"),
           colour = "#56B4E9", x = 3.05, y = -7.5) +
  ggplot2::annotate("text", hjust = 0, lineheight = 0.75,
           label = paste0(format(signif(pct_up, 3), nsmall = 1),
                          "%\nUp"),
           colour = "#D55E00", x = 3.05, y = 2.75) +
  ggplot2::xlab("Method") +
  ggplot2::ylab("Log2 Fold Change")

exp_violin_plot <- ggplot2::ggplot(plot_data_final,
                                   ggplot2::aes(x = group, y = log10(`Absolute Counts`+0.01), fill = group)) +
  ggplot2::geom_violin(alpha = 0.4) +
  ggplot2::scale_fill_manual(values = c("#CC79A7","#009E73")) +
  theme_hp() +
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 panel.grid.major.x = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_line(size = 0.25, color = "gray"),
                 legend.position = "none") +
  ggplot2::xlab("Yeast Condition") +
  ggplot2::ylab("log10 Absolute Counts")

g <- cowplot::plot_grid(exp_violin_plot, fc_violin_plot, labels = "AUTO", align = "vh", nrow = 1, hjust = -1)
cowplot::save_plot('../results/figure5.pdf', g, base_aspect_ratio = 2, base_height = 3.25)
