# most of these colors were pulled from the colorblind pallet
# take a look at mamabear::de_benchmark.R
method_colors <- c(DESeq2 = '#2166ac', # blue
  `DESeq2 +\nC.N.` = '#103356',
  `DESeq2 +\nRUVg` = '#5499de',
  edgeR = '#92c5de', # light blue
  `edgeR +\nC.N.` = '#4ea1ca',
  `edgeR +\nRUVg` = '#c4e0ed',
  limma = '#E69F00', # gold
  `limma +\nC.N.` = '#996900', # dark gold
  `sleuth LRT` = '#000000', # previously #999999
  `sleuth Wald` = '#999999',
  `sleuth-ALR\nLRT` = '#FF3030', # firebrick1
  `sleuth-ALR\nWald` = '#990000',
  `ALDEx2\noverlap`= '#009E73' # green
)

method_colors_old <- c('DESeq2' = "#0072B2",
  edgeR = '#CC79A7',
  limmaVoom = '#E69F00',
  ALDEx2 = '#56B4E9',
  sleuth = 'black',
  `sleuth-ALR` = 'firebrick1')

full_method_colors <- c(
  DESeq2 = '#2166ac',
  DESeq2_RUVg = '#2166ac',
  DESeq2_denom = '#2166ac',
  edgeR = '#92c5de',
  edgeR_denom = '#92c5de',
  edgeR_RUVg = '#92c5de',
  limmaVoom = '#E69F00',
  limmaVoom_denom = '#E69F00',
  ALDEx2.iqlr.overlap = '#009E73',
  ALDEx2.iqlr.wilcoxon = '#CC79A7',
  ALDEx2.iqlr.welch = '#FFC64D',
  ALDEx2.clr.overlap = '#009E73',
  ALDEx2.clr.wilcoxon = '#CC79A7',
  ALDEx2.clr.welch = '#FFC64D',
  ALDEx2.denom.overlap = '#009E73',
  ALDEx2.denom.wilcoxon = '#CC79A7',
  ALDEx2.denom.welch = '#FFC64D',
  sleuth.lrt = '#000000',
  sleuth.wt = '#000000',
  sleuthALR.lrt = '#FF3030',
  sleuthALR.wt = '#990000', # dark red
  sleuthALR.counts.lrt = '#FF66FF', # light purple
  sleuthALR.counts.wt = '#660066') # dark purple

method_kal_ltys <- rep(1, length(full_method_colors))
method_sal_ltys <- rep(2, length(full_method_colors))
names(method_kal_ltys) <- names(method_sal_ltys) <- names(full_method_colors)

#' Load the differential expression truth
#'
#' This function loads the isoform level information and gene level information.
#' The gene level differential expression is defined by at least 1 isoform of the gene being differentially expressed.
#'
#' @param sim_name the name of the simulation (with respect to the results location)
#' @param ttg the 'target_to_gene_mapping' coming from ensembl
#' @return a list containing the gene and isoform differential expression information
get_de_info <- function(sim_name, sim_number, ttg) {
  sims_dir <- file.path('../sims', sim_name)

  # de_info contains the 'true' fold change as well as which transcripts are DE
  all_sim <- readRDS(file.path(sims_dir, 'sims.rds'))
  de_info <- all_sim[[sim_number]]$info

  if (length(grep("_gene", colnames(de_info))) > 0) {
    tmp <- de_info
  } else {
    tmp <- dplyr::inner_join(de_info, ttg, by = 'target_id')
  }

  de_genes <- dplyr::select(tmp, ens_gene, is_de)
  de_genes <- dplyr::group_by(de_genes, ens_gene)
  de_genes <- dplyr::summarize(de_genes, is_de = any(is_de))

  de_genes <- dplyr::rename(de_genes, target_id = ens_gene)
  de_genes <- dplyr::mutate(de_genes, log_fc = NA)

  list(de_info = de_info, de_genes = de_genes)
}

#' DEPRECATED: use `load_union_counts` below
#' Get gene counts from featureCounts
#'
#' Assumes \code{featureCounts} has been run and the results are in the appropriate place.
#'
#' @param sim_name the name of the simulation (with respect to the results location)
#' @param sim_num the integer number of the simulation (with respect to \code{sim_name})
#' @return a matrix containing gene counts
get_gene_counts <- function(sim_name, sim_num = 1) {
  base_dir <- file.path('../results', sim_name)

  experiment_name <- paste0('exp_', sim_num)
  counts <- read.table(file.path(base_dir, experiment_name, 'counts.tsv'),
    header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  counts <- as.matrix(counts)
  colnames(counts) <- sub('X', 'sample_', colnames(counts))

  counts
}

#' Get gene counts from featureCounts
#'
#' Assumes \code{featureCounts} has been run and the results are in the appropriate place.
#'
#' @param sim_info a `sim_info` object resulting from `parse_simulation`
#' @param which_sample the integer number of the simulation (with respect to \code{sim_info$name})
#' @return a matrix containing gene counts
load_union_counts <- function(sim_info, which_sample) {
  which_sample <- as.character(as.integer(which_sample))
  path <- file.path('..', 'sims', sim_info$name, paste0('exp_', which_sample),
    'results', 'featureCounts', 'counts.tsv')

  counts <- data.table::fread(path, sep = '\t', data.table = FALSE, header = TRUE)
  rownames(counts) <- counts$gene_id
  counts$gene_id <- NULL
  setnames(counts, paste0('sample_', colnames(counts)))
  colnames(counts) <- sub('X', 'sample_', colnames(counts))
  counts <- as.matrix(counts)

  counts
}

readFeatureCounts <- function(path) {
  suppressWarnings( result <- data.table::fread(path, data.table = TRUE) )

  # pull out the columns 'gene_id' and 'count'
  result <- result[, c(1, 7), with = FALSE]
  data.table::setnames(result, c('gene_id', 'count'))

  result
}

load_union_counts_general <- function(file_paths, sample_names) {
  all_counts <- lapply(seq_along(file_paths),
    function(i) {
      counts <- readFeatureCounts(file_paths[i])
      data.table::setnames(counts, c('count'), sample_names[i])
    })
  counts_table <- Reduce(
    function(x, y) dplyr::inner_join(x, y, by = 'gene_id'),
    all_counts
    )

  counts_table <- data.frame(counts_table, stringsAsFactors = FALSE)
  rownames(counts_table) <- counts_table$gene_id
  counts_table$gene_id <- NULL

  counts_table <- as.matrix(counts_table)
  storage.mode(counts_table) <- 'integer'

  counts_table
}

readKalTpms <- function(path) {
  suppressWarnings( result <- data.table::fread(path, data.table = TRUE) )

  # pull out the columns 'target_id' and 'tpm'
  result <- result[, c(1, 5), with = FALSE]
  data.table::setnames(result, c('target_id', 'tpm'))

  result
}

load_union_tpms_general <- function(file_paths, sample_names, ttg) {
  ttg <- data.table::as.data.table(ttg)

  all_tpms <- lapply(seq_along(file_paths),
    function(i) {
      tpms <- readKalTpms(file_paths[i])
      tpms$sample <- sample_names[i]
      tmp <- merge(tpms, ttg, by.x = "target_id", by.y = "target_id_backup")
      tpm_gene <- tmp[, j = list(tpm = sum(tpm)),
                      by = list(ens_gene)]
      data.table::setnames(tpm_gene, c('gene_id', sample_names[i]))
      tpm_gene
    })
  tpm_table <- Reduce(
    function(x, y) dplyr::inner_join(x, y, by = 'gene_id'),
    all_tpms
    )

  tpm_table <- data.frame(tpm_table, stringsAsFactors = FALSE)
  rownames(tpm_table) <- tpm_table$gene_id
  tpm_table$gene_id <- NULL

  tpm_table <- as.matrix(tpm_table)

  tpm_table
}
