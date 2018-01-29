# This file contains functions to run several different differential expression
# methods, usually simply by providing the "count matrix"

library("data.table")

library("Biobase")
library("DESeq2")
library("edgeR")
library("limma")
library("sleuth")
library("sleuthALR")
library("ALDEx2")

get_human_gene_names <- function() {
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host = "dec2016.archive.ensembl.org")
  ttg <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "transcript_version",
      "ensembl_gene_id", "version", "external_gene_name"),
    mart = mart)
  ttg <- dplyr::mutate(ttg, target_id = paste(ensembl_transcript_id, transcript_version, sep = "."),
    ens_gene = paste(ensembl_gene_id, version, sep = "."))
  ttg <- dplyr::rename(ttg, ext_gene = external_gene_name)
  ttg <- dplyr::select(ttg, target_id, ens_gene, ext_gene)
  ttg
}

get_mouse_gene_names <- function() {
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "mmusculus_gene_ensembl",
    host = "dec2016.archive.ensembl.org")
  ttg <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "transcript_version",
    "ensembl_gene_id", "external_gene_name", "version"),
    mart = mart)
  ttg <- dplyr::rename(ttg, ext_gene = external_gene_name)
  ttg <- dplyr::mutate(ttg, ens_gene = paste(ensembl_gene_id, version, sep="."),
    target_id = paste(ensembl_transcript_id, transcript_version, sep="."))
  ttg <- dplyr::select(ttg, target_id, ens_gene, ext_gene)
  ttg
}


#' parse the simulation name
#'
#' Get back a named list with the simulation info
#'
#' @param sim_name a simulation name in the form 'isoform_3_3_20_1_1'
parse_simulation <- function(sim_name) {
  sim_name_split <- strsplit(sim_name, '_')[[1]]

  res <- list()
  res[['type']] <- sim_name_split[[1]]
  res[['a']] <- as.integer(sim_name_split[[2]])
  res[['b']] <- as.integer(sim_name_split[[3]])
  res[['n']] <- as.integer(sim_name_split[[4]])
  res[['seed']] <- as.integer(sim_name_split[[5]])
  res[['sf']] <- as.integer(sim_name_split[[6]])
  res[['name']] <- sim_name

  res
}

sleuth_filter <- function(mat, ...) {
  apply(mat, 1, sleuth::basic_filter, ...)
}

edgeR_filter <- function(mat, ...) {
  rowSums(cpm(mat) > 1) >= 3
}

DESeq2_filter <- function(mat, ...) {
  rowSums(mat) > 1
}

#' create benchmark objects against a reference
#'
#' Given several differential expression tables, compare them against a
#' reference in order to make benchmark objects using \code{new_de_benchmark}
#' @param results_list a named list with a set of results (differential expression tables)
#' @param de_info the "truth" that was used to simulate differential expression
#' @param reference a character string vector with a set of references
compare_reference <- function(results_list, de_info,
  reference = c('sleuth.wt', 'sleuth.lrt')) {

  other_methods <- names(results_list)[!(names(results_list) %in% reference)]
  if (length(other_methods) == 0) {
    stop('No other methods (or missing names) in results_list')
  }

  res <- lapply(other_methods,
    function(m) {
      ms <- c(reference, m)
      new_de_benchmark(results_list[ms], ms, de_info)
    })
  names(res) <- other_methods

  res
}

#' @param sim_name a simulation name such as 'isoform_3_3_20_1_1'
#' @param which_sample which sample (replication) to load (an integer from 1 to N)
#' @param method_filtering if \code{TRUE}, use the methods own filtering.
#' Otherwise, use the filtering provided from sleuth.
#' @param ... additional arguments passed to \code{run_sleuth}
#' NOTE: cuffdiff uses a filter method internally and it cannot be changed
load_isoform_results_intersect <- function(
  sim_name,
  which_sample,
  tool,
  method_label,
  method_fit_function,
  ...) {
  sim <- parse_simulation(sim_name)

  tool <- match.arg(tool, c("kallisto", "salmon"))
  if (which_sample > sim$n) {
    stop('which_sample must be less than the total number of replications: ', sim$n)
  }

  which_sample <- sprintf('%02d', which_sample)
  n <- sim$a + sim$b

  kal_dirs <- file.path('..', 'sims', sim_name, paste0("run", which_sample),
    paste0("sample_", sprintf('%02d', 1:n)), tool)
  sample_to_condition <- get_sample_to_condition(sim$a, sim$b, kal_dirs)

  if (method_label == 'sleuth') {
    sir <- run_sleuth(sample_to_condition, gene_mode = NULL, ...)
    all_results <- list(sleuth.lrt = sir$sleuth.lrt, sleuth.wt = sir$sleuth.wt)
    all_results
  } else if (method_label == 'sleuthALR') {
    all_results <- run_alr(sample_to_condition, gene_mode = NULL, ...)
    all_results
  } else if (method_label == 'ALDEx2') {
    # simply do this so we can read in the data
    so_data <- sleuth_prep(sample_to_condition, ~1, max_bootstrap = 2)
    obs_raw <- sleuth:::spread_abundance_by(so_data$obs_raw, "est_counts")
    rm(so_data)

    s_which_filter <- sleuth_filter(obs_raw)

    res <- runALDEx2(obs_raw, stc = sample_to_condition, match_filter = s_which_filter,
                     which_test = 'all', ...)
    all_results <- list(ALDEx2.overlap = res$overlap, ALDEx2.welsh = res$welsh,
                        ALDEx2.wilcoxon = res$wilcoxon)
    all_results
  } else {
    # simply do this so we can read in the data
    so_data <- sleuth_prep(sample_to_condition, ~1, max_bootstrap = 2)
    obs_raw <- sleuth:::spread_abundance_by(so_data$obs_raw, "est_counts")
    rm(so_data)

    s_which_filter <- sleuth_filter(obs_raw, ...)

    method_result <- method_fit_function(obs_raw, sample_to_condition,
      s_which_filter)

    all_results[[method_label]] <- method_result$result

    all_results
  }
}

###
# TODO: deprecate load_isoform_results_intersect and replace with
# load_isoform_results_intersect_df
###
#' @param sim_name a simulation name such as 'isoform_3_3_20_1_1'
#' @param which_sample which sample (replication) to load (an integer from 1 to N)
#' @param method_filtering if \code{TRUE}, use the methods own filtering.
#' Otherwise, use the filtering provided from sleuth.
#' @param ... additional arguments passed to \code{run_sleuth}
#' NOTE: cuffdiff uses a filter method internally and it cannot be changed
load_isoform_results_intersect_df <- function(
  sample_to_condition,
  method_label,
  method_fit_function,
  ...) {

  sample_to_condition <- as.data.frame(sample_to_condition,
    stringsAsFactors = FALSE)
  rownames(sample_to_condition) <- sample_to_condition$sample

  message('### Loading data with sleuth...')
  so_data <- sleuth_prep(sample_to_condition, ~1, max_bootstrap = 3)
  obs_raw <- sleuth:::spread_abundance_by(so_data$obs_raw, "est_counts")
  rm(so_data)

  obs_raw <- obs_raw[, rownames(sample_to_condition)]

  s_which_filter <- sleuth_filter(obs_raw, ...)

  message('### Running method: ', method_label)
  method_result <- method_fit_function(obs_raw, sample_to_condition,
    s_which_filter)
  sir <- run_sleuth(sample_to_condition, gene_mode = NULL)

  all_results <- Filter(is.data.frame, sir)

  message('### Running sleuth...')

  all_results[[method_label]] <- method_result$result

  all_results
}

load_gene_results_intersect <- function(
  sim_name,
  which_sample,
  method_label,
  method_fit_function,
  ...) {
  sim <- parse_simulation(sim_name)

  if (which_sample > sim$n) {
    stop('which_sample must be less than the total number of replications: ', sim$n)
  }

  which_sample <- as.integer(which_sample)
  n <- sim$a + sim$b

  sim_info <- get_de_info(sim_name, which_sample, transcript_gene_mapping)
  de_info <- sim_info$de_info
  de_genes <- sim_info$de_genes

  kal_dirs <- file.path('..', 'sims', sim_name, paste0("exp_", which_sample),
    1:n, "kallisto")
  sample_to_condition <- get_sample_to_condition(sim$a, sim$b, kal_dirs)

  # TODO: load the gene counts
  counts <- load_union_counts(sim, which_sample)

  # simply do this so we can read in the data
  so_data <- sleuth_prep(sample_to_condition, ~1, max_bootstrap = 3)
  obs_raw <- sleuth:::spread_abundance_by(so_data$obs_raw, "est_counts")

  # s_which_filter <- sleuth_filter(obs_raw, ...)
  tmp <- so_data$obs_raw
  tmp <- dplyr::group_by(tmp, target_id)
  tmp <- dplyr::summarize(tmp, pass_filter = sleuth::basic_filter(est_counts))
  tmp <- dplyr::inner_join(tmp, transcript_gene_mapping, by = 'target_id')
  sleuth_gene_filter <- dplyr::filter(tmp, pass_filter)
  sleuth_gene_filter <- dplyr::select(sleuth_gene_filter, target_id = ens_gene)
  sleuth_gene_filter <- dplyr::distinct(sleuth_gene_filter)

  s_which_filter <- rownames(counts) %in% sleuth_gene_filter$target_id
  names(s_which_filter) <- rownames(counts)

  message(paste0('### running method: ', method_label))
  method_result <- method_fit_function(counts, sample_to_condition,
    s_which_filter)
  # since the gene filter is derived from the isoform filter, just use the sleuth
  # filter for gene lifting
  # debugonce(sleuth_prep)

  message('### running gene lifting')
  slr <- run_sleuth(sample_to_condition, gene_mode = 'lift')
  slr <- Filter(is.data.frame, slr)
  names(slr) <- paste0(names(slr), '.lift')

  message('### running gene aggregation')
  sar <- run_sleuth(sample_to_condition, gene_mode = 'aggregate',
    gene_column = 'ens_gene')

  # TODO: adjust the fdr based off of the filtering scheme
  # e.g. take the intersection of the tests and recompute the fdr
  # This ensures that the calibration tests can be comparable

  slr <- Filter(is.data.frame, slr)
  # names(slr) <- paste0(names(slr), '.lift')
  sar <- Filter(is.data.frame, sar)
  names(sar) <- paste0(names(sar), '.agg')

  all_results <- c(slr, sar)
  # all_results <- sar
  all_results[[method_label]] <- method_result$result

  all_results
}


###
# these functions are used for the isoform level analysis
###
limma_filter_and_run <- function(counts, stc, match_filter) {
  which_targets <- DESeq2_filter(counts)
  match_filter <- match_filter & which_targets
  cds <- make_count_data_set(counts[match_filter, ], stc)

  res <- runVoom(cds, FALSE, FALSE)
  match_filter <- names(which(match_filter))
  list(result = res, filter = match_filter)
}

# NEW METHOD FOR ALDEx2
aldex2_filter_and_run <- function(counts, denom, stc, match_filter, which_test) {
  which_targets <- sleuth_filter(counts)
  match_filter <- match_filter & which_targets
  res <- runALDEx2(counts, stc$condition, denom = denom, FALSE, "t", which_test)
  match_filter <- names(which(match_filter))
  if (which_test == "all") {
    list(aldex2_overlap = res$overlap,
         aldex2_welsh = res$welsh,
         aldex2_wilcoxon = res$wilcoxon,
         filter = match_filter)
  } else {
    list(result = res, filter = match_filter)
  }
}

# DEPRECATED
# DESeq2_filter_and_run <- function(count_matrix, stc, sleuth_filter) {
#   # we should check if taking the intersection of results does better
#   count_matrix <- round(count_matrix)
#   mode(count_matrix) <- 'integer'
#   which_targets <- DESeq2_filter(count_matrix)
#   cds <- make_count_data_set(count_matrix[which_targets, ], stc)
#   res <- runDESeq2(cds, FALSE, FALSE)
#
#   sleuth_filter <- sleuth_filter & which_targets
#   sleuth_filter <- names(which(sleuth_filter))
#
#   list(result = res, filter = sleuth_filter)
# }

DESeq2_filter_and_run_intersect <- function(counts, stc, match_filter, # nolint
  is_counts = TRUE) {
  if (is_counts) {
    counts <- round(counts)
    mode(counts) <- 'integer'
    which_targets <- DESeq2_filter(counts)
    match_filter <- match_filter & which_targets
    cds <- make_count_data_set(counts[match_filter, ], stc)
  } else {
    cds <- DESeqDataSetFromTximport(counts, stc, ~condition)
    which_targets <- DESeq2_filter(counts(cds))
    match_filter <- match_filter & which_targets
    cds <- cds[match_filter, ]
  }
  res <- runDESeq2(cds, FALSE, FALSE, is_counts)

  match_filter <- names(which(match_filter))

  list(result = res, filter = match_filter)
}

edgeR_filter_and_run <- function(counts, stc, match_filter, is_counts = TRUE) {
  if (is_counts) {
    counts <- round(counts)
    mode(counts) <- 'integer'
    which_targets <- edgeR_filter(counts)
    match_filter <- match_filter & which_targets
    cds <- make_count_data_set(counts[match_filter, ], stc)
    design <- NULL
  } else {
    txi <- counts
    # below boilerplate taken from tximport vignette and modified to include filtering
    cts <- txi$counts
    which_targets <- edgeR_filter(cts)
    match_filter <- match_filter & which_targets

    cts <- cts[match_filter, ]
    normMat <- txi$length[match_filter, ]
    normMat <- normMat/exp(rowMeans(log(normMat)))
    o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))

    y <- DGEList(cts)
    y$offset <- t(t(log(normMat)) + o)
    cds <- y
    # y is now ready for estimate dispersion functions see edgeR User's Guide
    design <- model.matrix(~condition, stc)
    colnames(design)[2] <- "pData(e)$conditionB"
  }
  res <- runEdgeR(cds, FALSE, FALSE, is_counts, design)

  match_filter <- names(which(match_filter))

  list(result = res, filter = match_filter)
}

#' generate `sample_to_condition` that sleuth is expecting
#'
#' @param n_a the number of samples in condition A
#' @param n_a the number of samples in condition B
#' @param kal_dirs if not NULL, then add the appropriate `path` column
#' @return a \code{data.frame} in the proper sleuth form
get_sample_to_condition <- function(n_a, n_b, kal_dirs = NULL) {
  n <- n_a + n_b

  sample_to_condition <- data.frame(
    sample = paste0("sample_", sprintf('%02d', 1:n)),
    condition = factor(c(rep("A", n_a), rep("B", n_b))),
    stringsAsFactors = FALSE)

  if (!is.null(kal_dirs)) {
    sample_to_condition <- dplyr::mutate(sample_to_condition, path = kal_dirs)
  }
  rownames(sample_to_condition) <- sample_to_condition$sample

  sample_to_condition
}

#' Generate equal size factors
#'
#' Generate size factors all equal to 1. Helpful for simulations.
#'
#' @param x the count matrix
#' @return a vector of all ones
all_ones <- function(x) {
  p <- ncol(x)
  sf <- rep.int(1, p)
  names(sf) <- colnames(x)

  sf
}

#' @param counts \code{matrix} of counts with transcripts on the rows and
#' samples on the columns
#' @param sample_info \code{data.frame} of sample information with at least a
#' column called \code{condition}
make_count_data_set <- function(counts, sample_info) {
  ExpressionSet(counts, AnnotatedDataFrame(sample_info))
}

run_sleuth_prep <- function(sample_info, max_bootstrap = 30, gene_column = NULL,
  ...) {
  so <- sleuth_prep(sample_info, ~ condition, max_bootstrap = max_bootstrap,
    target_mapping = transcript_gene_mapping,
    ...)
  so <- sleuth_fit(so)

  so
}

# new method to run sleuthALR
run_alr <- function(sample_info,
  max_bootstrap = 30,
  gene_column = NULL,
  denom = NULL,
  ...) {
  so <- sleuthALR::make_lr_sleuth_object(sample_info,
    target_mapping = transcript_gene_mapping,
    beta = 'conditionB',
    denom_name = denom, aggregate_column = gene_column,
    max_bootstrap = 30,
    ...)
  lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt',
    show_all = FALSE)[, c('target_id', 'pval', 'qval')]
  wt <- sleuth_results(so, 'conditionB',
    show_all = FALSE)[, c('target_id', 'pval', 'qval')]
  res <- list(sleuthALR.lrt = lrt, sleuthALR.wt = wt)
  res
}

#' @param gene_mode if NULL, do isoform mode, if 'lift' do gene lifting, if 'aggregate', do gene aggregation
run_sleuth <- function(sample_info,
  max_bootstrap = 30,
  gene_mode = NULL,
  gene_column = NULL,
  ...) {

  so <- NULL
  if (!is.null(gene_column)) {
    so <- run_sleuth_prep(sample_info, max_bootstrap = max_bootstrap,
       aggregation_column = gene_column,
      ...)
  } else {
  so <- run_sleuth_prep(sample_info, max_bootstrap = max_bootstrap,
    ...)
  }
  so <- sleuth_wt(so, 'conditionB')
  so <- sleuth_fit(so, ~ 1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')

  res <- NULL
  if (is.null(gene_mode)) {
    lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt',
      show_all = FALSE)[, c('target_id', 'pval', 'qval')]
    wt <- sleuth_results(so, 'conditionB',
      show_all = FALSE)[, c('target_id', 'pval', 'qval')]
    res <- list(sleuth.lrt = lrt, sleuth.wt = wt)
  } else if (gene_mode == 'lift') {
    # test every isoform
    lrt <- get_gene_lift(so, 'reduced:full', test_type = 'lrt')
    wt <- get_gene_lift(so, 'conditionB', test_type = 'wt')
    res <- list(sleuth.lrt = lrt, sleuth.wt = wt)
  } else if (gene_mode == 'aggregate') {
    lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt',
      show_all = FALSE)[, c('target_id', 'pval', 'qval')]
    wt <- sleuth_results(so, 'conditionB',
      show_all = FALSE)[, c('target_id', 'pval', 'qval')]
    res <- list(sleuth.lrt = lrt, sleuth.wt = wt)
  } else {
    stop('Unrecognized mode for "run_sleuth"')
  }

  res$so <- so

  res
}

# if method_filtering is false, then use the sleuth filter
get_filtered_isoform_cds <- function(so, stc, method_filtering = FALSE) {
  pass_filt_names <- so$filter_df[['target_id']]

  obs_raw <- sleuth:::spread_abundance_by(so$obs_raw, "est_counts")
  if (!method_filtering) {
    obs_raw <- obs_raw[pass_filt_names,]
  }
  isoform_cds <- make_count_data_set(round(obs_raw), stc)

  isoform_cds
}

#' @param obj a sleuth object
#' @param ... additional arguments to \code{sleuth_gene_table}
#' @example \code{get_gene_lift(so, 'reduced:full', test_type = 'lrt')}
get_gene_lift <- function(obj, ...) {
  sgt <- sleuth_gene_table(obj, ...)
  sgt <- group_by(sgt, ens_gene)
  do(sgt, {
    min_index <- which.min(.$qval)
    result <- .[min_index, ]
    result <- dplyr::select(result, -target_id)
    result <- dplyr::rename(result, target_id = ens_gene)
    # dplyr::select(result, target_id, pval, qval)
    result
  })
}

rename_target_id <- function(df, as_gene = FALSE) {
  if (as_gene) {
    dplyr::rename(df, ens_gene = target_id)
  } else {
    df
  }
}

## NEW METHOD: Run ALDEx2
runALDEx2 <- function(counts, conditions = NULL, denom = "all", as_gene = TRUE, test = "t", statistic = "welsh") {
  mode(counts) <- "integer"
  counts_mat <- as.data.frame(counts)
  x <- ALDEx2::aldex.clr(reads = counts_mat, conds = conditions, mc.samples = 128, denom = denom,
                 verbose = TRUE)
  if (test == "t") {
    message('aldex.ttest: doing t-test')
    x_tt <- ALDEx2::aldex.ttest(x, conditions, paired.test = FALSE)
  } else if (test == "glm") {
    message('aldex.glm: doing Kruskal Wallace and glm test')
    x_tt <- ALDEx2::aldex.glm(x, conditions)
  }
  message('aldex.effect: calculating effect sizes')
  x_effect <- ALDEx2::aldex.effect(x, conditions, include.sample.summary = TRUE,
                                   verbose = TRUE)#, useMC = TRUE)
  result_df <- data.frame(x_effect, x_tt)
  result_df <- result_df[order(result_df$we.eBH, result_df$we.ep, result_df$overlap), ]
  statistic <- match.arg(statistic, c("all", "welsh", "wilcoxon", "overlap"))
  if (statistic == "welsh") {
    rename_target_id(
      data.frame(target_id = rownames(result_df),
        dplyr::select(result_df, pval = we.ep,
          qval = we.eBH,
          nonpara_pval = wi.ep,
          nonpara_qval = wi.eBH,
          overlap = overlap,
          effect = effect),
        stringsAsFactors = FALSE),
      as_gene = as_gene)
  } else if (statistic == "wilcoxon") {
    rename_target_id(
      data.frame(target_id = rownames(result_df),
        dplyr::select(result_df, pval = wi.ep,
          qval = wi.eBH,
          para_pval = we.ep,
          para_qval = we.eBH,
          overlap = overlap,
          effect = effect),
        stringsAsFactors = FALSE),
      as_gene = as_gene)
  } else if (statistic == "overlap") {
    result <- rename_target_id(
      data.frame(target_id = rownames(result_df),
        dplyr::select(result_df, qval = overlap,
          para_pval = we.ep,
          para_qval = we.eBH,
          nonpara_pval = wi.ep,
          nonpara_qval = wi.eBH,
          effect = effect),
        stringsAsFactors = FALSE),
      as_gene = as_gene)
    result$pval <- result_df$overlap
    result
  } else {
    result <- list()
    result$welsh <- rename_target_id(
      data.frame(target_id = rownames(result_df),
        dplyr::select(result_df, pval = we.ep,
          qval = we.eBH, effect = effect),
        stringsAsFactors = FALSE),
      as_gene = as_gene)
    result$wilcoxon <- rename_target_id(
      data.frame(target_id = rownames(result_df),
        dplyr::select(result_df, pval = wi.ep,
          qval = wi.eBH, effect = effect),
        stringsAsFactors = FALSE),
      as_gene = as_gene)
    result$overlap <- rename_target_id(
      data.frame(target_id = rownames(result_df),
        dplyr::select(result_df, qval = overlap,
          effect = effect),
        stringsAsFactors = FALSE),
      as_gene = as_gene)
    result$overlap$pval <- result_df$overlap
    result
  }
}

# The code below is a slightly modified version of the code from `DESeq2paper`
# http://www-huber.embl.de/DESeq2paper/
runDESeq2 <- function(e, as_gene = TRUE, compute_filter = FALSE, is_counts = TRUE) {

  dds <- NULL
  if (is_counts) {
    dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  } else {
    dds <- e
  }

  if (compute_filter) {
    # Section 1.3.6 in DESeq2 vignette
    # https://www.bioconductor.org/packages/3.3/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
    dds <- dds[ rowSums(counts(dds)) > 1, ]
  }
  dds <- DESeq(dds,quiet=TRUE)
  res <- results(dds)
  beta <- res$log2FoldChange
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  # pvals[rowSums(exprs(e)) == 0] <- NA
  padj[is.na(padj)] <- 1

  rename_target_id(
    data.frame(target_id = rownames(res),
      pval = pvals, qval = padj, beta = beta,
      stringsAsFactors = FALSE),
    as_gene = as_gene)
}

runEdgeR <- function(e, as_gene = TRUE, compute_filter = FALSE, is_counts = TRUE, design = NULL) {
  if (is_counts) {
    design <- model.matrix(~ pData(e)$condition)
    dgel <- DGEList(exprs(e))
  } else {
    dgel <- e
  }
  if (compute_filter) {
    # Section 2.6 in edgeR vignette
    # https://www.bioconductor.org/packages/3.3/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
    keep <- rowSums(cpm(dgel) > 1) >= 2
    dgel <- dgel[keep, , keep.lib.sizes=FALSE]
  }
  dgel <- calcNormFactors(dgel)
  dgel <- estimateGLMCommonDisp(dgel, design)
  dgel <- estimateGLMTrendedDisp(dgel, design)
  dgel <- estimateGLMTagwiseDisp(dgel, design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  # predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  # predbeta10 <- predFC(exprs(e), design, prior.count=10, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  predbeta <- predFC(dgel$counts, design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  predbeta10 <- predFC(dgel$counts, design, prior.count=10, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  # pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1

  rename_target_id(
    data.frame(
      target_id = rownames(edger.lrt$table),
      pval = pvals,
      qval = padj,
      beta = log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
      predbeta = predbeta[,"pData(e)$conditionB"],
      predbeta10 = predbeta10[,"pData(e)$conditionB"],
      stringsAsFactors = FALSE),
    as_gene = as_gene)
}

runEdgeRRobust <- function(e, as_gene = TRUE) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  # settings for robust from robinson_lab/edgeR_robust/robust_simulation.R
  dgel <- estimateGLMRobustDisp(dgel, design, maxit=6)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  # list(pvals=pvals, padj=padj, beta=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
  #      predbeta=predbeta[,"pData(e)$conditionB"])
  rename_target_id(
    data.frame(
      target_id = rownames(edger.lrt$table),
      pval = pvals,
      qval = padj,
      beta=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
      predbeta=predbeta[,"pData(e)$conditionB"],
      stringsAsFactors = FALSE),
    as_gene = as_gene)
}

runVoom <- function(e, as_gene = TRUE, compute_filter = FALSE) {
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  if (compute_filter) {
    # Section 2.6 in edgeR vignette
    # https://www.bioconductor.org/packages/3.3/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
    keep <- rowSums(cpm(dgel) > 1) >= 2
    dgel <- dgel[keep, , keep.lib.sizes=FALSE]
  }
  dgel <- calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  pvals <- tt$P.Value
  # pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1

  rename_target_id(data.frame(target_id = rownames(tt),
    pval = pvals,
    qval = padj,
    beta = tt$logFC,
    stringsAsFactors = FALSE),
  as_gene = as_gene)
}

# these methods used right now and need to be updated

runDSS <- function(e) {
  X <- as.matrix(exprs(e))
  colnames(X) <- NULL
  designs <- as.character(pData(e)$condition)
  seqData <- newSeqCountSet(X, designs)
  seqData <- estNormFactors(seqData)
  seqData <- estDispersion(seqData)
  result <- waldTest(seqData, "B", "A")
  result <- result[match(rownames(seqData),rownames(result)),]
  pvals <- result$pval
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=( log2(exp(1)) * result$lfc ))
}

runDESeq2Outliers <- function(e, retDDS=FALSE) {
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  ddsDefault <- DESeq(dds, quiet=TRUE)
  ddsNoRepl <- ddsDefault
  if (ncol(e) >= 14) {
    # insert original maximum Cook's distances
    # so the rows with replacement will be filtered
    # this avoid re-running with minReplicateForReplace=Inf
    mcols(ddsNoRepl)$maxCooks <- apply(assays(ddsNoRepl)[["cooks"]], 1, max)
  }
  resDefault <- results(ddsDefault)
  resNoFilt <- results(ddsDefault, cooksCutoff=FALSE)
  resNoRepl <- results(ddsNoRepl)
  resList <- list("DESeq2"=resDefault, "DESeq2-noFilt"=resNoFilt, "DESeq2-noRepl"=resNoRepl)
  resOut <- lapply(resList, function(res) {
    pvals <- res$pvalue
    padj <- res$padj
    pvals[is.na(pvals)] <- 1
    pvals[rowSums(exprs(e)) == 0] <- NA
    padj <- p.adjust(pvals,method="BH")
    padj[is.na(padj)] <- 1
    list(pvals=pvals, padj=padj)
  })
  return(resOut)
}
