.libPaths(c("~/R_library", .libPaths()))
args <- commandArgs(trailing = TRUE)

#' Prepare salmon for Sleuth
#' This is wrapper function for the \code{prepare_fish_for_sleuth} function
#' from the \code{wasabi} package.
#' 
#' @param in_dir, character with the relative or absolute path containing the
#'     salmon result directories
#' @param pattern, character with the regular expression to identify which
#'     directories need to be included for analysis
prepare_salmon <- function(in_dir=NULL, pattern=NULL) {
  salmon_dirs <- list.files(in_dir, pattern=pattern, full.names=T)
  wasabi::prepare_fish_for_sleuth(salmon_dirs)
}

samples <- paste("sample", sprintf('%02d', 1:10), sep = "_")
in_dirs <- file.path("../sims", args[1], args[2], samples)

prepare_salmon(in_dir = in_dirs, pattern = args[3])
