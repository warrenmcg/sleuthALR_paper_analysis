# This script is to move the files from the run directory "sims/[exp]/[run]"
# to the sample directory "sims/[exp]/[run]/[sample]".
# This conforms with the rest of the simulation pipeline, which expects the
# simulated data to be in sample directories

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop('Usage: Rscript move_sim_files.R DIRECTORY')
}

input_dir <- args[1]
current_dir <- getwd()
setwd(input_dir)
dirs <- list.files('.', pattern = 'run')
for (directory in dirs) {
  print(directory)
  file_names <- list.files(directory, 'fasta.gz')
  print(file_names)
  sample_names <- unique(gsub('_[12].fasta.gz', '', basename(file_names)))
  dummy <- lapply(sample_names, function(name) {
    dir.create(file.path(directory, name), showWarnings = F)
    files <- file_names[grepl(name, file_names)]
    file.rename(file.path(directory, files), file.path(directory, name, files))
  })
}
setwd(current_dir)

message("The script for renaming the simulated data files is complete")
