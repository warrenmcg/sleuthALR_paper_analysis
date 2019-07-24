.libPaths(c("~/R_library", .libPaths()))
library(polyester)
library(sleuth)
library(absSimSeq)

fasta_file <- '../../annotation/Homo_sapiens.gencode.v25.fa'
sleuth_file <- '../results/finn_sleuth.rds'
outdir <- '../sims/isoform_5_5_30_645175_1'
dir.create(outdir, showWarnings = F)

### These are the high-expression spike-ins
# Two previous references observed that the ratios of
# spike-ins between Mix1 and Mix2 were accurately estimated
# for spike-ins above a certain expression threshold.
# The threshold was selected to be log2 concentration of at least
# 3 attomoles.
# The two references:
# dx.doi.org/10.1038/nbt.2957 (Figure 4C)
# dx.doi.org/10.1038/ncomms6125 (Figure 5)
data('ERCC92_data', package = "absSimSeq")
log_means <- log2(rowMeans(ERCC92_data[,c(4:5)]))
denom <- names(log_means[log_means > 3])

### We will use the mean across all Finnish samples as the starting point for the control group
sample_index <- 'mean'

### Fifteen total simulation runs (see groups A-C below)
num_runs <- 30

### Five samples each in the control and "experimental" conditions
num_reps <- c(5,5)

### To simulate more realistic data, adding in some samples with GC bias
### See the `gc_bias` option in polyester documentation for more information
gc_bias <- rep(c(0, 0, 1, 1, 1), 2)

### The parameters for the truncated normal
# 1.25 = minimum fold-change (25% up; 20% (1/1.25) down)
# 2 = mean of the truncated normal
# 2 = standard deviation of the truncated normal
de_levels <- c(1.25, 2, 2)

### The DE type to indicate use of the truncated normal
de_type <- 'normal'

### probability of being differentially expressed
# 3 groups of simulations with 5 runs in each group:
# Group A -- "small" (intended to match assumptions made by previous normalization methods)
# Group B -- "down" (intended to violate assumptions made by previous methods)
# Group C -- "up" (intended to violate assumptions made by previous methods)
# Note: Groups B and C have the same probability for differential expression (20% of all transcripts)
###
de_probs <- c(rep(0.05, 5), rep(0.2, 20), rep (0.05, 5))

### probability of being up-regulated
# Group A, "small", has all 5 runs with 30% up so that the total RNA remains about the same
# Group B, "down", has all 5 runs with 10% up (i.e. 90% down)
# Group C, "up",  has all 5 runs wtih 90% up
###
dir_probs <- c(rep(0.3, 5), rep(0.1, 5), rep(0.9, 5), rep(0.5, 10), rep(0.9, 5))

symmetries <- c(rep(FALSE, 20), rep(TRUE, 5), rep(FALSE, 5))

### Mean library size per sample
# This is set to 30 million reads
###
mean_lib_size <- 30*10^6

### The seed for the full simulation
# Used code in 'testing.R' to determine this number
###
seed <- 645175

sims <- rep(TRUE, 30)

### Launching the actual full simulation
final_results <- absSimSeq::run_abs_simulation(sleuth_file, fasta_file, sample_index,
                                               outdir = outdir, denom = denom,
                                               num_runs = num_runs, num_reps = num_reps,
                                               seed = seed, de_levels = de_levels,
                                               de_type = de_type, mean_lib_size = mean_lib_size,
                                               de_probs = de_probs, dir_probs = dir_probs,
                                               symmetries = symmetries,
                                               polyester_sim = sims, gc_bias = gc_bias,
                                               num_cores = 15)

oracles <- absSimSeq::extract_oracles(final_results)

save(final_results, sample_index, de_probs, dir_probs, fasta_file,
     sleuth_file, outdir, gc_bias, num_runs, num_reps, denom,
     seed, de_type, de_levels, symmetries, mean_lib_size, sample_index,
     file = file.path('../results', 'polyester_ground_truth.RData'))
saveRDS(oracles, file = file.path(outdir, '../polyester_oracles.rds'))
