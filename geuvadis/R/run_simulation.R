.libPaths(c("~/R_library", .libPaths()))
fasta_file <- '../../annotation/Homo_sapiens.gencode.v25.fa'
sleuth_file <- '../results/finn_sleuth.rds'
outdir <- '../sims/isoform_5_5_15_387_1'
dir.create(outdir, showWarnings = F)

### This denominator was chosen for two reasons:
# 1) It is known to not be changing across all fifteen runs
# 2) Among the other genes that are not changing across all runs, 
#    it had the lowest average COV
###
denom <- 'ENST00000377532.7' # PER3-202

### This sample is the most deeply sequenced Finnish sample in the GEUVADIS dataset
sample_index <- 'mean'

### Fifteen total simulation runs (see groups A-C below)
num_runs <- 15

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
# Group A -- "few changes" (intended to match assumptions made by previous normalization methods)
# Group B -- "many changes down" (intended to violate assumptions made by previous methods)
# Group C -- "many changes up" (intended to violate assumptions made by previous methods)
# Note: Groups B and C have the same probability
###
de_probs <- c(rep(0.05, 5), rep(0.20, 10))

### probability of being up-regulated
# Group A has 3 runs of 50/50, and one run each with a skew of 90% down or up, respectively
# Group B has all 5 runs with 90% down
# Group C has all 5 runs wtih 90% up
###
dir_probs <- c(rep(0.3, 5), rep(0.1, 5), rep(0.9, 5))

### Mean library size per sample
# This is set to 25 million reads
###
mean_lib_size <- 30*10^6

### The seed for the full simulation
# Used code in 'testing.R' to determine this number
###
seed <- 645175

### Launching the actual full simulation
# On our machine, we need to use 8 cores because of the memory demands
# 8 cores consumes about 100-110 GB of RAM (out of 128 on our machine)
final_results <- absSimSeq::run_abs_simulation(fasta_file, sleuth_file, sample_index,
                                               outdir = outdir, denom = denom,
                                               num_runs = num_runs, num_reps = num_reps,
                                               seed = seed, de_levels = de_levels,
                                               de_type = de_type, mean_lib_size = mean_lib_size,
                                               de_probs = de_probs, dir_probs = dir_probs,
                                               polyester_sim = T, gc_bias = gc_bias,
                                               sleuth_save = T, num_cores = 8)

oracles <- absSimSeq::extract_oracles(final_results)

save(final_results, sample_index, de_probs, dir_probs, fasta_file,
     sleuth_file, outdir, gc_bias, num_runs, num_reps, denom,
     seed, de_type, de_levels, mean_lib_size, sample_index,
     file = file.path('../results', 'polyester_ground_truth.RData'))
saveRDS(oracles, file = file.path(outdir, '../polyester_oracles.rds'))
