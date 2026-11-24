# =============================================================================
# DiffBind Normalization Comparison: Simulation Analysis
# =============================================================================
# Description: Applies five normalization strategies across randomized sample
#              combinations (negative control experiment) and records the number
#              of differentially bound regions detected under each method.
# Input:       Randomized sample information list (sample_information_randn.rds)
# Output:      DiffBind analysis workspace and differential region counts table
# =============================================================================

library(DiffBind)
library(tidyverse)

rm(list = ls())

# --- Load data ---------------------------------------------------------------

sample_information_randn <- readRDS(
  file.path(working_dir,
            "R_analysis/02.Sample_information/data/sample_information_randn.rds")
)

# --- Initialise output lists -------------------------------------------------

nruns <- length(sample_information_randn)

dbObj_first_list          <- vector("list", nruns)
dbObj_counts_summit_list  <- vector("list", nruns)
dbObj_RPKMnorm_list       <- vector("list", nruns)
dbObj_dm6norm_list        <- vector("list", nruns)
dbObj_RLEbackground_list  <- vector("list", nruns)
dbObj_RLERiP_list         <- vector("list", nruns)
dbObj_greenlist_list      <- vector("list", nruns)

# --- Build DiffBind objects --------------------------------------------------

for (nrun in seq_len(nruns)) {
  dbObj_first_list[[nrun]] <- dba(sampleSheet = sample_information_randn[[nrun]])
}

# Count reads in consensus peak set (summit ± 200 bp; peaks in ≥ 2 samples)
dbObj_counts_summit_list <- lapply(dbObj_first_list, dba.count,
                                   summit = 200, minOverlap = 2)

# --- Normalization strategies ------------------------------------------------

# 1. Library-size normalization (DiffBind default / RPKM)
dbObj_RPKMnorm_list <- lapply(dbObj_counts_summit_list, dba.normalize)

# 2. Spike-in normalization using dm6 scaling factors
for (nrun in seq_len(nruns)) {
  dbObj_dm6norm_list[[nrun]] <- dba.normalize(
    dbObj_counts_summit_list[[nrun]],
    normalize = sample_information_randn[[nrun]]$dm6_reads_norm_scale
  )
}

# 3. RLE normalization using background reads
for (nrun in seq_len(nruns)) {
  dbObj_RLEbackground_list[[nrun]] <- dba.normalize(
    dbObj_counts_summit_list[[nrun]],
    normalize = "RLE", library = "background", background = TRUE
  )
}

# 4. RLE normalization using reads-in-peaks (RiP)
for (nrun in seq_len(nruns)) {
  dbObj_RLERiP_list[[nrun]] <- dba.normalize(
    dbObj_counts_summit_list[[nrun]],
    normalize = "RLE", library = "RiP"
  )
}

# 5. Greenlist normalization
for (nrun in seq_len(nruns)) {
  dbObj_greenlist_list[[nrun]] <- dba.normalize(
    dbObj_counts_summit_list[[nrun]],
    normalize = sample_information_randn[[nrun]]$greenlist_normalizer
  )
}

# --- Contrast and analyze ----------------------------------------------------

dbObj_dm6_contrast_list        <- lapply(dbObj_dm6norm_list,       dba.contrast)
dbObj_RPKMnorm_contrast_list   <- lapply(dbObj_RPKMnorm_list,      dba.contrast)
dbObj_RLEbackground_contrast_list <- lapply(dbObj_RLEbackground_list, dba.contrast)
dbObj_RLERiP_contrast_list     <- lapply(dbObj_RLERiP_list,        dba.contrast)
dbObj_greenlist_contrast_list  <- lapply(dbObj_greenlist_list,      dba.contrast)

dbObj_dm6_analyse_list         <- lapply(dbObj_dm6_contrast_list,         dba.analyze)
dbObj_RPKM_analyse_list        <- lapply(dbObj_RPKMnorm_contrast_list,    dba.analyze)
dbObj_RLEback_analyse_list     <- lapply(dbObj_RLEbackground_contrast_list, dba.analyze)
dbObj_RLERiP_analyse_list      <- lapply(dbObj_RLERiP_contrast_list,      dba.analyze)
dbObj_greenlist_analyse_list   <- lapply(dbObj_greenlist_contrast_list,    dba.analyze)

# --- Save workspace ----------------------------------------------------------

save.image(
  file.path(working_dir,
            "R_analysis/04.Negative_exp/data/diffbind_analyse_workspace.RData")
)

# --- Extract differential region counts --------------------------------------

FDR_cutoff <- 0.05

diff_nums_rec <- data.frame(matrix(
  nrow = nruns,
  ncol = 0
))

diffbind_results <- vector("list", nruns)

for (i in seq_len(nruns)) {

  # Sample combination label
  diff_nums_rec$rand_comb[i] <- paste(
    paste(sample_information_randn[[i]]$Condition[sample_information_randn[[i]]$Factor == "B"], collapse = ","),
    paste(sample_information_randn[[i]]$Condition[sample_information_randn[[i]]$Factor == "A"], collapse = ","),
    sep = "vs"
  ) |> gsub("million", "mill", x = _)

  # Total consensus peaks
  diff_nums_rec$num_of_peaks[i] <- dbObj_counts_summit_list[[i]] |>
    dba.peakset(bRetrieve = TRUE) |>
    as.data.frame() |>
    nrow()

  # Differential regions per normalization method (FDR < 0.05)
  diffbind_results[[i]][["RPKM_diff_regions"]] <-
    dba.report(dbObj_RPKM_analyse_list[[i]],    th = FDR_cutoff) |> as.data.frame()
  diffbind_results[[i]][["RLEback_diff_regions"]] <-
    dba.report(dbObj_RLEback_analyse_list[[i]], th = FDR_cutoff) |> as.data.frame()
  diffbind_results[[i]][["dm6scaled_diff_regions"]] <-
    dba.report(dbObj_dm6_analyse_list[[i]],     th = FDR_cutoff) |> as.data.frame()
  diffbind_results[[i]][["RLERiP_diff_regions"]] <-
    dba.report(dbObj_RLERiP_analyse_list[[i]],  th = FDR_cutoff) |> as.data.frame()
  diffbind_results[[i]][["greenlist_diff_regions"]] <-
    dba.report(dbObj_greenlist_analyse_list[[i]], th = FDR_cutoff) |> as.data.frame()

  # Record counts
  diff_nums_rec$RPKM_diff_regions[i]      <- nrow(diffbind_results[[i]][["RPKM_diff_regions"]])
  diff_nums_rec$RLEback_diff_regions[i]   <- nrow(diffbind_results[[i]][["RLEback_diff_regions"]])
  diff_nums_rec$dm6scaled_diff_regions[i] <- nrow(diffbind_results[[i]][["dm6scaled_diff_regions"]])
  diff_nums_rec$RLERiP_diff_regions[i]    <- nrow(diffbind_results[[i]][["RLERiP_diff_regions"]])
  diff_nums_rec$greenlist_diff_regions[i] <- nrow(diffbind_results[[i]][["greenlist_diff_regions"]])
}

# --- Save results ------------------------------------------------------------

saveRDS(
  diff_nums_rec,
  file.path(working_dir,
            "R_analysis/04.Negative_exp/data/diff_nums_rec.rds")
)
