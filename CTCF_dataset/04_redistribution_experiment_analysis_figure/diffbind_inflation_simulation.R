# =============================================================================
# DiffBind Normalization Comparison: Artificial Signal Inflation Simulation
# =============================================================================
# Description: Simulates differential binding by artificially inflating read
#              counts at randomly selected peaks in groups A and B. Evaluates
#              the sensitivity and specificity of five normalization methods
#              across randomized sample combinations.
# Input:       Randomized sample information list (sample_information_randn.rds)
# Output:      Per-run results, experiment summary table, and TP/FP/FN/TN
#              records saved per inflation factor
# =============================================================================

library(DiffBind)
library(tidyverse)

rm(list = ls())

# --- Paths -------------------------------------------------------------------

working_dir <- "/path/to/project"

# --- Load data ---------------------------------------------------------------

sample_information_randn <- readRDS(
  file.path(working_dir,
            "R_analysis/02.Sample_information/data/sample_information_randn.rds")
)

nruns <- length(sample_information_randn)

# --- Simulation parameters ---------------------------------------------------

frac_inflate_A <- 0.05  # Fraction of peaks inflated in group A
frac_inflate_B <- 0.05  # Fraction of peaks inflated in group B

# --- Simulation function -----------------------------------------------------

analysing_samples <- function(nrun_choose, inflation_factor) {

  sample_information_chosen <- sample_information_randn[[nrun_choose]]
  dbObj_first               <- dba(sampleSheet = sample_information_chosen)

  # Count reads in consensus peak set (summit ± 200 bp; peaks in ≥ 2 samples)
  dbObj_counts_summit_noninflate <- dba.count(dbObj_first, summit = 200, minOverlap = 2)
  dbObj_counts_summit_inflate    <- dbObj_counts_summit_noninflate

  # --- Inflate randomly selected peaks ---------------------------------------

  num_peaks <- dbObj_counts_summit_noninflate |>
    dba.peakset(bRetrieve = TRUE) |>
    as.data.frame() |>
    nrow()

  num_peaks_inflate_A <- floor(num_peaks * frac_inflate_A)
  num_peaks_inflate_B <- floor(num_peaks * frac_inflate_B)

  peak_idx <- sample.int(num_peaks, num_peaks_inflate_A + num_peaks_inflate_B)
  n_peaks_inflate_A <- peak_idx[seq_len(num_peaks_inflate_A)]
  n_peaks_inflate_B <- peak_idx[(num_peaks_inflate_A + 1):length(peak_idx)]

  score_cols <- c("Score", "RPKM", "Reads", "cRPKM", "cReads")

  for (n_peak_A in n_peaks_inflate_A) {
    for (s in 1:3) {
      dbObj_counts_summit_inflate[["peaks"]][[s]][n_peak_A, score_cols] <-
        dbObj_counts_summit_noninflate[["peaks"]][[s]][n_peak_A, score_cols] * inflation_factor
    }
  }

  for (n_peak_B in n_peaks_inflate_B) {
    for (s in 4:6) {
      dbObj_counts_summit_inflate[["peaks"]][[s]][n_peak_B, score_cols] <-
        dbObj_counts_summit_noninflate[["peaks"]][[s]][n_peak_B, score_cols] * inflation_factor
    }
  }

  dbObj_counts_summit_inflate <- dba(dbObj_counts_summit_inflate)

  # --- Build per-sample count records ----------------------------------------

  build_counts_rec <- function(dbObj) {
    sample_labels <- c("A1", "A2", "A3", "B1", "B2", "B3")
    sample_list   <- lapply(seq_along(sample_labels), function(s) {
      df <- as.data.frame(dbObj[["peaks"]][[s]])
      colnames(df)[4:8] <- paste(sample_labels[s], colnames(df)[4:8], sep = "_")
      df
    })
    rec <- Reduce(function(x, y) full_join(x, y, by = c("Chr", "Start", "End"),
                                           relationship = "one-to-one"),
                  sample_list)
    rec[, -grep("cRPKM|cReads", names(rec))]
  }

  inflate_rec    <- build_counts_rec(dbObj_counts_summit_inflate)
  noninflate_rec <- build_counts_rec(dbObj_counts_summit_noninflate)

  # Additional reads introduced by inflation
  inflate_counts    <- inflate_rec[,    grep("Reads", colnames(inflate_rec))]
  noninflate_counts <- noninflate_rec[, grep("Reads", colnames(noninflate_rec))]
  diff_counts       <- inflate_counts - noninflate_counts

  # --- Normalization and analysis --------------------------------------------

  # 1. Library-size normalization
  # Correct library sizes to account for artificially added reads
  dbObj_LS_tmp      <- dba.normalize(dbObj_counts_summit_inflate)
  new_lib_sizes     <- (dbObj_LS_tmp[["norm"]][["DESeq2"]][["lib.sizes"]] +
                          colSums(diff_counts)) |> unname()
  dbObj_LS_inflate    <- dba.normalize(dbObj_counts_summit_inflate,
                                       normalize = new_lib_sizes / mean(new_lib_sizes)) |>
    dba.contrast() |> dba.analyze() |> dba.report(th = 1) |> as.data.frame()
  dbObj_LS_noninflate <- dba.normalize(dbObj_counts_summit_noninflate) |>
    dba.contrast() |> dba.analyze() |> dba.report(th = 1) |> as.data.frame()

  # 2. RLE normalization using reads-in-peaks (RiP)
  dbObj_RLE_RiP_inflate <- dba.normalize(dbObj_counts_summit_inflate,
                                         normalize = "RLE", library = "RiP") |>
    dba.contrast() |> dba.analyze() |> dba.report(th = 1) |> as.data.frame()
  dbObj_RLE_RiP_noninflate <- dba.normalize(dbObj_counts_summit_noninflate,
                                            normalize = "RLE", library = "RiP") |>
    dba.contrast() |> dba.analyze() |> dba.report(th = 1) |> as.data.frame()

  # 3. Spike-in normalization using dm6 scaling factors
  dbObj_spikein_inflate <- dba.normalize(dbObj_counts_summit_inflate,
                                         normalize = sample_information_chosen$dm6_reads_norm_scale) |>
    dba.contrast() |> dba.analyze() |> dba.report(th = 1) |> as.data.frame()
  dbObj_spikein_noninflate <- dba.normalize(dbObj_counts_summit_noninflate,
                                            normalize = sample_information_chosen$dm6_reads_norm_scale) |>
    dba.contrast() |> dba.analyze() |> dba.report(th = 1) |> as.data.frame()

  # 4. Greenlist normalization
  dbObj_greenlist_inflate <- dba.normalize(dbObj_counts_summit_inflate,
                                           normalize = sample_information_chosen$greenlist_normalizer) |>
    dba.contrast() |> dba.analyze() |> dba.report(th = 1) |> as.data.frame()
  dbObj_greenlist_noninflate <- dba.normalize(dbObj_counts_summit_noninflate,
                                              normalize = sample_information_chosen$greenlist_normalizer) |>
    dba.contrast() |> dba.analyze() |> dba.report(th = 1) |> as.data.frame()

  # 5. RLE normalization using background reads
  # Add inflated peak reads to background bins before renormalizing
  dbObj_background_inflate <- dba.normalize(dbObj_counts_summit_inflate,
                                            normalize = "RLE", library = "background",
                                            background = TRUE)

  background_bins   <- as.data.frame(dbObj_background_inflate[["norm"]][["background"]][["binned"]]@rowRanges)
  background_counts <- as.data.frame(
    dbObj_background_inflate[["norm"]][["background"]][["binned"]]@assays@data@listData[["counts"]]
  )
  background_gr <- GRanges(data.frame(background_bins, background_counts))

  inflate_reads    <- inflate_rec[,    grep("Reads", colnames(inflate_rec))]
  noninflate_reads <- noninflate_rec[, grep("Reads", colnames(noninflate_rec))]

  reads_diff <- data.frame(inflate_rec[, c("Chr", "Start", "End")],
                           inflate_reads - noninflate_reads)
  reads_diff    <- reads_diff[rowSums(reads_diff[, 4:9]) > 0, ] |> GRanges()
  reads_diff_df <- as.data.frame(reads_diff)

  overlap_rec      <- as.data.frame(findOverlaps(reads_diff, background_gr))
  reads_diff_mat   <- as.data.frame(reads_diff)[, 6:11]

  background_counts_inflate <- background_counts
  background_counts_inflate[overlap_rec$subjectHits, ] <-
    background_counts[overlap_rec$subjectHits, ] +
    reads_diff_mat[overlap_rec$queryHits, ]

  dbObj_background_inflate[["norm"]][["background"]][["binned"]]@assays@data@listData[["counts"]] <-
    as.matrix(background_counts_inflate)

  dbObj_RLEbackground_inflate <- dba.normalize(dbObj_background_inflate,
                                               normalize = "RLE", library = "background",
                                               background = TRUE) |>
    dba.contrast() |> dba.analyze() |> dba.report(th = 1) |> as.data.frame()
  dbObj_RLEbackground_noninflate <- dba.normalize(dbObj_counts_summit_noninflate,
                                                  normalize = "RLE", library = "background",
                                                  background = TRUE) |>
    dba.contrast() |> dba.analyze() |> dba.report(th = 1) |> as.data.frame()

  # --- Return results --------------------------------------------------------

  results <- list(
    dbObj_first                    = dbObj_first,
    reads_diff                     = reads_diff_df,
    A_peaks_inflated               = n_peaks_inflate_A,
    B_peaks_inflated               = n_peaks_inflate_B,
    noninflate_rec                 = noninflate_rec,
    inflate_rec                    = inflate_rec,
    dbObj_counts_summit_noninflate = dbObj_counts_summit_noninflate,
    dbObj_counts_summit_inflate    = dbObj_counts_summit_inflate,
    LS_noninflate_df               = dbObj_LS_noninflate,
    LS_inflate_df                  = dbObj_LS_inflate,
    RLE_RiP_noninflate_df          = dbObj_RLE_RiP_noninflate,
    RLE_RiP_inflate_df             = dbObj_RLE_RiP_inflate,
    spikein_noninflate_df          = dbObj_spikein_noninflate,
    spikein_inflate_df             = dbObj_spikein_inflate,
    greenlist_noninflate_df        = dbObj_greenlist_noninflate,
    greenlist_inflate_df           = dbObj_greenlist_inflate,
    RLE_BG_noninflate_df           = dbObj_RLEbackground_noninflate,
    RLE_BG_inflate_df              = dbObj_RLEbackground_inflate
  )

  return(results)
}

# --- Run simulation across inflation factors ---------------------------------

pval_threshold  <- 0.05
rand_runs       <- paste0("rand", seq_len(nruns), ".rds")
norm_methods    <- c("LS", "RLE_RiP", "spikein", "greenlist", "RLE_BG")

for (inflation_factor_choose in c(2)) {

  results_dir <- file.path(working_dir,
                           "R_analysis/05.Artificially_inflate/results",
                           paste0("inflation_factor_", inflation_factor_choose))

  print(paste("Inflation factor:", inflation_factor_choose))

  inflated_DESeq_record  <- list()
  experiment_summary_nrand <- data.frame(rand_runs = paste0("rand_", seq_len(nruns)))

  for (nrun in seq_len(nruns)) {
    print(paste("Run:", nrun))

    randn    <- analysing_samples(nrun, inflation_factor = inflation_factor_choose)
    rand_int <- paste0("rand", nrun)

    saveRDS(randn, file.path(results_dir, "individual_runs", rand_runs[nrun]))

    # --- Experiment metadata -------------------------------------------------

    experiment_summary_nrand[nrun, "groupA_samples"] <-
      randn[["dbObj_first"]][["masks"]][["A"]][randn[["dbObj_first"]][["masks"]][["A"]]] |>
      names() |> paste(collapse = ", ")

    experiment_summary_nrand[nrun, "groupB_samples"] <-
      randn[["dbObj_first"]][["masks"]][["B"]][randn[["dbObj_first"]][["masks"]][["B"]]] |>
      names() |> paste(collapse = ", ")

    experiment_summary_nrand[nrun, "inflation_factor"]      <- inflation_factor_choose
    experiment_summary_nrand[nrun, "npeaks_consensus"]      <- nrow(randn[["noninflate_rec"]])
    experiment_summary_nrand[nrun, "npeaks_inflated_groupA"] <- length(randn[["A_peaks_inflated"]])
    experiment_summary_nrand[nrun, "npeaks_inflated_groupB"] <- length(randn[["B_peaks_inflated"]])

    # --- Differential peak counts and TP/FP/FN/TN per method -----------------

    for (norm_method in norm_methods) {

      ni_df <- randn[[paste0(norm_method, "_noninflate_df")]]
      i_df  <- randn[[paste0(norm_method, "_inflate_df")]]

      # Non-inflated counts
      experiment_summary_nrand[nrun, paste0(norm_method, "_NI_diffpeaks")] <-
        sum(ni_df$FDR < pval_threshold)
      experiment_summary_nrand[nrun, paste0(norm_method, "_NI_uppeaks")] <-
        sum(ni_df$FDR < pval_threshold & ni_df$Fold > 0)
      experiment_summary_nrand[nrun, paste0(norm_method, "_NI_downpeaks")] <-
        sum(ni_df$FDR < pval_threshold & ni_df$Fold < 0)

      # Inflated counts
      experiment_summary_nrand[nrun, paste0(norm_method, "_I_diffpeaks")] <-
        sum(i_df$FDR < pval_threshold)
      experiment_summary_nrand[nrun, paste0(norm_method, "_I_uppeaks")] <-
        sum(i_df$FDR < pval_threshold & i_df$Fold > 0)
      experiment_summary_nrand[nrun, paste0(norm_method, "_I_downpeaks")] <-
        sum(i_df$FDR < pval_threshold & i_df$Fold < 0)

      # TP/FP/FN/TN classification
      DESeq_inflate <- left_join(i_df, randn[["reads_diff"]])
      DESeq_inflate$A_inflate <- rowSums(DESeq_inflate[, c("A1_Reads", "A2_Reads", "A3_Reads")]) > 0
      DESeq_inflate$B_inflate <- rowSums(DESeq_inflate[, c("B1_Reads", "B2_Reads", "B3_Reads")]) > 0
      DESeq_inflate$inflate   <- !is.na(DESeq_inflate$A_inflate | DESeq_inflate$B_inflate)

      DESeq_inflate$TP <- DESeq_inflate$FDR < 0.05 &  DESeq_inflate$inflate
      DESeq_inflate$FP <- DESeq_inflate$FDR < 0.05 & !DESeq_inflate$inflate
      DESeq_inflate$FN <- DESeq_inflate$FDR > 0.05 &  DESeq_inflate$inflate
      DESeq_inflate$TN <- DESeq_inflate$FDR > 0.05 & !DESeq_inflate$inflate

      inflated_DESeq_record[[rand_int]][[paste0(norm_method, "_DESeq_inflate")]] <- DESeq_inflate

      experiment_summary_nrand[nrun, paste0(norm_method, "_I_TP")] <- sum(na.omit(DESeq_inflate$TP))
      experiment_summary_nrand[nrun, paste0(norm_method, "_I_FP")] <- sum(na.omit(DESeq_inflate$FP))
      experiment_summary_nrand[nrun, paste0(norm_method, "_I_FN")] <- sum(na.omit(DESeq_inflate$FN))
      experiment_summary_nrand[nrun, paste0(norm_method, "_I_TN")] <- sum(na.omit(DESeq_inflate$TN))
    }
  }

  # --- Save results ----------------------------------------------------------

  saveRDS(experiment_summary_nrand,
          file.path(results_dir, "experiment_summary.rds"))
  saveRDS(inflated_DESeq_record,
          file.path(results_dir, "inflated_DESeq_record.rds"))
}
