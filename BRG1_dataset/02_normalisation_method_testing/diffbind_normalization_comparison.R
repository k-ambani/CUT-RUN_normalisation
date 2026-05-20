# =============================================================================
# DiffBind Normalization Comparison: 2v2 Analysis
# =============================================================================
# Description: Compares four normalization strategies for CUT&RUN data using
#              DiffBind: library-size (RPKM), spike-in (sacCer3), RLE with
#              background reads, and RLE with reads-in-peaks (RiP).
# Input:       BAM files (unfiltered), peak files, and spike-in counts
# Output:      DiffBind analysis objects saved as .RData workspace
# =============================================================================

# --- Environment setup -------------------------------------------------------


library(DiffBind)
library(tidyverse)

rm(list = ls())

# --- Paths -------------------------------------------------------------------

working_dir <- paste0(
  "/researchers/krutika.ambani/Goel_lab_members/Krutika_Ambani/",
  "CnR_Normalisation_paper/goel_lab_datasets/20241230_BRG1_with_replicates"
)

# --- Load data ---------------------------------------------------------------

sample_information <- readRDS(
  file.path(working_dir, "R_analysis/1.Sample_information/data/sample_information.rds")
)

# Convert replicate labels (e.g. "R1") to integers
sample_information$Replicate <- sample_information$Replicate |>
  gsub("R", "", x = _) |>
  as.numeric()

summary_bam_counts <- readRDS(
  file.path(working_dir, "R_analysis/1.Sample_information/data/spike_in_counts.rds")
)

# Compute per-sample sacCer3 spike-in scaling factors
summary_bam_counts$sacCer3_reads_norm_scale <-
  summary_bam_counts$nohg38_sacCer3 / mean(summary_bam_counts$nohg38_sacCer3)

# Subset to first four samples (2v2 comparison)
sample_information  <- sample_information[1:4, ]
summary_bam_counts  <- summary_bam_counts[1:4, ]

# --- Build DiffBind object ---------------------------------------------------

dbObj_first <- dba(sampleSheet = sample_information)

# Count reads in consensus peak set (summit ± 200 bp; peaks called in ≥ 2 samples)
dbObj_counts_summit <- dba.count(dbObj_first, summit = 200, minOverlap = 2)

# --- Normalization strategies ------------------------------------------------

# 1. Library-size normalization (DiffBind default / RPKM)
dbObj_RPKMnorm <- dbObj_counts_summit |>
  dba.normalize() |>
  dba.contrast(contrast = c("Treatment", "Abema", "DMSO")) |>
  dba.analyze()

# 2. Spike-in normalization using sacCer3 scaling factors
dbObj_sacCer3norm <- dbObj_counts_summit |>
  dba.normalize(normalize = summary_bam_counts$sacCer3_reads_norm_scale) |>
  dba.contrast(contrast = c("Treatment", "Abema", "DMSO")) |>
  dba.analyze()

# 3. RLE normalization using background reads
dbObj_RLEbackground <- dbObj_counts_summit |>
  dba.normalize(normalize = "RLE", library = "background", background = TRUE) |>
  dba.contrast(contrast = c("Treatment", "Abema", "DMSO")) |>
  dba.analyze()

# 4. RLE normalization using reads-in-peaks (RiP)
dbObj_RLERiP <- dbObj_counts_summit |>
  dba.normalize(normalize = "RLE", library = "RiP") |>
  dba.contrast(contrast = c("Treatment", "Abema", "DMSO")) |>
  dba.analyze()

# --- Save workspace ----------------------------------------------------------

save.image(
  file.path(
    working_dir,
    "R_analysis/3.Diffbind_comparison_of_interest/data/diffbind_analyse_workspace_2v2.RData"
  )
)
