# =============================================================================
# DiffBind Normalization Comparison: 2v2 Analysis
# =============================================================================
# Description: Compares five normalization strategies for CUT&RUN data using
#              DiffBind: library-size (RPKM), spike-in (sacCer3), RLE with
#              background reads, RLE with reads-in-peaks (RiP), and greenlist.
# =============================================================================

library(DiffBind)
library(tidyverse)
library(utils)

rm(list = ls())

# --- Paths -------------------------------------------------------------------

working_dir <- "/path/to/project"

peaks_dir <- file.path(working_dir, "results/macs2/with_input/macs2_q0.05/narrowPeak")

# --- BAM files ---------------------------------------------------------------

hg38_bamReads_dir <- file.path(
  working_dir, "results/bam_files/final_bams/hg38",
  list.files(
    path    = file.path(working_dir, "results/bam_files/final_bams/hg38"),
    pattern = "_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam$"
  )
)

# Remove IgG negative controls
index             <- !grepl("IgG", hg38_bamReads_dir)
hg38_bamReads_dir <- hg38_bamReads_dir[index]

# --- Peak files --------------------------------------------------------------

Peaks_dir_macs2_with_input <- file.path(
  peaks_dir,
  list.files(path = peaks_dir, pattern = "*narrowPeak")
)

# --- Verify BAM and peak files are in the same order -------------------------

peak_names <- Peaks_dir_macs2_with_input |>
  gsub("_peaks.narrowPeak", "", x = _) |>
  gsub(paste0(peaks_dir, "/"), "", x = _)

bam_names <- hg38_bamReads_dir |>
  gsub(file.path(working_dir, "results/bam_files/final_bams/hg38/"), "", x = _) |>
  gsub("_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam", "", x = _)

table(peak_names == bam_names)

# --- Build sample sheet ------------------------------------------------------

SampleID <- gsub("-native_S.*", "", bam_names)

sample_information <- data.frame(
  SampleID = SampleID,
  bamReads = hg38_bamReads_dir,
  Peaks    = Peaks_dir_macs2_with_input,
  PeakCaller = "narrow"
) |>
  separate(SampleID,
           into    = c("Sample_no", "Tissue", "Condition", "Treatment", "Replicate", "Antibody"),
           sep     = "-",
           remove  = FALSE) |>
  filter(grepl("-BRG1", SampleID)) |>
  mutate(Factor = paste(Condition, Treatment, sep = "_"))

# --- Spike-in normalization factors ------------------------------------------

summary_bam_counts <- readRDS(
  file.path(working_dir, "R_analysis/1.Sample_information/data/spike_in_counts.rds")
)

summary_bam_counts$sacCer3_reads_norm_scale <-
  summary_bam_counts$nohg38_sacCer3 / mean(summary_bam_counts$nohg38_sacCer3)

# --- Greenlist normalization factors -----------------------------------------

greenlist <- read.table(
  file.path(working_dir,
            "R_analysis/5.Diffbind_comparisons_greenlist_normalisation_incl",
            "greenlist_normfacs/glist_sizeFactors.tsv"),
  sep = "\t", header = TRUE
) |> as.data.frame()

sample_df       <- data.frame(sample = rownames(greenlist))
rownames(greenlist) <- seq_len(nrow(greenlist))

sample_i <- separate(sample_df, col = sample,
                     into   = c("Sample_no", "Tissue", "Condition", "Treatment", "Replicate", "Antibody"),
                     sep    = "-",
                     remove = TRUE) |>
  mutate(SampleID = paste(Sample_no, Tissue, Condition, Treatment, Replicate, Antibody, sep = "-"))

greenlist <- data.frame(sample_i, greenlist)[, c("SampleID", "sf", "normalizer")]

sample_information_greenlist <- data.frame(
  greenlist,
  greenlist_normaliser = 1 / greenlist$normalizer
)

sample_information <- full_join(sample_information, sample_information_greenlist)

# --- Subset to samples of interest (2v2) -------------------------------------

sample_information <- sample_information[1:4, ]
summary_bam_counts <- summary_bam_counts[1:4, ]

# --- Build DiffBind object ---------------------------------------------------

dbObj_first         <- dba(sampleSheet = sample_information)
dbObj_counts_summit <- dba.count(dbObj_first, summit = 200, minOverlap = 2)
dbObj_counts_summit$config$cores <- 32

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

# 5. Greenlist normalization
dbObj_greenlistnorm <- dbObj_counts_summit |>
  dba.normalize(normalize = sample_information$greenlist_normaliser) |>
  dba.contrast(contrast = c("Treatment", "Abema", "DMSO")) |>
  dba.analyze()

# --- Save workspace ----------------------------------------------------------

save.image(
  file.path(working_dir,
            "R_analysis/5.Diffbind_comparisons_greenlist_normalisation_incl",
            "data/diffbind_analyse_workspace_2v2.RData")
)
