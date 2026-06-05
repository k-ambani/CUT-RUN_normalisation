# =============================================================================
# Sample Information and DiffBind Normalization Comparison: 2v2 Combinations
# =============================================================================
# Description: Builds the sample information table (BAM paths, peak paths,
#              greenlist and spike-in scaling factors), then applies six
#              normalization strategies to all pairwise 2v2 sample combinations.
#              Normalization is performed once on the full dataset; contrasts
#              are specified per combination.
# Output:      sample_information.rds, DiffBind analysis workspace
# =============================================================================

library(DiffBind)
library(tidyverse)

rm(list = ls())

# --- Paths -------------------------------------------------------------------

working_dir <- "/path/to/project"

bam_dir   <- file.path(working_dir, "results/bam_files/final_bams/hg38")
peaks_dir <- file.path(working_dir, "results/macs2/withinput/macs2_q0.05/narrowPeak")

# =============================================================================
# PART 1: Build sample information
# =============================================================================

# --- BAM files ---------------------------------------------------------------

hg38_bamReads_dir <- file.path(
  bam_dir,
  list.files(path = bam_dir, pattern = "_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam$")
)

# Remove IgG negative controls
hg38_bamReads_dir <- hg38_bamReads_dir[!grepl("IgG", hg38_bamReads_dir)]

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
  gsub(paste0(bam_dir, "/"), "", x = _) |>
  gsub("_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam", "", x = _)

table(peak_names == bam_names)

# --- Build sample sheet ------------------------------------------------------

# Standardise sample IDs
clean_sampleID <- function(x) {
  x |>
    gsub("_S.*", "", x = _) |>
    gsub("MCF7-CSS-R1-ER",   "MCF7-CSS-noE2-R1-ER",       x = _) |>
    gsub("MCF7-CSS-R2-ER",   "MCF7-CSS-noE2-R2-ER",       x = _) |>
    gsub("MCF7-Complete",    "MCF7-Complete-sgControl",    x = _) |>
    gsub("sgER-Complete",    "MCF7-Complete-sgER",         x = _)
}

SampleID <- clean_sampleID(bam_names)

sample_i <- data.frame(
  SampleID = SampleID,
  bamReads = hg38_bamReads_dir,
  Peaks    = Peaks_dir_macs2_with_input,
  PeakCaller = "narrow"
) |>
  separate(SampleID,
           into   = c("Tissue", "Condition", "Treatment", "Replicate", "Antibody"),
           sep    = "-",
           remove = FALSE) |>
  mutate(Condition = factor(Condition, levels = c("CSS", "Complete"))) |>
  arrange(Condition)

# --- Greenlist normalization factors -----------------------------------------

greenlist <- read.table(
  file.path(working_dir,
            "R_analysis/1.1Greenlist_norm/greenlist_norm/glist_sizeFactors.tsv"),
  sep = "\t", header = TRUE
)

greenlist_SampleID <- rownames(greenlist) |>
  gsub("_onlyhg38_chrmrm_blrm_duprm_unmappedrm_multimaprm", "", x = _) |>
  sub("_.*", "", x = _) |>
  clean_sampleID()

sample_information_greenlist <- data.frame(
  SampleID             = greenlist_SampleID,
  greenlist_normaliser = 1 / greenlist$normalizer
)

sample_i <- full_join(sample_i, sample_information_greenlist)

# --- Spike-in normalization factors ------------------------------------------

summary_bam_counts <- read_tsv(
  file.path(working_dir, "results/summary_bam_counts.txt")
)

col_names          <- unlist(strsplit(names(summary_bam_counts), split = " "))
summary_bam_counts <- separate(summary_bam_counts, 1, into = col_names, sep = " ")

summary_bam_counts$SampleID <- summary_bam_counts$basename |>
  gsub("_onlyhg38_chrmrm_blrm_duprm_unmappedrm_multimaprm", "", x = _) |>
  sub("_.*", "", x = _) |>
  clean_sampleID()

summary_bam_counts <- summary_bam_counts |>
  mutate(
    dm6_intermediate_bam_counts    = as.numeric(dm6_intermediate_bam_counts),
    sacCer3_intermediate_bam_counts = as.numeric(sacCer3_intermediate_bam_counts),
    dm6_reads_norm_scale            = dm6_intermediate_bam_counts    / mean(dm6_intermediate_bam_counts),
    sacCer3_reads_norm_scale        = sacCer3_intermediate_bam_counts / mean(sacCer3_intermediate_bam_counts)
  )

sample_information <- full_join(sample_i, summary_bam_counts)

# --- Save sample information -------------------------------------------------

saveRDS(
  sample_information,
  file.path(working_dir, "R_analysis/1.Sample_information/data/sample_information.rds")
)

# =============================================================================
# PART 2: DiffBind 2v2 normalization comparison
# =============================================================================

sample_information$Replicate <- as.numeric(gsub("R", "", sample_information$Replicate))
sample_information$Factor    <- sample_information$Treatment

# --- Define 2v2 sample combinations ------------------------------------------

nruns <- 6

rand_combos_2 <- data.frame(matrix(nrow = nruns, ncol = 4))
names(rand_combos_2) <- c("A1", "A2", "B1", "B2")

rand_combos_2[1, ] <- c(3, 4, 1, 2)
rand_combos_2[2, ] <- c(5, 6, 1, 2)
rand_combos_2[3, ] <- c(7, 8, 1, 2)
rand_combos_2[4, ] <- c(3, 4, 5, 6)
rand_combos_2[5, ] <- c(7, 8, 3, 4)
rand_combos_2[6, ] <- c(7, 8, 5, 6)

sample_information_combs <- vector("list", nruns)
contrast_list             <- vector("list", nruns)

for (nrun in seq_len(nruns)) {
  sample_information_combs[[nrun]] <- sample_information[as.numeric(rand_combos_2[nrun, ]), ]
  sample_information_combs[[nrun]]$Factor <- sample_information_combs[[nrun]]$Treatment
  contrast_list[[nrun]] <- unique(sample_information_combs[[nrun]]$Treatment)
}

# --- Build DiffBind object ---------------------------------------------------

dbObj_first <- dba(sampleSheet = sample_information, minOverlap = 1)

# Count reads in consensus peak set (summit ± 200 bp; peaks in ≥ 1 sample)
dbObj_counts_summit <- dba.count(dbObj_first, summit = 200, minOverlap = 1)

# --- Normalization strategies (applied to full dataset) ----------------------

# 1. Library-size normalization (DiffBind default / RPKM)
dbObj_RPKMnorm      <- dba.normalize(dbObj_counts_summit)

# 2. Spike-in normalization using sacCer3 scaling factors
dbObj_sacCer3norm   <- dba.normalize(dbObj_counts_summit,
                                     normalize = sample_information$sacCer3_reads_norm_scale)

# 3. Spike-in normalization using dm6 scaling factors
dbObj_dm6norm       <- dba.normalize(dbObj_counts_summit,
                                     normalize = sample_information$dm6_reads_norm_scale)

# 4. RLE normalization using background reads
dbObj_RLEbackground <- dba.normalize(dbObj_counts_summit,
                                     normalize = "RLE", library = "background",
                                     background = TRUE)

# 5. RLE normalization using reads-in-peaks (RiP)
dbObj_RLERiP        <- dba.normalize(dbObj_counts_summit,
                                     normalize = "RLE", library = "RiP")

# 6. Greenlist normalization
dbObj_greenlist     <- dba.normalize(dbObj_counts_summit,
                                     normalize = sample_information$greenlist_normaliser)

# --- Contrast and analyze per combination ------------------------------------

dbObj_RPKMnorm_list      <- vector("list", nruns)
dbObj_sacCer3norm_list   <- vector("list", nruns)
dbObj_dm6norm_list       <- vector("list", nruns)
dbObj_RLEbackground_list <- vector("list", nruns)
dbObj_RLERiP_list        <- vector("list", nruns)
dbObj_greenlist_list     <- vector("list", nruns)

for (nrun in seq_len(nruns)) {
  contrast <- c("Factor", contrast_list[[nrun]][2], contrast_list[[nrun]][1])

  dbObj_RPKMnorm_list[[nrun]]      <- dbObj_RPKMnorm      |> dba.contrast(contrast = contrast) |> dba.analyze()
  dbObj_sacCer3norm_list[[nrun]]   <- dbObj_sacCer3norm   |> dba.contrast(contrast = contrast) |> dba.analyze()
  dbObj_dm6norm_list[[nrun]]       <- dbObj_dm6norm       |> dba.contrast(contrast = contrast) |> dba.analyze()
  dbObj_RLEbackground_list[[nrun]] <- dbObj_RLEbackground |> dba.contrast(contrast = contrast) |> dba.analyze()
  dbObj_RLERiP_list[[nrun]]        <- dbObj_RLERiP        |> dba.contrast(contrast = contrast) |> dba.analyze()
  dbObj_greenlist_list[[nrun]]     <- dbObj_greenlist      |> dba.contrast(contrast = contrast) |> dba.analyze()
}

# --- Save workspace ----------------------------------------------------------

save.image(
  file.path(working_dir, "R_analysis/2.Diffbind/data/diffbind_analyse_workspace_2v2.RData")
)
