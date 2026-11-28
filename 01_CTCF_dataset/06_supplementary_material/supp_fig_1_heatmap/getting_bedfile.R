library(rtracklayer)
library(tidyverse)
library(DiffBind)
library(GenomicFeatures)




# Loading data

sample_information_all <- readRDS("/path/to/project/R_analysis_final/02.Sample_information/data/sample_information_all.rds")



dbObj_first <- dba(sampleSheet = sample_information_all)
dbObj_counts_summit <- dba.count(dbObj_first, summit=200,minOverlap = 4)

# Retrieve the consensus peakset as a GRanges object
consensus_peaks <- dba.peakset(dbObj_counts_summit, bRetrieve = TRUE)

# Export directly as BED
rtracklayer::export.bed(consensus_peaks, "/path/to/project/R_analysis_final/03.Deeptools/bedfile/consensus_peakset.bed")







