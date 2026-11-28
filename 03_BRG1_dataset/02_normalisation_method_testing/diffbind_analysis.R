
library(DiffBind)
library(tidyverse)
library(utils)

rm(list=ls())


############################################################################################
# RUNNING 2V2 combinations
############################################################################################


#-------------------------------------------------------------------------------------------------------#
# Creating sample information
#-------------------------------------------------------------------------------------------------------#

working_dir <- "/path/to/project"
### -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# bam files directory
hg38_bamReads_dir <- file.path(
  paste(working_dir, "/results/bam_files/final_bams/hg38", sep = ""),
  list.files(
    path = paste(working_dir, "/results/bam_files/final_bams/hg38", sep = ""),
    pattern = "_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam$"
  )
)

# need to remove the IgG  negative controls
hg38_bamReads_dir %>% grep("IgG", ., value = TRUE)
index <- !hg38_bamReads_dir %>% grepl("IgG", .)
hg38_bamReads_dir <- hg38_bamReads_dir[index]

### -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# peaks
peaks_dir = "/path/to/project/results/macs2/with_input/macs2_q0.05/narrowPeak"
# Peaks directory
Peaks_dir_macs2_with_input <- file.path(
  paste(peaks_dir),
  list.files(
    path = peaks_dir,
    pattern = "Par.*\\.*narrowPeak"
  )
)

### -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# check peaks and bam samples are in the same order
peak_names <- Peaks_dir_macs2_with_input %>%
  gsub("_peaks.narrowPeak", "", .) %>%
  gsub("/path/to/project/results/macs2/with_input/macs2_q0.05/narrowPeak/", "", .)

bam_names <- hg38_bamReads_dir %>%
  gsub("/path/to/project/results/bam_files/final_bams/hg38/", "", .) %>%
  gsub("_chrmrm_blrm_duprm_unmappedrm_multimaprm.bam", "", .)

# same order, same name
table(peak_names == bam_names)

# SampleID
SampleID <- bam_names

SampleID <- SampleID %>% gsub("-native_S.*", "", .)

sample_information <- data.frame(
  SampleID = SampleID,
  bamReads = hg38_bamReads_dir,
  Peaks = Peaks_dir_macs2_with_input
)

# narrow because looking at transcription factors.
sample_information$PeakCaller <- "narrow"

sample_information <- sample_information %>% separate(SampleID, into = c("Sample_no","Tissue", "Condition", "Treatment", "Replicate", "Antibody"), sep = "-", remove = FALSE)
sample_information <- sample_information[grepl("-BRG1",sample_information$SampleID),]

sample_information$Factor <- paste(sample_information$Condition, sample_information$Treatment,sep="_")




#---------------------------------------#
# Spike-in norm factors
#---------------------------------------#

summary_bam_counts <- read_tsv(paste(working_dir,"/results/summary_hg38_sacCer3_bam_counts.txt",sep=""))

col_names <- names(summary_bam_counts) %>% strsplit(split=" ") %>% unlist()
summary_bam_counts <- separate(summary_bam_counts,1,into=col_names,sep=" ")

SampleID <- summary_bam_counts$basename %>% gsub("_onlyhg38_chrmrm_blrm_duprm_unmappedrm_multimaprm","",.)
SampleID <- SampleID  %>% sub("_.*", "", .)
SampleID <- SampleID %>% gsub("-native", "", .)

summary_bam_counts$SampleID <- SampleID


summary_bam_counts$sacCer3_intermediate_bam_counts <- as.numeric(summary_bam_counts$sacCer3_intermediate_bam_counts)

summary_bam_counts$sacCer3_reads_norm_scale <- summary_bam_counts$sacCer3_intermediate_bam_counts/mean(summary_bam_counts$sacCer3_intermediate_bam_counts)
summary_bam_counts <- summary_bam_counts[1:4,]

sample_information <- full_join(sample_information,summary_bam_counts)




#---------------------------------------#
# Greenlist norm factors
#---------------------------------------#
greenlist<- read.table("/path/to/project/R_analysis/5.Diffbind_comparisons_greenlist_normalisation_incl/greenlist_normfacs/glist_sizeFactors.tsv", sep = '\t', header = TRUE) %>% as.data.frame()
sample_df <- data.frame(sample = rownames(greenlist))
rownames(greenlist) = 1:4

sample_df <- separate(sample_df, col = sample, into = c("Sample_no","Tissue", "Condition", "Treatment", "Replicate", "Antibody"), sep = "-", remove = TRUE)
sample_df$SampleID <- paste(sample_df$Sample_no, sample_df$Tissue, sample_df$Condition, sample_df$Treatment, sample_df$Replicate,sample_df$Antibody, sep="-")
greenlist <- data.frame(sample_df,greenlist)
greenlist <- greenlist[,c("SampleID","sf","normalizer")]

sample_information_greenlist <- data.frame(greenlist, 
                                           greenlist_normaliser = 1/greenlist$normalizer)


sample_information <- full_join(sample_information,sample_information_greenlist)

#---------------------------------------#
# Creating dba objects
#---------------------------------------#

# Creating initial dba objects
dbObj_first <- dba(sampleSheet = sample_information)


dbObj_counts_summit <- dba.count(dbObj_first, summit=200,minOverlap = 2)
dbObj_counts_summit$config$cores <- 32

#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING library size NORMALISATION
#--------------------------------------#
# diffbind default params
  dbObj_RPKMnorm <- dba.normalize(dbObj_counts_summit) %>%
    dba.contrast(.,contrast = c("Treatment","Abema","DMSO")) %>%
    dba.analyze()


#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING sacCer3 spike-in NORMALISATION
#--------------------------------------#

dbObj_sacCer3norm <- dba.normalize(dbObj_counts_summit, normalize = sample_information$sacCer3_reads_norm_scale) %>%
  dba.contrast(.,contrast =c("Treatment","Abema","DMSO")) %>%
  dba.analyze()

#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING background reads and RLE
#--------------------------------------#

  dbObj_RLEbackground <-dba.normalize(dbObj_counts_summit, normalize = "RLE", library = "background", background=TRUE) %>%
    dba.contrast(.,contrast = c("Treatment","Abema","DMSO")) %>%
    dba.analyze()

#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING Reads in Peak and RLE
#--------------------------------------#

  dbObj_RLERiP <- dba.normalize(dbObj_counts_summit, normalize = "RLE", library = "RiP")  %>%
    dba.contrast(.,contrast = c("Treatment","Abema","DMSO")) %>%
    dba.analyze()

  
  #---------------------------------------#
  # ANALYSING BY CELLLINE-ANTIBODY USING Reads in Peak and RLE
  #--------------------------------------#
  
  
  dbObj_greenlistnorm <- dba.normalize(dbObj_counts_summit, normalize = sample_information$greenlist_normaliser) %>%
    dba.contrast(.,contrast =c("Treatment","Abema","DMSO")) %>%
    dba.analyze()
  
  
  
save.image(file="/path/to/project/R_analysis_reviewerresponse/1.Diffbind_analysis/data/diffbind_analyse_workspace_2v2.RData")

