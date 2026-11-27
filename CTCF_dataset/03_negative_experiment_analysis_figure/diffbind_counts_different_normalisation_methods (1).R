library(DiffBind)
library(tidyverse)
library(BiocParallel)

# Check what SLURM gave you
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 4))
print(ncores)
register(MulticoreParam(workers = min(ncores, 8)))  # Cap at 8 for safety
set.seed(10)
rm(list=ls())


#---------------------------------------#
# Loading data
#---------------------------------------#

#loading sample_information workspace (spike-in and greenlist scale factors included )
#load("/data/gpfs/projects/punim2745/03.CUTRUN_normalisation_KA/Henikoff_simulations/R_analysis_final/02.Sample_information/data/sample_information_workspace.RData")
sample_information_randn <- readRDS("/data/gpfs/projects/punim2745/03.CUTRUN_normalisation_KA/Henikoff_simulations/R_analysis_final/02.Sample_information/data/sample_information_randn.rds")
#---------------------------------------#
# Initialising
#---------------------------------------#

#output objects
nruns = length(sample_information_randn)
dbObj_first_list <- vector(mode = "list", length = nruns)
dbObj_counts_summit_list<- vector(mode="list", length = nruns)

dbObj_dm6norm_list<- vector(mode="list", length = nruns)
dbObj_RPKMnorm_list <- vector(mode = "list", length=nruns)
dbObj_RLEbackground_list <- vector(mode="list", length=nruns)
dbObj_RLERiP_list <- vector(mode="list", length=nruns)
dbObj_greenlist_list <- vector(mode="list", length=nruns)



#---------------------------------------#
# Starting main experiment
#---------------------------------------#

for (nrun in 1:nruns){

  # Creating initial dba objects
 dbObj_first_list[[nrun]] <- dba(sampleSheet = sample_information_randn[[nrun]])

}



#---------------------------------------#
# Creating dba counts object
#---------------------------------------#
dbObj_counts_summit_list <- lapply(dbObj_first_list, dba.count, summit=200,minOverlap = 2)



#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING library size NORMALISATION
#--------------------------------------#
# diffbind default params
dbObj_RPKMnorm_list <- lapply(dbObj_counts_summit_list,dba.normalize)

#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING spike-in NORMALISATION
#--------------------------------------#

for (nrun in 1:nruns){
  dbObj_dm6norm_list[[nrun]] <- dba.normalize(dbObj_counts_summit_list[[nrun]], normalize = sample_information_randn[[nrun]]$dm6_reads_norm_scale)

}


#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING background reads and RLE
#--------------------------------------#


for (nrun in 1:nruns){
  dbObj_counts_summit_list[[nrun]]$config$cores <- 32
  dbObj_RLEbackground_list[[nrun]] <- dba.normalize(dbObj_counts_summit_list[[nrun]], normalize = "RLE", library = "background", background=TRUE)
}

#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING Reads in Peak and RLE
#--------------------------------------#

for (nrun in 1:nruns){
  #dbObj_counts_summit_list[[nrun]]$config$cores <- 32   #https://support.bioconductor.org/p/p133921/
  dbObj_RLERiP_list[[nrun]] <- dba.normalize(dbObj_counts_summit_list[[nrun]], normalize = "RLE", library = "RiP")
}


#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING greenlist normalisation method
#--------------------------------------#

for (nrun in 1:nruns){
  dbObj_greenlist_list[[nrun]] <- dba.normalize(dbObj_counts_summit_list[[nrun]], normalize = sample_information_randn[[nrun]]$greenlist_normalizer)
}



#---------------------------------------#
# Specifying contrasts and analysing and saving workspace
#--------------------------------------#


dbObj_dm6_contrast_list <- lapply(dbObj_dm6norm_list,dba.contrast)
dbObj_RPKMnorm_contrast_list <- lapply(dbObj_RPKMnorm_list,dba.contrast)
dbObj_RLEbackground_contrast_list <- lapply(dbObj_RLEbackground_list,dba.contrast)
dbObj_RLERiP_contrast_list <- lapply(dbObj_RLERiP_list,dba.contrast)
dbObj_greenlist_contrast_list <- lapply(dbObj_greenlist_list, dba.contrast)

dbObj_dm6_analyse_list <- lapply(dbObj_dm6_contrast_list,dba.analyze)
dbObj_RPKM_analyse_list <- lapply(dbObj_RPKMnorm_contrast_list,dba.analyze)
dbObj_RLEback_analyse_list <- lapply(dbObj_RLEbackground_contrast_list,dba.analyze)
dbObj_RLERiP_analyse_list <- lapply(dbObj_RLERiP_contrast_list,dba.analyze)
dbObj_greenlist_analyse_list <- lapply(dbObj_greenlist_contrast_list, dba.analyze)



save.image(file="/data/gpfs/projects/punim2745/03.CUTRUN_normalisation_KA/Henikoff_simulations/R_analysis_response_to_reviewers/05.Negative_experiment/data/diffbind_analyse_workspace.RData")


#---------------------------------------#
# Downstream analysis
#--------------------------------------#

### initialising results list
# Preallocating nested list
num_outer_elements <- nruns
num_inner_elements <- 0

# Preallocate the nested list
template_nested <- vector("list", length = num_outer_elements)

# Loop through the outer list and initialize inner lists
for (i in 1:num_outer_elements) {
  template_nested[[i]] <- vector("list", length = num_inner_elements)
}



FDR_cutoff = 0.05


diffbind_results_all <- template_nested
diffbind_results <- template_nested

colnames_diff_nums_rec <- c("rand_comb", "num_of_peaks", "RPKM_diff_regions", "RLEback_diff_regions", "dm6scaled_diff_regions", "greenlist_diff_regions")
diff_nums_rec <- as.data.frame(matrix(nrow = nruns, ncol = length(colnames_diff_nums_rec)))
colnames(diff_nums_rec) <- colnames_diff_nums_rec


for (i in 1:nruns){
  
  diff_nums_rec$rand_comb[i] <-   paste(paste(sample_information_randn[[i]]$Condition[sample_information_randn[[i]]$Factor=="B"],collapse = ","),paste(sample_information_randn[[i]]$Condition[sample_information_randn[[i]]$Factor=="A"],collapse = ","),sep="vs") %>% gsub("million","mill",.)
  
  
  diff_nums_rec$num_of_peaks[i] <- dbObj_counts_summit_list[[i]] %>% dba.peakset(.,bRetrieve = TRUE) %>% as.data.frame() %>% nrow
  
  
  # Library size
  diffbind_results[[i]][["RPKM_diff_regions"]] <- dbObj_RPKM_analyse_list[[i]] %>% dba.report(.,th=0.05) %>% as.data.frame()
  
  diff_nums_rec$RPKM_diff_regions[i] <- diffbind_results[[i]][["RPKM_diff_regions"]] %>% nrow
  
  # RLE Background
  diffbind_results[[i]][["RLEback_diff_regions"]] <- dbObj_RLEback_analyse_list[[i]] %>% dba.report(,th=0.05) %>% as.data.frame()
  
  diff_nums_rec$RLEback_diff_regions[i] <- diffbind_results[[i]][["RLEback_diff_regions"]] %>% nrow
  
  
  # Spike-in
  diffbind_results[[i]][["dm6scaled_diff_regions"]] <- dbObj_dm6_analyse_list[[i]] %>% dba.report(.,th=0.05) %>% as.data.frame()
  
  diff_nums_rec$dm6scaled_diff_regions[i] <- diffbind_results[[i]][["dm6scaled_diff_regions"]] %>% nrow
  
  # RLE RiP
  diffbind_results[[i]][["RLERiP_diff_regions"]] <- dbObj_RLERiP_analyse_list[[i]] %>% dba.report(.,th=0.05) %>% as.data.frame()
  diff_nums_rec$RLERiP_diff_regions[i] <- diffbind_results[[i]][["RLERiP_diff_regions"]] %>% nrow
  
  # Greenlist
  diffbind_results[[i]][["greenlist_diff_regions"]] <- dbObj_greenlist_analyse_list[[i]] %>% dba.report(.,th=0.05) %>% as.data.frame()
  diff_nums_rec$greenlist_diff_regions[i] <- diffbind_results[[i]][["greenlist_diff_regions"]] %>% nrow
  
}



saveRDS(diff_nums_rec, "/data/gpfs/projects/punim2745/03.CUTRUN_normalisation_KA/Henikoff_simulations/R_analysis_response_to_reviewers/05.Negative_experiment/data/diff_nums_rec.rds")

