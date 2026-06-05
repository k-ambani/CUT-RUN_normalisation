
library(DiffBind)
library(tidyverse)

rm(list=ls())


############################################################################################
# RUNNING 2V2 combinations
############################################################################################


#-------------------------------------------------------------------------------------------------------#
# Bam files not filtered
#-------------------------------------------------------------------------------------------------------#

working_dir <-  "/data/gpfs/projects/punim2745/03.CUTRUN_normalisation_KA/ER_dataset/3.realign_to_human_yeast/"

#loading sample_information workspace
sample_information <- readRDS(file=paste(working_dir,"/R_analysis/1.Sample_information/data/sample_information.rds",sep=""))
sample_information$Replicate <- sample_information$Replicate %>% gsub("R","",.) %>% as.numeric()
sample_information$Factor <- sample_information$Treatment


#---------------------------------------#
# Initialising
#---------------------------------------#

#output objects
nruns = 6
dbObj_first_list <- vector(mode = "list", length = nruns)
dbObj_counts_summit_list<- vector(mode="list", length = nruns)
dbObj_dm6norm_list<- vector(mode="list", length = nruns)
dbObj_sacCer3norm_list<- vector(mode="list", length = nruns)
dbObj_RPKMnorm_list <- vector(mode = "list", length=nruns)
dbObj_RLEbackground_list <- vector(mode="list", length=nruns)
dbObj_RLERiP_list <- vector(mode="list", length=nruns)
dbObj_greenlist_list <- vector(mode="list",length=nruns)

#---------------------------------------#
# Cleaning and reshaping bam count reads for spike-in normalisation
#---------------------------------------#


sample_information_combs = list()
summary_bam_counts_combs = list()
contrast_list = list()



rand_combos_2 <- data.frame(matrix(nrow=nruns,ncol=4))
names(rand_combos_2) <- c("A1","A2","B1","B2")

rand_combos_2[1,1:4] <- c(3,4,1,2)
rand_combos_2[2,1:4] <- c(5,6,1,2)
rand_combos_2[3,1:4] <- c(7,8,1,2)
rand_combos_2[4,1:4] <- c(3,4,5,6)
rand_combos_2[5,1:4] <- c(7,8,3,4)
rand_combos_2[6,1:4] <- c(7,8,5,6)




for (nrun in 1:nruns){
  sample_information_combs[[nrun]] <- sample_information[rand_combos_2[nrun,] %>% as.numeric(),]
  sample_information_combs[[nrun]]$Factor <-  sample_information_combs[[nrun]][["Treatment"]]
  contrast_list[[nrun]] <- sample_information_combs[[nrun]][["Treatment"]] %>% unique()

 
}

# Creating initial dba objects
dbObj_first_list <- dba(sampleSheet = sample_information,minOverlap = 1)


#---------------------------------------#
# Creating dba counts object
#---------------------------------------#
dbObj_counts_summit <- dba.count(dbObj_first_list, summit=200,minOverlap = 1)



#---------------------------------------#
# Normalisations whole experiment
#---------------------------------------#

dbObj_RPKMnorm <- dba.normalize(dbObj_counts_summit)
dbObj_sacCer3norm <- dba.normalize(dbObj_counts_summit, normalize = sample_information$sacCer3_reads_norm_scale)
dbObj_dm6norm <- dba.normalize(dbObj_counts_summit, normalize = sample_information$dm6_reads_norm_scale)
dbObj_RLEbackground <- dba.normalize(dbObj_counts_summit, normalize = "RLE", library = "background", background=TRUE)
dbObj_RLERiP <- dba.normalize(dbObj_counts_summit, normalize = "RLE", library = "RiP") 
dbObj_greenlist <-dba.normalize(dbObj_counts_summit, normalize = sample_information$greenlist_normaliser)

#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING library size NORMALISATION
#--------------------------------------#
# diffbind default params
for (nrun in 1:nruns){
  dbObj_RPKMnorm_list[[nrun]] <- dbObj_RPKMnorm %>%
    dba.contrast(.,contrast = c("Factor",contrast_list[[nrun]][2],contrast_list[[nrun]][1])) %>%
    dba.analyze()

}
#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING sacCer3 spike-in NORMALISATION
#--------------------------------------#

for (nrun in 1:nruns){
dbObj_sacCer3norm_list[[nrun]] <- dbObj_sacCer3norm %>%
  dba.contrast(.,contrast = c("Factor",contrast_list[[nrun]][2],contrast_list[[nrun]][1])) %>%
  dba.analyze()

}
#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING dm6 spike-in NORMALISATION
#--------------------------------------#

for (nrun in 1:nruns){
  dbObj_dm6norm_list[[nrun]] <- dbObj_dm6norm %>%
    dba.contrast(.,contrast = c("Factor",contrast_list[[nrun]][2],contrast_list[[nrun]][1])) %>%
    dba.analyze()
  
}

#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING background reads and RLE
#--------------------------------------#


for (nrun in 1:nruns){
  dbObj_RLEbackground_list[[nrun]] <-dbObj_RLEbackground %>%
    dba.contrast(.,contrast = c("Factor",contrast_list[[nrun]][2],contrast_list[[nrun]][1])) %>%
    dba.analyze()
}

#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING Reads in Peak and RLE
#--------------------------------------#

for (nrun in 1:nruns){
  dbObj_RLERiP_list[[nrun]] <- dbObj_RLERiP %>%
    dba.contrast(.,contrast = c("Factor",contrast_list[[nrun]][2],contrast_list[[nrun]][1])) %>%
    dba.analyze()
}

#---------------------------------------#
# ANALYSING BY CELLLINE-ANTIBODY USING greenlist NORMALISATION
#--------------------------------------#

for (nrun in 1:nruns){
  dbObj_greenlist_list[[nrun]] <- dbObj_greenlist %>%
    dba.contrast(.,contrast = c("Factor",contrast_list[[nrun]][2],contrast_list[[nrun]][1])) %>%
    dba.analyze()
  
}

save.image(file=paste(working_dir,"/R_analysis/2.Diffbind/data/diffbind_analyse_workspace_2v2.RData",sep=""))




