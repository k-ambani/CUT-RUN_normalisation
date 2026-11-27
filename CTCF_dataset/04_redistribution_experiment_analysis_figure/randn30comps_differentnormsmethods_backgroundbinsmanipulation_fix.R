
library(DiffBind)
library(tidyverse)

rm(list=ls())


#---------------------------------------#
# Loading data
#---------------------------------------#

sample_information_randn <- readRDS("/data/gpfs/projects/punim2745/03.CUTRUN_normalisation_KA/Henikoff_simulations/R_analysis_final/02.Sample_information/data/sample_information_randn.rds")



nruns = length(sample_information_randn)
set.seed(10)

#---------------------------------------#
### simulation dials
#---------------------------------------#
#inflation_factor=5
frac_inflate_A = 0.05
frac_inflate_B = 0.05

#---------------------------------------#
### simulation function
#---------------------------------------#



analysing_samples <- function(nrun_choose,inflation_factor) {

sample_information_chosen = sample_information_randn[[nrun_choose]]
dbObj_first <- dba(sampleSheet = sample_information_chosen)

# Creating dba counts object
dbObj_counts_summit_noninflate<- dba.count(dbObj_first, summit=200,minOverlap = 2)
dbObj_counts_summit_inflate <- dbObj_counts_summit_noninflate

dbObj_counts_summit_noninflate$config$cores <- 32
dbObj_counts_summit_inflate$config$cores <- 32


# Inflating random peaks
num_peaks <- dbObj_counts_summit_noninflate %>% dba.peakset(bRetrieve = TRUE) %>% as.data.frame() %>% nrow


num_peaks_inflate_A <- floor(num_peaks*frac_inflate_A)
num_peaks_inflate_B <- floor(num_peaks*frac_inflate_B)

n_peaks_inflate_A_B =sample.int(num_peaks,num_peaks_inflate_A + num_peaks_inflate_B)

n_peaks_inflate_A <- n_peaks_inflate_A_B[1:num_peaks_inflate_A]
n_peaks_inflate_B <- n_peaks_inflate_A_B[(num_peaks_inflate_A+1):length(n_peaks_inflate_A_B)]
  

for(n_peak_A in n_peaks_inflate_A){
  
  dbObj_counts_summit_inflate[["peaks"]][[1]][n_peak_A,c("Score","RPKM","Reads","cRPKM","cReads")] <- dbObj_counts_summit_noninflate[["peaks"]][[1]][n_peak_A,c("Score","RPKM","Reads","cRPKM","cReads")]*inflation_factor
  dbObj_counts_summit_inflate[["peaks"]][[2]][n_peak_A,c("Score","RPKM","Reads","cRPKM","cReads")] <- dbObj_counts_summit_noninflate[["peaks"]][[2]][n_peak_A,c("Score","RPKM","Reads","cRPKM","cReads")]*inflation_factor
  dbObj_counts_summit_inflate[["peaks"]][[3]][n_peak_A,c("Score","RPKM","Reads","cRPKM","cReads")] <- dbObj_counts_summit_noninflate[["peaks"]][[3]][n_peak_A,c("Score","RPKM","Reads","cRPKM","cReads")]*inflation_factor
  
}

for(n_peak_B in n_peaks_inflate_B){
  
  dbObj_counts_summit_inflate[["peaks"]][[4]][n_peak_B,c("Score","RPKM","Reads","cRPKM","cReads")] <- dbObj_counts_summit_noninflate[["peaks"]][[4]][n_peak_B,c("Score","RPKM","Reads","cRPKM","cReads")]*inflation_factor
  dbObj_counts_summit_inflate[["peaks"]][[5]][n_peak_B,c("Score","RPKM","Reads","cRPKM","cReads")] <- dbObj_counts_summit_noninflate[["peaks"]][[5]][n_peak_B,c("Score","RPKM","Reads","cRPKM","cReads")]*inflation_factor
  dbObj_counts_summit_inflate[["peaks"]][[6]][n_peak_B,c("Score","RPKM","Reads","cRPKM","cReads")] <- dbObj_counts_summit_noninflate[["peaks"]][[6]][n_peak_B,c("Score","RPKM","Reads","cRPKM","cReads")]*inflation_factor
  
}

dbObj_counts_summit_inflate<- dba(dbObj_counts_summit_inflate) 



# Saving inflated and noninflated peakset to for sanity check

inflate_df <- dbObj_counts_summit_inflate %>% dba.peakset(bRetrieve = TRUE) %>% as.data.frame()
noninflate_df <- dbObj_counts_summit_noninflate %>% dba.peakset(bRetrieve = TRUE) %>% as.data.frame() 

A1_inflate_counts <- dbObj_counts_summit_inflate[["peaks"]][[1]] %>% as.data.frame()
colnames(A1_inflate_counts)[4:8] <- paste("A1",colnames(A1_inflate_counts)[4:8],sep="_")

A2_inflate_counts <- dbObj_counts_summit_inflate[["peaks"]][[2]] %>% as.data.frame()
colnames(A2_inflate_counts)[4:8] <- paste("A2",colnames(A2_inflate_counts)[4:8],sep="_")

A3_inflate_counts <- dbObj_counts_summit_inflate[["peaks"]][[3]] %>% as.data.frame()
colnames(A3_inflate_counts)[4:8] <- paste("A3",colnames(A3_inflate_counts)[4:8],sep="_")


B1_inflate_counts <- dbObj_counts_summit_inflate[["peaks"]][[4]] %>% as.data.frame()
colnames(B1_inflate_counts)[4:8] <- paste("B1",colnames(B1_inflate_counts)[4:8],sep="_")

B2_inflate_counts <- dbObj_counts_summit_inflate[["peaks"]][[5]] %>% as.data.frame()
colnames(B2_inflate_counts)[4:8] <- paste("B2",colnames(B2_inflate_counts)[4:8],sep="_")

B3_inflate_counts <- dbObj_counts_summit_inflate[["peaks"]][[6]] %>% as.data.frame()
colnames(B3_inflate_counts)[4:8] <- paste("B3",colnames(B3_inflate_counts)[4:8],sep="_")




inflate_rec <- full_join(A1_inflate_counts,A2_inflate_counts, by = c("Chr", "Start", "End"), relationship = "one-to-one")

inflate_rec <- full_join(inflate_rec,A3_inflate_counts, by = c("Chr", "Start", "End"), relationship = "one-to-one")

inflate_rec <- full_join(inflate_rec,B1_inflate_counts, by = c("Chr", "Start", "End"), relationship = "one-to-one")

inflate_rec <- full_join(inflate_rec,B2_inflate_counts, by = c("Chr", "Start", "End"), relationship = "one-to-one")

inflate_rec <- full_join(inflate_rec,B3_inflate_counts, by = c("Chr", "Start", "End"), relationship = "one-to-one")

columns_to_remove <- c(grep("cRPKM", names(inflate_rec)),grep("cReads", names(inflate_rec)))
inflate_rec <- inflate_rec[,-columns_to_remove]



## For non-inflate


A1_noninflate_counts <- dbObj_counts_summit_noninflate[["peaks"]][[1]] %>% as.data.frame()
colnames(A1_noninflate_counts)[4:8] <- paste("A1",colnames(A1_noninflate_counts)[4:8],sep="_")

A2_noninflate_counts <- dbObj_counts_summit_noninflate[["peaks"]][[2]] %>% as.data.frame()
colnames(A2_noninflate_counts)[4:8] <- paste("A2",colnames(A2_noninflate_counts)[4:8],sep="_")

A3_noninflate_counts <- dbObj_counts_summit_noninflate[["peaks"]][[3]] %>% as.data.frame()
colnames(A3_noninflate_counts)[4:8] <- paste("A3",colnames(A3_noninflate_counts)[4:8],sep="_")


B1_noninflate_counts <- dbObj_counts_summit_noninflate[["peaks"]][[4]] %>% as.data.frame()
colnames(B1_noninflate_counts)[4:8] <- paste("B1",colnames(B1_noninflate_counts)[4:8],sep="_")

B2_noninflate_counts <- dbObj_counts_summit_noninflate[["peaks"]][[5]] %>% as.data.frame()
colnames(B2_noninflate_counts)[4:8] <- paste("B2",colnames(B2_noninflate_counts)[4:8],sep="_")

B3_noninflate_counts <- dbObj_counts_summit_noninflate[["peaks"]][[6]] %>% as.data.frame()
colnames(B3_noninflate_counts)[4:8] <- paste("B3",colnames(B3_noninflate_counts)[4:8],sep="_")




noninflate_rec <- full_join(A1_noninflate_counts,A2_noninflate_counts, by = c("Chr", "Start", "End"), relationship = "one-to-one")


noninflate_rec <- full_join(noninflate_rec,A3_noninflate_counts, by = c("Chr", "Start", "End"), relationship = "one-to-one")

noninflate_rec <- full_join(noninflate_rec,B1_noninflate_counts, by = c("Chr", "Start", "End"), relationship = "one-to-one")

noninflate_rec <- full_join(noninflate_rec,B2_noninflate_counts, by = c("Chr", "Start", "End"), relationship = "one-to-one")

noninflate_rec <- full_join(noninflate_rec,B3_noninflate_counts, by = c("Chr", "Start", "End"), relationship = "one-to-one")

columns_to_remove <- c(grep("cRPKM", names(noninflate_rec)),grep("cReads", names(noninflate_rec)))
noninflate_rec <- noninflate_rec[,-columns_to_remove]



#---------------------------------------#
# Getting additional number of counts from inflation
#---------------------------------------#

inflate_counts <- inflate_rec[,grep("Read", colnames(inflate_rec))]

noninflate_counts <- noninflate_rec[,grep("Read", colnames(noninflate_rec))]


diff_counts <- inflate_counts - noninflate_counts


#---------------------------------------#
# Analysing with Library Size
#---------------------------------------#

dbObj_LS_noninflate <- dba.normalize(dbObj_counts_summit_noninflate)  %>% dba.contrast() %>% dba.analyze() %>% dba.report(.,th=1) %>% as.data.frame()

get_lib_sizes <- dba.normalize(dbObj_counts_summit_inflate) 
non_inflate_lib_sizes <- get_lib_sizes[["norm"]][["DESeq2"]][["lib.sizes"]]
new_lib_sizes <- c(non_inflate_lib_sizes + colSums(diff_counts)) %>% unname()


dbObj_LS_inflate <- dba.normalize(dbObj_counts_summit_inflate, normalize = new_lib_sizes/mean(new_lib_sizes)) %>% dba.contrast() %>% dba.analyze() %>% dba.report(.,th=1) %>% as.data.frame()





#---------------------------------------#
# Analysing with greenlist normalisation method
#--------------------------------------#

dbObj_greenlist_inflate <- dba.normalize(dbObj_counts_summit_inflate, normalize = sample_information_chosen$greenlist_normalizer)  %>% dba.contrast() %>% dba.analyze() %>% dba.report(.,th=1) %>% as.data.frame()
dbObj_greenlist_noninflate <- dba.normalize(dbObj_counts_summit_noninflate, normalize = sample_information_chosen$greenlist_normalizer)  %>% dba.contrast() %>% dba.analyze() %>% dba.report(.,th=1) %>% as.data.frame()


#---------------------------------------#
# Background normalisation - preparing dba object
#---------------------------------------#

# Creating background bins sub-object - doesn't matter which normalize (RLE or TMM) parameter, just getting the background bins counts
dbObj_background_inflate <- dba.normalize(dbObj_counts_summit_inflate, normalize = "RLE", library = "background", background=TRUE) 


background_bins <- dbObj_background_inflate[["norm"]][["background"]][["binned"]]@rowRanges  %>% as.data.frame()

background_counts <- dbObj_background_inflate[["norm"]][["background"]][["binned"]]@assays@data@listData[["counts"]] %>% as.data.frame()



background_gr <- data.frame(background_bins,background_counts) %>% GRanges()



raw_counts_mat <- dbObj_background_inflate[["norm"]][["background"]][["binned"]]@assays@data@listData
#granges_rawcounts <- dbObj_LS_inflate %>% dba.peakset(.,bRetrieve = TRUE) %>% as.data.frame()


noninflate_reads <- noninflate_rec[,grep("Reads",colnames(noninflate_rec))]

inflate_reads <- inflate_rec[,grep("Reads",colnames(inflate_rec))]



reads_diff = inflate_reads - noninflate_reads
reads_diff <- data.frame(inflate_rec[,c("Chr", "Start", "End")],reads_diff)
reads_diff = reads_diff[rowSums(reads_diff[,c(4:9)])>0,] %>% GRanges()
reads_diff_df <- reads_diff %>% as.data.frame()


overlap_rec <- findOverlaps(reads_diff,background_gr) %>% as.data.frame()

reads_diff_mat <- reads_diff%>% as.data.frame()
reads_diff_mat <- reads_diff_mat[,6:11]

background_counts_inflate <- background_counts

background_counts_inflate[overlap_rec$subjectHits,] <- background_counts[overlap_rec$subjectHits,] + reads_diff_mat[overlap_rec$queryHits,]


dbObj_background_inflate[["norm"]][["background"]][["binned"]]@assays@data@listData[["counts"]] <- background_counts_inflate %>% as.matrix # inserting inflated background reads into dba object


#---------------------------------------#
# Analysing with RLE background normalisation
#---------------------------------------#

#Renormalising with new background bins and then carrying out downstream analysis
dbObj_RLEbackground_inflate <-  dba.normalize(dbObj_background_inflate, normalize = "RLE", library = "background", background=TRUE)  %>% dba.contrast() %>% dba.analyze() %>% dba.report(.,th=1) %>% as.data.frame() 


dbObj_RLEbackground_noninflate <- dba.normalize(dbObj_counts_summit_noninflate, normalize = "RLE", library = "background", background=TRUE)  %>% dba.contrast() %>% dba.analyze() %>% dba.report(.,th=1) %>% as.data.frame()


#---------------------------------------#
# Analysing with RLE RiP
#---------------------------------------#

dbObj_RLE_RiP_inflate <- dba.normalize(dbObj_counts_summit_inflate, normalize = "RLE", library = "RiP") %>% dba.contrast() %>% dba.analyze() %>% dba.report(.,th=1) %>% as.data.frame()
dbObj_RLE_RiP_noninflate <- dba.normalize(dbObj_counts_summit_noninflate, normalize = "RLE", library = "RiP")  %>% dba.contrast() %>% dba.analyze() %>% dba.report(.,th=1) %>% as.data.frame()

#---------------------------------------#
# Analysing with spike-in normalisation
#---------------------------------------#


dbObj_spikein_noninflate <- dba.normalize(dbObj_counts_summit_noninflate, normalize = sample_information_chosen$dm6_reads_norm_scale)  %>% dba.contrast() %>% dba.analyze() %>% dba.report(.,th=1) %>% as.data.frame()


# Recalculating dm6 scale factors because of the extra reads due to inflating peaks
# new-spikein reads  = new library size * old-spikein reads / old library size

sample_information_chosen$new_dm6_intermediate_bam_counts <- (new_lib_sizes)* sample_information_chosen$dm6_intermediate_bam_counts/ non_inflate_lib_sizes
sample_information_chosen$new_dm6_reads_norm_scale = sample_information_chosen$new_dm6_intermediate_bam_counts/mean(sample_information_chosen$new_dm6_intermediate_bam_counts)

dbObj_spikein_inflate <- dba.normalize(dbObj_counts_summit_inflate, normalize = sample_information_chosen$new_dm6_reads_norm_scale)  %>% dba.contrast() %>% dba.analyze() %>% dba.report(.,th=1) %>% as.data.frame()


#---------------------------------------#
# Putting results in a list
#---------------------------------------#

results <- list(dbObj_first, reads_diff_df, n_peaks_inflate_A, n_peaks_inflate_B,
                noninflate_rec, inflate_rec,
                dbObj_counts_summit_noninflate, dbObj_counts_summit_inflate, 
                dbObj_LS_noninflate, dbObj_LS_inflate, 
                dbObj_RLEbackground_noninflate, dbObj_RLEbackground_inflate,
                dbObj_greenlist_noninflate, dbObj_greenlist_inflate,
                dbObj_RLE_RiP_noninflate, dbObj_RLE_RiP_inflate,
                dbObj_spikein_noninflate, dbObj_spikein_inflate
                )

names(results) <- c("dbObj_first", "reads_diff", "A_peaks_inflated", "B_peaks_inflated",
                    "noninflate_rec", "inflate_rec", 
                    "dbObj_counts_summit_noninflate", "dbObj_counts_summit_inflate", 
                "LS_noninflate_df", "LS_inflate_df",
                "RLE_BG_noninflate_df", "RLE_BG_inflate_df",
                "greenlist_noninflate_df", "greenlist_inflate_df",
                "RLE_RiP_noninflate_df", "RLE_RiP_inflate_df",
                "spikein_noninflate_df", "spikein_inflate_df"
                )


return(results)
  
}



rand_runs <- paste("rand",1:nruns,".rds",sep="")

pval_threshold <- 0.05


for(inflation_factor_choose in c(2,5,10,20)){

    
  working_dir <-  paste("/data/gpfs/projects/punim2745/03.CUTRUN_normalisation_KA/Henikoff_simulations/R_analysis_response_to_reviewers/02.Artificially_inflate_incl_spikein/results/inflation_factor",as.character(inflation_factor_choose),sep="_")

  print(inflation_factor_choose)
  
  
  inflated_DESeq_record <- list()
  
  experiment_summary_nrand <- data.frame(rand_runs=paste("rand",c(1:nruns),sep="_"))
  
for(nrun in 1:nruns){
  print(nrun)
  randn <- analysing_samples(nrun,inflation_factor = inflation_factor_choose)
  #randn = rand1
saveRDS(randn,file = paste(working_dir,"individual_runs",rand_runs[nrun],sep="/"))


#---------------------------------------#
# Summary of results table
#---------------------------------------#


rand_int <- paste("rand",nrun,sep="")

# Saving meta information about experiment
experiment_summary_nrand[nrun,"groupA_samples"] <- randn[["dbObj_first"]][["masks"]][["A"]][randn[["dbObj_first"]][["masks"]][["A"]]] %>% names %>% paste(.,collapse = ", ")
experiment_summary_nrand[nrun,"groupB_samples"] <- randn[["dbObj_first"]][["masks"]][["B"]][randn[["dbObj_first"]][["masks"]][["B"]]] %>% names %>% paste(.,collapse = ", ")
experiment_summary_nrand[nrun,"inflation_factor"] <- inflation_factor_choose
experiment_summary_nrand[nrun,"npeaks_consensus"] <- randn[["noninflate_rec"]] %>% nrow
experiment_summary_nrand[nrun,"npeaks_inflated_groupA"] <- randn[["A_peaks_inflated"]] %>% length
experiment_summary_nrand[nrun,"npeaks_inflated_groupB"] <- randn[["B_peaks_inflated"]] %>% length

norm_methods <- c("LS", "RLE_BG", "greenlist", "RLE_RiP", "spikein")

for(norm_method in norm_methods){
  # Diff peaks
  experiment_summary_nrand[nrun,paste0(norm_method,"_NI_diffpeaks")] <- c(randn[[paste0(norm_method,"_noninflate_df")]]$FDR < pval_threshold) %>% sum
  experiment_summary_nrand[nrun,paste0(norm_method,"_NI_uppeaks")] <- (c(randn[[paste0(norm_method,"_noninflate_df")]]$FDR < pval_threshold) &  c(randn[[paste0(norm_method,"_noninflate_df")]]$Fold>0)) %>% sum
  experiment_summary_nrand[nrun,paste0(norm_method,"_NI_downpeaks")] <- (c(randn[[paste0(norm_method,"_noninflate_df")]]$FDR < pval_threshold) &  c(randn[[paste0(norm_method,"_noninflate_df")]]$Fold<0)) %>% sum
  
  experiment_summary_nrand[nrun,paste0(norm_method,"_I_diffpeaks")] <- c(randn[[paste0(norm_method,"_inflate_df")]]$FDR < pval_threshold) %>% sum
  experiment_summary_nrand[nrun,paste0(norm_method,"_I_uppeaks")] <- (c(randn[[paste0(norm_method,"_inflate_df")]]$FDR < pval_threshold) &  c(randn[[paste0(norm_method,"_inflate_df")]]$Fold>0)) %>% sum
  experiment_summary_nrand[nrun,paste0(norm_method,"_I_downpeaks")] <- (c(randn[[paste0(norm_method,"_inflate_df")]]$FDR < pval_threshold) &  c(randn[[paste0(norm_method,"_inflate_df")]]$Fold<0)) %>% sum
                                                                                                                                                
  # TRUE/FALSE POSITIVES/NEGATIVES
  DESeq_inflate <- left_join(randn[[paste0(norm_method,"_inflate_df")]], randn[["reads_diff"]])
  DESeq_inflate$A_inflate <- c(DESeq_inflate[,c("A1_Reads","A2_Reads","A3_Reads")] %>% rowSums())>0
  DESeq_inflate$B_inflate <- c(DESeq_inflate[,c("B1_Reads","B2_Reads","B3_Reads")] %>% rowSums())>0
  DESeq_inflate$inflate <- !is.na(DESeq_inflate$A_inflate | DESeq_inflate$B_inflate)
  
  DESeq_inflate$TP <- c(DESeq_inflate$FDR<0.05 & DESeq_inflate$inflate)
  DESeq_inflate$FP <- c(DESeq_inflate$FDR<0.05 & !DESeq_inflate$inflate)
  DESeq_inflate$FN <- c(DESeq_inflate$FDR>0.05 & DESeq_inflate$inflate)
  DESeq_inflate$TN <- c(DESeq_inflate$FDR>0.05 & !DESeq_inflate$inflate)
  
  inflated_DESeq_record[[rand_int]][[paste0(norm_method,"_DESeq_inflate")]] <- DESeq_inflate
  
  experiment_summary_nrand[nrun,paste0(norm_method,"_I_TP")] <- DESeq_inflate$TP %>% na.omit %>% sum
  experiment_summary_nrand[nrun,paste0(norm_method,"_I_FP")] <- DESeq_inflate$FP %>% na.omit %>% sum
  experiment_summary_nrand[nrun,paste0(norm_method,"_I_FN")] <- DESeq_inflate$FN %>% na.omit %>% sum
  experiment_summary_nrand[nrun,paste0(norm_method,"_I_TN")] <- DESeq_inflate$TN %>% na.omit %>% sum
  
  
}

}
  saveRDS(experiment_summary_nrand,file = paste(working_dir,"experiment_summary.rds",sep="/"))
  saveRDS(inflated_DESeq_record,file = paste(working_dir,"inflated_DESeq_record.rds",sep="/"))
  
}
