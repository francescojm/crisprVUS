library(pheatmap)
library(ggplot2)
library(openxlsx)
library(ReactomePA)
library(clusterProfiler)
library(tidyverse)

path_data<-"/data"
path_results<-"/results/20220208/_DR_plots"
home<-"/Volumes/home/VUS/VUS/"

setwd(paste(home, "/", path_results, sep=""))

##load data
tissues_drugs<-gsub("_DR_validation.RData", "", list.files(pattern="validation.RData"))

results_drugs<-list()
ind_iter<-0
for(ctiss in tissues_drugs){
  ind_iter<-ind_iter+1
  load(paste(ctiss, "_DR_validation.RData", sep=""))
  results_drugs[[ind_iter]] <- rRES
}
names(results_drugs)<-tissues_drugs


genes_tot<-c()
vars_tot<-c()
tiss_tot<-c()
drugs_tot<-c()
for(i in 1:length(results_drugs)){
  vars_list<-strsplit(results_drugs[[i]]$Var[results_drugs[[i]]$bestRR<1.6], " \\| ")
  genes<-rownames(results_drugs[[i]])[results_drugs[[i]]$bestRR<1.6]
  drugs<-results_drugs[[i]]$bestRR_drug_name[results_drugs[[i]]$bestRR<1.6]
  names(vars_list)<-genes
  
  genes_tot<-c(genes_tot, unlist(lapply(genes, function(x){rep(x,length(unlist(vars_list[x])))})))
  vars_tot<-c(vars_tot, gsub("p.", "", unlist(vars_list)))
  tiss_tot<-c(tiss_tot, rep(tissues_drugs[i], length(unlist(vars_list))))
  names(vars_list)<-drugs
  drugs_tot<-c(drugs_tot, unlist(lapply(drugs, function(x){rep(x,length(unlist(vars_list[x])))})))
  
}

summary_drugs<-data.frame(genes=genes_tot, vars=vars_tot, tiss=tiss_tot, drugs_tot)
summary_drugs$var_id<-paste(summary_drugs$genes, summary_drugs$vars, sep="-")


View(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id),])
View(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id & summary_vars_bytiss$tissue=="ALL" & summary_vars_bytiss$database=="COSMIC"),])

print(paste("Number of SAMs that are also reliable DAMs", nrow(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id & summary_vars_bytiss$tissue=="ALL" & summary_vars_bytiss$database=="COSMIC"),])))

print(paste("Number of SAMs that are also reliable DAMs in at least one patient", nrow(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id & summary_vars_bytiss$tissue=="ALL" & summary_vars_bytiss$database=="COSMIC" & summary_vars_bytiss$num_patients>0 ),])))

setwd("/Volumes/home/VUS/VUS/results/20220208")

summary_vars_bytiss$driver<-ifelse(summary_vars_bytiss$Gene %in% driver_genes, "Driver", "Non-driver")
summary_vars_bytiss$SIFT<-vars_vep$transcript_consequences_sift_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]
summary_vars_bytiss$PolyPhen<-vars_vep$transcript_consequences_polyphen_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]

summary_vars_bytiss$SIFT<-gsub("\\(.*\\)", "", summary_vars_bytiss$SIFT)
summary_vars_bytiss$PolyPhen<-gsub("\\(.*\\)", "", summary_vars_bytiss$PolyPhen)
summary_vars_bytiss$tot_patients<-summary_vars_bytiss$Num_Intogen+summary_vars_bytiss$Num_COSMIC
summary_vars_bytiss<-summary_vars_bytiss[,c("ID","Gene","Var","driver","Perc_tiss", "Tiss", "tot_patients", "Num_Intogen", "Num_COSMIC",
                                      "Perc_Intogen", "Perc_COSMIC",
                                      "Num_tiss_Intogen", "Num_tiss_COSMIC", 
                                      "Type_Intogen", "Type_COSMIC",
                                      "Match_tiss_Intogen", "Match_tiss_COSMIC",
                                      "Match_tiss_Intogen_name", "Match_tiss_COSMIC_name", "Matching_inboth", "Matching_union",
                                      "SIFT", "PolyPhen", "tissue")]
summary_vars_bytiss<-summary_vars_bytiss[order(summary_vars_bytiss$tot_patients, decreasing = T),]

write.xlsx(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id),], "vars_depANDdrug.xlsx")

print(paste("Number of SAMs that are also reliable DAMs in at least one patient (same type)", nrow(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id & summary_vars_bytiss$tissue=="ALL_MATCHING" & summary_vars_bytiss$database=="COSMIC" & summary_vars_bytiss$num_patients>0 ),])))

print(paste("Number of SAMs that are also reliable DAMs in at least one patient (same type) AND predicted deleterious", nrow(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id & summary_vars_bytiss$tissue=="ALL_MATCHING" & summary_vars_bytiss$database=="COSMIC" & summary_vars_bytiss$num_patients>0 &
                                                                                                                                                         ((summary_vars_bytiss$SIFT %in% c("deleterious", "deleterious_low_confidence"))|(summary_vars_bytiss$PolyPhen %in% c("possibly_damaging", "probably_damaging")))),])))
