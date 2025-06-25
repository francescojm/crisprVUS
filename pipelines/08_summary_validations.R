library(pheatmap)
library(ggplot2)
library(openxlsx)
library(ReactomePA)
library(clusterProfiler)
library(tidyverse)

path_data<-"/data"
path_results<-"/results/20250221/_DR_plots"
home<-"E:/VUS_2024build"

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

save(summary_drugs, file=paste(home, path_results, "/summary_drugs.RData", sep=""))


                                                                                                                                                         ((summary_vars_bytiss$SIFT %in% c("deleterious", "deleterious_low_confidence"))|(summary_vars_bytiss$PolyPhen %in% c("possibly_damaging", "probably_damaging")))),])))
