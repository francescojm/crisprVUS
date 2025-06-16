setwd("/home/aurora.savino/VUS/VUS/results/20220208/rand_hits/")
load("/home/aurora.savino/VUS/VUS/results/20220208/hits_nop.RData")

rand_files<-list.files(pattern="rand")
genes_rand<-matrix(nrow=length(hits), ncol=length(rand_files), dimnames=list(c(hits), c(gsub("_results_rand.RData", "", rand_files))))

for(file in rand_files){
  print(file)
  load(file)
  
  genes<-c()
  for(i in 1:1000){
    print(i)
    genes<-c(genes,RESTOT_rand[[i]][[1]]$GENE[which(RESTOT_rand[[i]][[1]]$rank_ratio<1.6 & RESTOT_rand[[i]][[1]]$medFitEff< -.5)])
  }
  
  genes_rand[names(table(genes)), gsub("_results_rand.RData", "", file)]<-table(genes)
}
genes_rand[is.na(genes_rand)]<-0
save(genes_rand, file="genes_rand.RData")


###################
library(pheatmap)
library(ggplot2)
library(openxlsx)
library(ReactomePA)
library(clusterProfiler)
library(tidyverse)

path_data<-"/data"
path_results<-"/results/20220208"
home<-"/home/aurora.savino/VUS/VUS"

setwd(paste(home, "/", path_results, sep=""))

##load dei dati
tissues<-gsub("_results.RData", "", list.files(pattern="results.RData"))

results<-list()
ind_iter<-0
for(ctiss in tissues){
  ind_iter<-ind_iter+1
  load(paste(ctiss, "_results.RData", sep=""))
  results[[ind_iter]] <- RESTOT
}
names(results)<-tissues

summary_hits<-matrix(0, nrow=length(hits), ncol=length(tissues))
rownames(summary_hits)<-hits
colnames(summary_hits)<-tissues

for(dg in hits){
  for(ctiss in tissues){
    if(dg %in% results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5]){
      summary_hits[dg,ctiss]<-1
    }
    
  }
}

summary_pval<-genes_rand/1000
summary_fdr<-apply(summary_pval,2,function(x){p.adjust(x,method="fdr")})
hist(-log10(summary_fdr[summary_hits==1]), breaks=100)

hits_filt<-rownames(summary_fdr)[which(summary_fdr<0.1 & summary_hits==1, arr.ind=T)[,1]]

sort(table(hits_filt))

for(ctiss in tissues){
   load(paste(ctiss, "_results.RData", sep=""))
   RESTOT$pval_rand<-NA
   incommon<-intersect(rownames(summary_pval), RESTOT$GENE)
  RESTOT[match(incommon,RESTOT$GENE), "pval_rand"]<-summary_pval[incommon,ctiss]
  save(RESTOT, file=paste(ctiss, "_results_ext.RData", sep=""))
}


