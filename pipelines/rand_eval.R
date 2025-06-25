library(tidyverse)
#set path
path_data<-"/data"
path_results<-"/results/20250221"
home<-"E:/VUS_2024build"

setwd(paste(home, "/", path_results, sep=""))

load("hits_nop.RData")

###############################
### load randomization results
################################

rand_files<-list.files(pattern="rand", recursive = T)
genes_rand<-matrix(nrow=length(hits_nop), ncol=length(rand_files), dimnames=list(c(hits_nop), c(gsub("_results_rand.RData", "", rand_files))))

for(file in rand_files){
  print(file)
  load(file)
  
  genes<-c()
  for(i in 1:1000){
    print(i)
    genes<-c(genes,RESTOT_rand[[i]][[1]]$GENE[which(RESTOT_rand[[i]][[1]]$rank_ratio<1.6 & RESTOT_rand[[i]][[1]]$medFitEff< -.5)])
  }
  
  genes_rand[intersect(names(table(genes)), hits_nop), gsub("VUS_rand/20250221/|_results_rand.RData", "", file)]<-table(genes)[intersect(names(table(genes)), hits_nop)]
}

genes_rand[is.na(genes_rand)]<-0
save(genes_rand, file="genes_rand.RData")


##load main results
tissues<-gsub("_results.RData", "", list.files(pattern="results.RData"))

results<-list()
ind_iter<-0
for(ctiss in tissues){
  ind_iter<-ind_iter+1
  load(paste(ctiss, "_results.RData", sep=""))
  results[[ind_iter]] <- RESTOT
}
names(results)<-tissues

###############################
### compute the empirical pvalue based on the randomizations
################################

summary_hits<-matrix(0, nrow=length(hits_nop), ncol=length(tissues))
rownames(summary_hits)<-hits_nop
colnames(summary_hits)<-tissues

for(dg in hits_nop){
  for(ctiss in tissues){
    if(dg %in% results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5]){
      summary_hits[dg,ctiss]<-1
    }
    
  }
}

summary_pval<-genes_rand/1000

for(ctiss in tissues){
   load(paste(ctiss, "_results.RData", sep=""))
   RESTOT$pval_rand<-NA
   incommon<-intersect(rownames(summary_pval), RESTOT$GENE)
  RESTOT[match(incommon,RESTOT$GENE), "pval_rand"]<-summary_pval[incommon,ctiss]
  save(RESTOT, file=paste(ctiss, "_results_ext.RData", sep=""))
}

##load DAM data with p
tissues<-gsub("_results_ext.RData", "", list.files(pattern="results_ext.RData"))

results<-list()
ind_iter<-0
for(ctiss in tissues){
  ind_iter<-ind_iter+1
  load(paste(ctiss, "_results_ext.RData", sep=""))
  results[[ind_iter]] <- RESTOT
}
names(results)<-tissues

hits<-c()
for(ctiss in tissues){
  hits<-c(hits, unique(results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2 ]))
}

hits<-unique(hits)
save(hits, file="hits.RData")


