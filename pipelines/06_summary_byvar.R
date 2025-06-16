#Create a table indicating, for each variant: 
#Number of mutated cell lines
#Number of tumor with the gene as hit
#Number of patients with mut in Intogen
#Number of patients with mut in COSMIC

library(tidyverse)

path_data<-"/data"
path_results<-"/results/20220208"
home<-"/home/aurora.savino/VUS/VUS"

setwd(paste(home, "/", path_results, sep=""))

##load data
tissues<-gsub("_results_ext.RData", "", list.files(pattern="results_ext.RData"))

results<-list()
ind_iter<-0
for(ctiss in tissues){
  ind_iter<-ind_iter+1
  load(paste(ctiss, "_results_ext.RData", sep=""))
  results[[ind_iter]] <- RESTOT
}
names(results)<-tissues

load(file="hits.RData")

intogen_freqs<-read.csv(paste(home, "/", path_data, "/intogen/intogen_summary_5Jul.txt", sep=""), sep="\t")
intogen_freqs<-intogen_freqs[which(intogen_freqs$COUNT>0),]

vars<-read.csv("allvars.csv", header=F)[-1,2]

COSMIC<-read_tsv(paste(home, path_data,"/COSMIC/v96/CosmicGenomeScreensMutantExport.tsv.gz", sep=""))
COSMIC<-data.frame(COSMIC)
ind_cl<-which(COSMIC[,"Sample.Type"] %in% c("cell-line","xenograft"))
COSMIC[,"Gene.name"]<-sub("_.*", "", COSMIC[,"Gene.name"])
COSMIC$Mutation.AA<-gsub("p.", "", COSMIC$Mutation.AA)

tot_COSMIC<-length(unique(COSMIC$ID_tumour[-ind_cl]))

##to undiform COSMIC nomenclature
##remove a letter in COSMIC before fs*
COSMIC$Mutation.AA[grep("fs\\*",COSMIC$Mutation.AA)]<-gsub("[[:upper:]]fs","fs", COSMIC$Mutation.AA[grep("fs\\*",COSMIC$Mutation.AA)]) 
#remove uppercase letters after "del" in hits, but also in intogen
vars[grep("del[[:upper:]]+$",vars)]<-gsub("del[[:upper:]]+$","del",vars[grep("del[[:upper:]]+$",vars)])


summary_vars<-matrix(nrow=length(vars), ncol=10)
rownames(summary_vars)<-vars

all_tiss_intogen <- vector("list", length(vars))
names(all_tiss_intogen) <- vars

all_tiss_cosmic <- vector("list", length(vars))
names(all_tiss_cosmic) <- vars


for(dg in vars){
  gene<-unlist(strsplit(dg, "-"))[[1]]
  mut<-unlist(strsplit(dg, "-"))[[2]]
  summary_vars[dg,1]<-gene
  summary_vars[dg,2]<-mut
  #tessuti where it is a hit
  tissues_hit<-c()
  for(ctiss in tissues){
    if(gene %in% results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2]){
      sel_var<-unlist(strsplit(results[[ctiss]]$var[results[[ctiss]]$GENE==gene & results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2], "\\ | "))
    if(mut %in% gsub("p.", "", unlist(sel_var))){
tissues_hit<-c(tissues_hit, ctiss)
      }   
    }
  }
  summary_vars[dg,3]<-length(tissues_hit)/length(tissues)
  summary_vars[dg,4]<-paste(tissues_hit, collapse=" | ")
  
  #intogen num
  if(length(grep(dg, intogen_freqs$MUTATION))>0){
      ind_dg<-grep(dg, intogen_freqs$MUTATION)
      ind_cohort<-which(intogen_freqs$COHORT=="ALL COHORTS" & intogen_freqs$CANCER_TYPE=="ALL CANCER_TYPE")
      
      summary_vars[dg,5]<-max(intogen_freqs$PERCENTAGE[intersect(ind_dg, ind_cohort)])
  summary_vars[dg,7]<-max(intogen_freqs$COUNT[intersect(ind_dg, ind_cohort)])
  summary_vars[dg,9]<-paste(setdiff(unique(intogen_freqs$CANCER_TYPE[ind_dg]), "ALL CANCER_TYPE"), collapse=" | ")

  ct_counts<-c()
  for(ct in setdiff(unique(intogen_freqs$CANCER_TYPE[ind_dg]), "ALL CANCER_TYPE")){
    
    ind_cohort<-which(intogen_freqs$CANCER_TYPE==ct)
    ct_counts<-c(ct_counts, sum(intogen_freqs$COUNT[intersect(ind_dg, ind_cohort)]))
  }
  if(length(setdiff(unique(intogen_freqs$CANCER_TYPE[ind_dg]), "ALL CANCER_TYPE"))>0){
  names(ct_counts)<-setdiff(unique(intogen_freqs$CANCER_TYPE[ind_dg]), "ALL CANCER_TYPE")
  all_tiss_intogen[[dg]]<-ct_counts
  }
  }

  #COSMIC num
  ind<-which(COSMIC[,"Gene.name"]==gene & COSMIC$Mutation.AA==mut)
  ind_sel<-setdiff(ind, ind_cl)
  summary_vars[dg,10]<-paste(unique(COSMIC$Primary.site[ind_sel]), collapse=" | ")
  summary_vars[dg,6]<-length(unique(COSMIC$ID_tumour[ind_sel]))/tot_COSMIC
  summary_vars[dg,8]<-length(unique(COSMIC$ID_tumour[ind_sel]))

  ct_counts<-c()
  for(ct in unique(COSMIC$Primary.site[ind_sel])){
    ind_tiss<-which(COSMIC$Primary.site==ct)
    ct_counts<-c(ct_counts, length(unique(COSMIC$ID_tumour[intersect(ind_sel, ind_tiss)])))
  }
  if(length(unique(COSMIC$Primary.site[ind_sel]))>0){
  names(ct_counts)<-unique(COSMIC$Primary.site[ind_sel])
  all_tiss_cosmic[[dg]]<-ct_counts
  }
  }

summary_vars<-data.frame(
  Gene=summary_vars[,1],
  Var=summary_vars[,2],
Perc_tiss=as.numeric(summary_vars[,3]),
Tiss=summary_vars[,4],
Perc_Intogen=as.numeric(summary_vars[,5]),
Perc_COSMIC=as.numeric(summary_vars[,6]),
Num_Intogen=as.numeric(summary_vars[,7]),
Num_COSMIC=as.numeric(summary_vars[,8]),
Type_Intogen=(summary_vars[,9]),
Type_COSMIC=(summary_vars[,10])
)
ind_remove<-which(duplicated(vars))

summary_vars<-summary_vars[-ind_remove,]
rownames(summary_vars)<-vars[-ind_remove]
all_tiss_intogen<-all_tiss_intogen[-ind_remove]
all_tiss_cosmic<-all_tiss_cosmic[-ind_remove]

save(summary_vars, file="summary_byvar.RData")
save(all_tiss_intogen, file="all_tiss_intogen_byvar.RData")
save(all_tiss_cosmic, file="all_tiss_comsic_byvar.RData")
