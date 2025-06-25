#Create a table indicating, for each variant: 
#Number of mutated cell lines
#Number of tumor with the gene as hit
#Number of patients with mut in Intogen
#Number of patients with mut in COSMIC

library(tidyverse)

path_data<-"/data"
path_results<-"/results/20250221"
home<-"E:/VUS_2024build"

###load patient data
load(file=paste(home, "/", path_results, "/intogen_counts_aggr.RData", sep=""))

COSMIC<-read_tsv(paste(home, path_data,"/raw/Cosmic_GenomeScreensMutant_Tsv_v101_GRCh38/Cosmic_GenomeScreensMutant_v101_GRCh38.tsv.gz", sep=""))
COSMIC<-data.frame(COSMIC)
COSMICsamples<-read_tsv(paste(home, path_data,"/raw/Cosmic_Sample_Tsv_v101_GRCh38/Cosmic_Sample_v101_GRCh38.tsv.gz", sep=""))
COSMICsamples<-data.frame(COSMICsamples)
COSMICclassification<-read_tsv(paste(home, path_data,"/raw/Cosmic_Classification_Tsv_v101_GRCh38/Cosmic_Classification_v101_GRCh38.tsv.gz", sep=""))
COSMICclassification<-data.frame(COSMICclassification)
COSMICprimarysite<-(COSMICclassification$PRIMARY_SITE[match(COSMIC$COSMIC_PHENOTYPE_ID, COSMICclassification$COSMIC_PHENOTYPE_ID)])


##load DAMs data

load(file=paste(home, "/", path_results, "/_allDAMs.RData", sep=""))
vars<-paste(allDAMs$GENE, gsub("p.", "", allDAMs$var), sep="-")
  
COSMIC$MUTATION_AA<-gsub("p.", "", COSMIC$MUTATION_AA)

##to uniform COSMIC nomenclature
##remove a letter in COSMIC before fs*
COSMIC$MUTATION_AA[grep("fs\\*",COSMIC$MUTATION_AA)]<-gsub("[[:upper:]]fs","fs", COSMIC$MUTATION_AA[grep("fs\\*",COSMIC$MUTATION_AA)]) 
#remove uppercase letters after "del" in DAMs, but also in intogen (no del mapped)
vars[grep("del[[:upper:]]+$",vars)]<-gsub("del[[:upper:]]+$","del",vars[grep("del[[:upper:]]+$",vars)])
allDAMs$var<-vars


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
  #tissues where it is a hit
  tissues_hit<-allDAMs$ctype[which(allDAMs$var==dg)]
  summary_vars[dg,3]<-length(tissues_hit)/36 #percentage of tissues
  summary_vars[dg,4]<-paste(tissues_hit, collapse=" | ")
  
  #intogen num
  if(length(which(intogen_counts_aggr$var==dg))>0){
      ind_dg<-which(intogen_counts_aggr$var==dg)
  
   summary_vars[dg,7]<-sum(intogen_counts_aggr[ind_dg, -74])
  summary_vars[dg,9]<-paste(colnames(intogen_counts_aggr[,-74])[which(intogen_counts_aggr[ind_dg,-74]>0)], collapse=" | ")

  ct_counts<-c()
  for(ct_ind in which(intogen_counts_aggr[ind_dg,-74]>0)){
    
      ct_counts<-c(ct_counts, intogen_counts_aggr[ind_dg,-74][ct_ind])
  }
  if(sum(intogen_counts_aggr[ind_dg,-74])){
  names(ct_counts)<-colnames(intogen_counts_aggr[,-74])[which(intogen_counts_aggr[ind_dg,-74]>0)]
  all_tiss_intogen[[dg]]<-ct_counts
  }
  }

  #COSMIC num
  ind_sel<-which(COSMIC$GENE_SYMBOL==gene & COSMIC$MUTATION_AA==mut)
  summary_vars[dg,10]<-paste(unique(COSMICprimarysite[ind_sel]), collapse=" | ")
  summary_vars[dg,8]<-length(unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))

  ct_counts<-c()
  for(ct in unique(COSMICprimarysite[ind_sel])){
    ind_tiss<-which(COSMICprimarysite==ct)
    ct_counts<-c(ct_counts, length(unique(COSMIC$COSMIC_SAMPLE_ID[intersect(ind_sel, ind_tiss)])))
  }
  if(length(unique(COSMICprimarysite[ind_sel]))>0){
  names(ct_counts)<-unique(COSMICprimarysite[ind_sel])
  all_tiss_cosmic[[dg]]<-ct_counts
  }
  }

summary_vars<-data.frame(
  Gene=summary_vars[,1],
  Var=summary_vars[,2],
Perc_tiss=as.numeric(summary_vars[,3]),
Tiss=summary_vars[,4],

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

save(summary_vars, file=paste(home, "/", path_results,"/summary_byvar.RData", sep=""))
save(all_tiss_intogen, file=paste(home, "/", path_results,"/all_tiss_intogen_byvar.RData", sep=""))
save(all_tiss_cosmic, file=paste(home, "/", path_results,"/all_tiss_comsic_byvar.RData", sep=""))
