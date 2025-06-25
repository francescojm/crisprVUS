path_data<-"/data"
path_results<-"/results/20250221"
home<-"E:/VUS_2024build"

setwd(paste(home, "/", path_results, sep=""))

load(file="allDAMs_positions.RData")

intogen_counts<-read.csv(paste(home, "/", path_data, "/raw/only_included_VUS_location_GRCh38_ctype.csv", sep=""), sep="\t")

matches_ind<-match(sub(":[^:]+$", "", intogen_counts[,1]), paste(gsub("chr", "", allvars$chromosome),allvars$position, sep=":"))

allvars$intogen<-NA
allvars$intogen[na.omit(matches_ind)]<-intogen_counts[which(!is.na(matches_ind)),1]

intogen_counts$gene<-allvars$gene_symbol[matches_ind]
intogen_counts$protein_mutation<-allvars$protein_mutation[matches_ind]
intogen_counts<-intogen_counts[which(!is.na(matches_ind)),]

intogen_counts$var<-paste(intogen_counts$gene, gsub("p.", "", intogen_counts$protein_mutation), sep="-")
intogen_counts<-intogen_counts[,c(2:74,84)]
library(dplyr)
intogen_counts<- intogen_counts %>% distinct()

#if various cdna with the same AA change -> sum the rows
intogen_counts_aggr<-as.data.frame(matrix(0, ncol = ncol(intogen_counts), nrow = length(unique(intogen_counts$var))))
colnames(intogen_counts_aggr)<-colnames(intogen_counts)
intogen_counts_aggr$var<-unique(intogen_counts$var)

for(v in unique(intogen_counts$var)){
  if(length(which(intogen_counts$var==v)>1)){
  intogen_counts_aggr[which(intogen_counts_aggr$var==v),-ncol(intogen_counts_aggr)]<-colSums(intogen_counts[which(intogen_counts$var==v),-ncol(intogen_counts)])
  } else {
    intogen_counts_aggr[which(intogen_counts_aggr$var==v),-ncol(intogen_counts_aggr)]<-intogen_counts[which(intogen_counts$var==v),-ncol(intogen_counts)]
}
  }


save(intogen_counts_aggr,file="intogen_counts_aggr.RData")

