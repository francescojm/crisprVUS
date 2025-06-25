library(reshape2)
library(ggplot2)
library(openxlsx)

#load("hits.RData")

#####################################################################
################################### BY VAR
#########################################################################

#######################
##matching cancer types
##################
home<-"E:/VUS_2024build/"
path_data<-"/data"
path_results<-"/results/20250221"

##driver genes
inTOgen_drivers<-read.table(paste(home, '/data/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
### inTOgene drivers downloaded from https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip on 20241002
driver_genes<-unique(inTOgen_drivers$SYMBOL)

setwd(paste(home, "/", path_results, sep=""))

load("summary_byvar.RData")

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

#load("results_with_expr.RData")

#library(tidyverse)
#gene_annot <- read_csv(paste(home, "/",path_data, "/gene_identifiers_20191101.csv", sep=""))

#CMP_annot <- read_csv(paste(home, "/",path_data,"/model_list_20210611.csv", sep="")) # from https://cog.sanger.ac.uk/cmp/download/model_list_20210611.csv

mapping<-read.csv(paste(home, path_data,'/raw/intOGen ctype mapping_AS.csv',sep=''),header = TRUE,row.names = 1, sep=";")

cancer_match_long_CMP<-c()
cancer_match_long_into<-c()
for(i in 1:nrow(mapping)){
  cancer_match_long_into<-c(cancer_match_long_into, unlist(strsplit(mapping[i,1], " \\| ")))
  cancer_match_long_CMP<-c(cancer_match_long_CMP, rep(rownames(mapping)[i],length(unlist(strsplit(mapping[i,1], " \\| ")))))
}
mapping_intogen<-cbind(cancer_match_long_CMP, cancer_match_long_into)

cancer_match_long_CMP<-c()
cancer_match_long_into<-c()
for(i in 1:nrow(mapping)){
  cancer_match_long_into<-c(cancer_match_long_into, unlist(strsplit(mapping[i,2], " \\| ")))
  cancer_match_long_CMP<-c(cancer_match_long_CMP, rep(rownames(mapping)[i],length(unlist(strsplit(mapping[i,2], " \\| ")))))
}
mapping_cosmic<-cbind(cancer_match_long_CMP, cancer_match_long_into)

load(file=paste(home, path_results, "summary_drivers_bin.RData", sep=""))

###########################
########### for each variant, indicate the number of tissues in intogen, cosmic
############################

#################
#### each cancer type separately
###################
load("summary_byvar.RData")
summary_vars$ID<-paste(summary_vars$Gene, summary_vars$Var, sep="-")
summary_vars$Num_tiss<-summary_vars$Perc_tiss*36
summary_vars$Num_Intogen[is.na(summary_vars$Num_Intogen)]<-0

summary_vars$Num_tiss_Intogen<-NA
summary_vars$Num_tiss_COSMIC<-NA
summary_vars$Match_tiss_Intogen<-NA
summary_vars$Match_tiss_COSMIC<-NA
summary_vars$Matching_inboth<-NA
summary_vars$Matching_union<-NA
summary_vars$Match_tiss_Intogen_name<-NA
summary_vars$Match_tiss_COSMIC_name<-NA

load("all_tiss_comsic_byvar.RData")
load("all_tiss_intogen_byvar.RData")

for(i in 1:nrow(summary_vars)){
  type_tiss<-unique(unlist(strsplit(summary_vars[i, "Tiss"], " \\| ")))
  type_intogen<-na.omit(unique(unlist(strsplit(summary_vars[i, "Type_Intogen"], " \\| "))))
  type_cosmic<-na.omit(unique(unlist(strsplit(summary_vars[i, "Type_COSMIC"], " \\| "))))
  
  summary_vars$Num_tiss_Intogen[i]<-length(type_intogen)
  summary_vars$Num_tiss_COSMIC[i]<-length(type_cosmic)
  
  matching_intogen<-unique(mapping_intogen[which(mapping_intogen[,2] %in% type_intogen),1])
  
  matching_cosmic<-unique(mapping_cosmic[which(mapping_cosmic[,2] %in% type_cosmic),1])
  
  summary_vars$Match_tiss_Intogen[i]<-length(intersect(type_tiss, matching_intogen))
  summary_vars$Match_tiss_COSMIC[i]<-length(intersect(type_tiss, matching_cosmic))
  summary_vars$Matching_inboth[i]<-length(Reduce(intersect, list(matching_intogen, type_tiss, matching_cosmic)))
  summary_vars$Matching_union[i]<-length(union(intersect(matching_intogen, type_tiss), intersect(matching_cosmic, type_tiss)))
  summary_vars$Match_tiss_Intogen_name[i]<-paste(intersect(type_tiss, matching_intogen), collapse="|")
  summary_vars$Match_tiss_COSMIC_name[i]<-paste(intersect(type_tiss, matching_cosmic), collapse="|")
}

summary_vars_bytiss<-rbind(summary_vars,summary_vars)
summary_vars_bytiss$tissue<-"ALL"
summary_vars_bytiss$num_patients<-c(summary_vars[,c("Num_COSMIC")], summary_vars[,c( "Num_Intogen")])
summary_vars_bytiss$database<-c(rep("COSMIC", nrow(summary_vars)), rep("Intogen", nrow(summary_vars)))

for(i in 1:nrow(summary_vars)){
  type_tiss<-unique(unlist(strsplit(summary_vars[i, "Tiss"], " \\| ")))
  type_intogen<-na.omit(unique(unlist(strsplit(summary_vars[i, "Type_Intogen"], " \\| "))))
  type_cosmic<-na.omit(unique(unlist(strsplit(summary_vars[i, "Type_COSMIC"], " \\| "))))
  
  matching_intogen<-unique(mapping_intogen[which(mapping_intogen[,2] %in% type_intogen),1])
  matching_cosmic<-unique(mapping_cosmic[which(mapping_cosmic[,2] %in% type_cosmic),1])
  
  ###matched tissues cosmic
  if(length(intersect(type_tiss, matching_cosmic))>0){
    num_cosmic_bytiss<-c()
    for(tiss in intersect(type_tiss, matching_cosmic)){
      tiss_cosmic<-mapping_cosmic[mapping_cosmic[,1]==tiss,2]
      num_cosmic_bytiss<-c(num_cosmic_bytiss, sum(unlist(all_tiss_cosmic[[summary_vars$ID[i]]][tiss_cosmic]), na.rm=T))
    }
    num_cosmic_bytiss<-c(num_cosmic_bytiss, sum(num_cosmic_bytiss))
    names(num_cosmic_bytiss)<-c(intersect(type_tiss, matching_cosmic), "ALL_MATCHING")
    toadd<-data.frame(summary_vars[i,], num_patients=num_cosmic_bytiss)
    toadd$tissue<-names(num_cosmic_bytiss)
    toadd$database<-"COSMIC"
    summary_vars_bytiss<-rbind(summary_vars_bytiss, toadd)
  } else {
    toadd<-data.frame(summary_vars[i,], num_patients=0)
    toadd$tissue<-"ALL_MATCHING"
    toadd$database<-"COSMIC"
    summary_vars_bytiss<-rbind(summary_vars_bytiss, toadd)
  }
  
  ###matched tissues intogen
  if(length(intersect(type_tiss, matching_intogen))>0){
    num_intogen_bytiss<-c()
    for(tiss in intersect(type_tiss, matching_intogen)){
      tiss_intogen<-mapping_intogen[mapping_intogen[,1]==tiss,2]
      num_intogen_bytiss<-c(num_intogen_bytiss, sum(unlist(all_tiss_intogen[[summary_vars$ID[i]]][tiss_intogen]), na.rm=T))
    }
    num_intogen_bytiss<-c(num_intogen_bytiss, sum(num_intogen_bytiss))
    names(num_intogen_bytiss)<-c(intersect(type_tiss, matching_intogen), "ALL_MATCHING")
    toadd<-data.frame(summary_vars[i,], num_patients=num_intogen_bytiss)
    toadd$tissue<-names(num_intogen_bytiss)
    toadd$database<-"Intogen"
    summary_vars_bytiss<-rbind(summary_vars_bytiss, toadd)
  } else {
    toadd<-data.frame(summary_vars[i,], num_patients=0)
    toadd$tissue<-"ALL_MATCHING"
    toadd$database<-"Intogen"
    summary_vars_bytiss<-rbind(summary_vars_bytiss, toadd)
  }
  
}

#####save list of DAMs  with patient anno
write.xlsx(summary_vars_bytiss, file=paste(home, path_results, "/summary_DAMs_bytiss.xlsx", sep=""))

novel_patients<-data.frame(tissue=NA, num_patients=NA, database=NA, gene=NA)
for(r in 1:nrow(summary_drivers_bin)){
  gene<-rownames(summary_drivers_bin)[r]
  patients_df<-summary_vars_bytiss[summary_vars_bytiss$Gene==gene & summary_vars_bytiss$tissue %in% names(which(summary_drivers_bin[gene,]==0)),c("tissue", "num_patients","database")]
  
  if(nrow(patients_df)>0){
    patients_df$gene<-gene
    novel_patients<-rbind.data.frame(novel_patients, patients_df)
  }
}

library(ggplot2)

n_gene<-by(novel_patients, novel_patients$gene, function(x){sum(x$num_patients)})
novel_patients$gene<-factor(novel_patients$gene, levels=names(sort(n_gene, decreasing=T)))

CL_colors<-read.xlsx(paste(home, path_data, "/raw/CL_tissue_ctype_colors.xlsx", sep=""), sheet = 2)
CL_colors_v<-CL_colors$Color_hex_code
names(CL_colors_v)<-CL_colors$Cancer.Type
pdf(paste(home, path_results, "/exploration/figures/Drivers_barplot_novelpatients.pdf", sep=""),8, 5)
ggplot(novel_patients[-1,], aes(x=gene, y=num_patients, fill=tissue))+geom_bar(stat="identity")+scale_fill_manual(values=CL_colors_v)+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
