library(stringr)
library(tidyverse)
pathdata <- "data"

load('results/20250221/_totalTestedVariants.RData')
###loading input data
gene_annot <- read_csv("data/raw/gene_identifiers_20241212.csv")
### gene_identifiers_20191101 downloaded from https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20241212.csv on 20250221

cl_variants <- read_csv('data/raw/mutations_all_20241212.csv')
### mutations_all_20230202 downloaded from  https://cog.sanger.ac.uk/cmp/download/mutations_all_20241212.zip on 20250221

cl_variants <- cbind(cl_variants,gene_annot$hgnc_symbol[match(cl_variants$gene_id,gene_annot$gene_id)])
colnames(cl_variants)[ncol(cl_variants)]<-'gene_symbol_2023'


CMP_annot <- read_csv(paste(pathdata,"/raw/model_list_20241120.csv", sep="")) 
### model_list_20240110.csv downloaded from https://cog.sanger.ac.uk/cmp/download/model_list_20241120.csv on 20250129

scaled_depFC<-read.csv(paste(pathdata,'/raw/CRISPRGeneEffect.csv', sep=""), row.names = 1)
colnames(scaled_depFC)<-gsub("\\..*","",colnames(scaled_depFC))
###scaled essentiality matrices downloaded from https://depmap.org/portal/data_page/?tab=allData on 20250129 (24Q4)

toremove<-which(is.na(CMP_annot$model_id[match(rownames(scaled_depFC),CMP_annot$BROAD_ID)]))
scaled_depFC<-scaled_depFC[-toremove,]
rownames(scaled_depFC)<-CMP_annot$model_id[match(rownames(scaled_depFC),CMP_annot$BROAD_ID)]

scaled_depFC<-t(scaled_depFC)

#remove genes with missing values
scaled_depFC<-scaled_depFC[-which(rowSums(is.na(scaled_depFC))>0),]
#create a binarized matrix (essential vs non-essential)
bdep<-scaled_depFC
bdep[scaled_depFC<=(-0.5)]<-1
bdep[scaled_depFC>(-0.5)]<-0

inTOgen_drivers<-read.table("data/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", sep='\t',stringsAsFactors = FALSE,header=TRUE)
load('results/20250221/_allHits.RData')
load('results/20250221/_allDAM_bearing_genes.RData')
load('results/20250221/_allDAMs.RData')
driver_genes<-inTOgen_drivers$SYMBOL

mapping<-read.csv("data/raw/intOGen ctype mapping_AS.csv",header = TRUE,row.names = 1, sep=";")

act_driver_genes<-inTOgen_drivers$SYMBOL[inTOgen_drivers$ROLE=="Act"]

####cell lines with co-occurrent DAM and essential driver mutation or only with DAM and not an essential driver mutation
co_occurrent<-c()
non_co_occurrent<-c()
ct_co_occurrent<-c()
ct_non_co_occurrent<-c()
DAM_bearing_lines<-unique(unlist(strsplit(allHits$ps_cl[-which(allHits$GENE %in% driver_genes)], ", ") ))
for(l in DAM_bearing_lines){
  ess_genes<-names(which(bdep[,l]==1))
  ct<-CMP_annot$cancer_type[CMP_annot$model_id==l]
  ct_into<-strsplit(mapping[ct, 1], " \\| ")
  act_genes<-inTOgen_drivers$SYMBOL[inTOgen_drivers$ROLE=="Act" & inTOgen_drivers$CANCER_TYPE %in% ct_into]
  
  other_Mutated_Genes<-cl_variants$gene_symbol_2023[which(cl_variants$model_id==l)]
  other_Mutated_Act_Genes<-intersect(other_Mutated_Genes, act_genes)
  other_Mutated_Act_Ess_Genes<-intersect(other_Mutated_Act_Genes, ess_genes)
  
  if(length(other_Mutated_Act_Ess_Genes)>0){
    co_occurrent<-c(co_occurrent, l)
    ct_co_occurrent<-c(ct_co_occurrent, ct)
  } else {
    non_co_occurrent<-c(non_co_occurrent, l)
    ct_non_co_occurrent<-c(ct_non_co_occurrent, ct)
  }
}

print(paste(length(co_occurrent)/(length(co_occurrent)+length(non_co_occurrent))*100, "% of the cell lines with a DAM have a co-occurrent essential known driver mutated"))


df_lines<-data.frame(line_index=c(co_occurrent, non_co_occurrent), tiss=c(ct_co_occurrent, ct_non_co_occurrent), 
                     type=c(rep("co-occurrent", rep(length(co_occurrent))), rep("non-co-occurrent", rep(length(non_co_occurrent)))))
df_lines$tiss<-factor(df_lines$tiss, levels=c(names(sort(table(df_lines$tiss), decreasing=T))))


pdf("results/20250221/co_occurrence_lines.pdf", 10, 7)
ggplot(df_lines, aes(x=tiss, fill=type))+geom_bar()+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

df_tiss<-data.frame(table(ct_non_co_occurrent))
rownames(df_tiss)<-df_tiss[,1]
df_tiss$cooccurr<-0
df_tiss[names(table(ct_co_occurrent)), "cooccurr"]<-c(table(ct_co_occurrent))
colnames(df_tiss)<-c("tissue", "non_cooccurr", "cooccurr")
df_tiss$tissue<-rownames(df_tiss)
df_tiss$non_cooccurr[is.na(df_tiss$non_cooccurr)]<-0
df_tiss$ratio<-df_tiss$non_cooccurr/(df_tiss$cooccurr+df_tiss$non_cooccurr)

ggplot(df_tiss, aes(x=tissue, y=ratio))+geom_bar(stat="identity")+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(paste(sum(df_tiss$ratio==1), "number of tissues with DAMs exclusively in cell lines without a co-occurrent mutated essential driver"))

pdf("results/20250221/co_occurrence_lines_hist.pdf", 7, 5)
ggplot(df_tiss, aes(x=ratio))+geom_histogram()+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

write.csv(df_tiss, file="results/20250221/co_occurrence_lines_pertiss.csv", quote = F, row.names = F)

####DAMs co-occurrent or not co-occurrent with an essential driver mutation (separately for each cell line with the DAM, so each DAM can be computed multiple times)
co_occurrent_d<-c()
non_co_occurrent_d<-c()
ct_co_occurrent_d<-c()
ct_non_co_occurrent_d<-c()

for(i in 1:nrow(allDAMs)){
  ct<-allDAMs$ctype[i]
  ind<-which(cl_variants$gene_symbol_2023==allDAMs$GENE[i] & cl_variants$protein_mutation==allDAMs$var[i] & cl_variants$model_id %in% CMP_annot$model_id[CMP_annot$cancer_type==ct])
  cl<-cl_variants$model_id[ind]
  cl<-intersect(cl, colnames(bdep))
  
  for(l in cl){
    ess_genes<-names(which(bdep[,l]==1))
    ct_into<-strsplit(mapping[ct, 1], " \\| ")
    act_genes<-inTOgen_drivers$SYMBOL[inTOgen_drivers$ROLE=="Act" & inTOgen_drivers$CANCER_TYPE %in% ct_into]
    
    other_Mutated_Genes<-cl_variants$gene_symbol_2023[which(cl_variants$model_id==l)]
    other_Mutated_Act_Genes<-intersect(other_Mutated_Genes, act_genes)
    other_Mutated_Act_Ess_Genes<-intersect(other_Mutated_Act_Genes, ess_genes)
    
    if(length(other_Mutated_Act_Ess_Genes)>0){
      co_occurrent_d<-c(co_occurrent_d, i)
      ct_co_occurrent_d<-c(ct_co_occurrent_d, ct)
    } else {
      non_co_occurrent_d<-c(non_co_occurrent_d, i)
      ct_non_co_occurrent_d<-c(ct_non_co_occurrent_d, ct)
    }
    
  }
  }

length(co_occurrent_d)/(length(co_occurrent_d)+length(non_co_occurrent_d))

df_lines<-data.frame(line_index=c(co_occurrent_d, non_co_occurrent_d), tiss=c(ct_co_occurrent_d, ct_non_co_occurrent_d), 
                     type=c(rep("co-occurrent", rep(length(co_occurrent_d))), rep("non-co-occurrent", rep(length(non_co_occurrent_d)))))
df_lines$tiss<-factor(df_lines$tiss, levels=c(names(sort(table(df_lines$tiss), decreasing=T))))

pdf("results/20250221/co_occurrence_inlines_DAMs.pdf", 10, 7)
ggplot(df_lines, aes(x=tiss, fill=type))+geom_bar()+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

df_tiss_d<-data.frame(table(ct_non_co_occurrent_d))
rownames(df_tiss_d)<-df_tiss_d[,1]
df_tiss_d$cooccurr<-0
df_tiss_d[names(table(ct_co_occurrent_d)), "cooccurr"]<-c(table(ct_co_occurrent_d))
colnames(df_tiss_d)<-c("tissue", "non_cooccurr", "cooccurr")
df_tiss_d$tissue<-rownames(df_tiss_d)
df_tiss_d$non_cooccurr[is.na(df_tiss_d$non_cooccurr)]<-0
df_tiss_d$ratio<-df_tiss_d$non_cooccurr/(df_tiss_d$cooccurr+df_tiss_d$non_cooccurr)

ggplot(df_tiss_d, aes(x=tissue, y=ratio))+geom_bar(stat="identity")+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(paste(sum(df_tiss_d$ratio==1), "number of tissues with DAMs exclusively in cell lines without a co-occurrent mutated essential driver"))

pdf("results/20250221/co_occurrence_inlines_DAMs_hist.pdf", 7, 5)
ggplot(df_tiss_d, aes(x=ratio))+geom_histogram()+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

write.csv(df_tiss_d, file="results/20250221/co_occurrence_inlines_DAMs_pertiss.csv", quote = F, row.names = F)

