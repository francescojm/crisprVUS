library(stringr)
library(tidyverse)
library(openxlsx)

####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250808_bugFixed_and_RR_th.1.71_wr/'

clc<-read.xlsx(paste(pathdata,'/raw/CL_tissue_ctype_colors.xlsx',sep=''),sheet = 2,rowNames = TRUE)

tract<-read.table(paste(pathdata,'/raw/Tractability_pacini_et_al_2024.txt',sep=''),sep="\t",
                  header=TRUE,row.names = NULL)

tractable_targets<-tract$id[which(tract$min_bucket==1)]

load(paste(resultPath,'_totalTestedVariants.RData',sep=''))

###loading input data
gene_annot <- read_csv(paste(pathdata,"/raw/gene_identifiers_20241212.csv",sep=''))
### gene_identifiers_20191101 downloaded from https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20241212.csv on 20250221

cl_variants <- read_csv(paste(pathdata,'/raw/mutations_all_20241212.csv',sep=''))
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

inTOgen_drivers<-read.table(paste(pathdata,"/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", sep=''),sep='\t',stringsAsFactors = FALSE,header=TRUE)
load(paste(resultPath,'/_allHits.RData',sep=''))
load(paste(resultPath,'/_allDAM_bearing_genes.RData',sep=''))
load(paste(resultPath,'/_allDAMs.RData',sep=''))
driver_genes<-inTOgen_drivers$SYMBOL

mapping<-read.csv(paste(pathdata,"/raw/intOGen ctype mapping_AS.csv",sep=''),header = TRUE,row.names = 1, sep=";")

act_driver_genes<-inTOgen_drivers$SYMBOL[inTOgen_drivers$ROLE=="Act"]

###########################################
#### Line-level analysis
############################################

####cell lines with co-occurrent DAM and essential driver mutation or only with DAM and not an essential driver mutation
co_occurrent<-c()
non_co_occurrent<-c()
ct_co_occurrent<-c()
ct_non_co_occurrent<-c()

####selecting only cell lines with DAMs in an unreported DAMbgs
DAM_bearing_lines<-unique(unlist(strsplit(allDAMs$ps_cl[-which(allDAMs$GENE %in% driver_genes)], ", ") ))
known_DAM_bearing_lines<-unique(unlist(strsplit(allDAMs$ps_cl[which(allDAMs$GENE %in% driver_genes)], ", ")))


RES<-do.call(rbind,lapply(DAM_bearing_lines,function(l){
  print(l)
  ess_genes<-names(which(bdep[,l]==1))
  
  ct<-CMP_annot$cancer_type[CMP_annot$model_id==l]
  
  ct_into<-unlist(strsplit(mapping[ct, 1], " \\| "))
  
  act_genes<-inTOgen_drivers$SYMBOL[which(inTOgen_drivers$ROLE=="Act" & is.element(inTOgen_drivers$CANCER_TYPE,ct_into))]
  
  all_Mutated_Genes<-cl_variants$gene_symbol_2023[which(cl_variants$model_id==l)]
  
  all_Mutated_act_Genes<-intersect(all_Mutated_Genes, act_genes)
  
  Act_Ess_Drivers<-intersect(all_Mutated_act_Genes, ess_genes)
  
  idcell<-grep(l,allDAMs$ps_cl)
  gg<-allDAMs$GENE[idcell]
  vars<-allDAMs$var[idcell]
  
  observed_DAMs<-paste(gg,vars,sep='')
  
  DAMs_in_unreported_DAMbgs<-paste(observed_DAMs[which(!is.element(allDAMs$GENE[idcell],driver_genes))],collapse=', ')
  
  druggable_DAMs_in_unreported_DAMbgs<-paste(intersect(setdiff(gg,driver_genes),tractable_targets),collapse=', ')
  
  DAMs_in_known_DAMbgs<-paste(observed_DAMs[which(is.element(allDAMs$GENE[idcell],act_genes))],collapse=', ')
  
  druggable_known_DAMs<-paste(intersect(intersect(gg,act_genes),tractable_targets),collapse=', ')
  
  druggable_Act_Ess_Drivers<-paste(intersect(Act_Ess_Drivers,tractable_targets),collapse=', ')
  
  Act_Ess_Drivers<-paste(Act_Ess_Drivers,collapse=', ')
  
  c(ct,l,DAMs_in_unreported_DAMbgs,druggable_DAMs_in_unreported_DAMbgs,DAMs_in_known_DAMbgs,druggable_known_DAMs,Act_Ess_Drivers,druggable_Act_Ess_Drivers)
}))

colnames(RES)<-c('cancer_type','cell line','DAMs in unreported DAMbgs','of which druggable','DAMs in known DAMbgs','of which druggable','GoF mutated and essential drivers','of which druggable')

DAMcanONC_co_occ<-RES
save(DAMcanONC_co_occ,file=paste(resultPath,'_DAM_canonical_oncogenic_addition_co_occurrence.RData',sep=''))

write.table(RES, file=paste(resultPath,'_DAM_canonical_oncogenic_addition_co_occurrence.tsv',sep=''),
          quote = F, row.names = F,sep='\t')

nn<-100*length(which(RES[,5]==''))/nrow(RES)

print(paste(nn,'% of CCLs from a given cancer-type, with an DAM in an unreported DAMbgs, lacks DAMs in cancer-type specific known GoF driver genes',sep=''))

nn<-100*length(which(RES[,7]==''))/nrow(RES)

print(paste(nn,'% of CCLs from a given cancer-type, with an DAM in an unreported DAMbgs, lacks cancer-type specific known GoF driver genes that are mutated and essential',sep=''))

occs<-unlist(lapply(RES[,5],function(x){x==''}))+0

utiss<-unique(RES[,1])

coc_results<-do.call(rbind,lapply(utiss, function(x){
  temp<-occs[which(RES[,1]==x)]
  n<-length(temp)
  s<-sum(temp)
  o<-n-s
  c(o,s)
}))

rownames(coc_results)<-utiss
colnames(coc_results)<-c('co-occurring oncAdd','other')
print(paste('median across cancer types =', median(100*coc_results[,2]/rowSums(coc_results))))

print(sort(100*coc_results[,2]/rowSums(coc_results)))

coc_results<-coc_results[order(rowSums(coc_results),decreasing=TRUE),]

pdf(paste(resultPath,'_figures_source/unreportedDAMs_oncogenicAddiction_coOcc.pdf',sep=''),11,7)
par(mar=c(16,4,2,2))
barplot(t(coc_results),border=FALSE,las=2,main=paste(sum(c(coc_results)),'cell lines with DAM in unreported DAMbgs'),ylab='n. cell lines',
        col=c('darkgray','blue'))

legend('topright',c('co-occurrent oncogenic addiction','others'),fill=c('darkgray','blue'),border=FALSE)
dev.off()



occs<-unlist(lapply(RES[,7],function(x){x==''}))+0

coc_results<-do.call(rbind,lapply(utiss, function(x){
  temp<-occs[which(RES[,1]==x)]
  n<-length(temp)
  s<-sum(temp)
  o<-n-s
  c(o,s)
}))

rownames(coc_results)<-utiss
colnames(coc_results)<-c('co-occurring oncAdd','other')
print(paste('median across cancer types =', median(100*coc_results[,2]/rowSums(coc_results))))

print(sort(100*coc_results[,2]/rowSums(coc_results)))

coc_results<-coc_results[order(rowSums(coc_results),decreasing=TRUE),]

pdf(paste(resultPath,'_figures_source/unreportedDAMs_oncogenicAddiction_coOcc_relaxed_criterion.pdf',sep=''),11,7)
par(mar=c(16,4,2,2))
barplot(t(coc_results),border=FALSE,las=2,main=paste(sum(c(coc_results)),'cell lines with DAM in unreported DAMbgs'),ylab='n. cell lines',
        col=c('darkgray','blue'))

legend('topright',c('co-occurrent oncogenic addiction','others'),fill=c('darkgray','blue'),border=FALSE)
dev.off()


