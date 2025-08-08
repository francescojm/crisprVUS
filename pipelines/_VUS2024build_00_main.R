set.seed(123)
library(CELLector)
library(tidyverse)

####setting paths
pathdata <- "data"
pathscript <- "../crisprVUS/pipelines"
resultPath<-'results/20250807_RR_th1.71/'

###loading input data
gene_annot <- read_csv(paste(pathdata, "/raw/gene_identifiers_20241212.csv", sep=""))
### gene_identifiers_20191101 downloaded from https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20241212.csv on 20250221

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

load(paste(pathdata,'/Robj/ADaM.RData', sep=""))
load(paste(pathdata,'/Robj/FiPer_outputs.RData', sep=""))
### Rbjects precomputed as in Vinceti et al, BMC Genomics, 2021

cl_variants <- read_csv(paste(pathdata,'/raw/mutations_all_20241212.csv', sep=""))
### mutations_all_20230202 downloaded from  https://cog.sanger.ac.uk/cmp/download/mutations_all_20241212.zip on 20250221

cl_variants <- cbind(cl_variants,gene_annot$hgnc_symbol[match(cl_variants$gene_id,gene_annot$gene_id)])
colnames(cl_variants)[ncol(cl_variants)]<-'gene_symbol_2023'


###print a few figures
print(paste(nrow(CMP_annot),'annotated models in the Cell Models Passports'))
print(paste('of which',length(which(is.element(CMP_annot$model_id,cl_variants$model_id))),'with mutation data'))

#select only cell lines with data in both CRISPR screens and CMP annotation file
CMP_annot<-CMP_annot[which(is.element(CMP_annot$model_id,colnames(scaled_depFC))),]
print(paste('of which',nrow(CMP_annot),'with high quality CRISPR data'))

#select tissues
tissues<-CMP_annot$cancer_type
st<-summary(as.factor(tissues))

tissues<-sort(setdiff(tissues,names(which(st<5))))
tissues<-setdiff(tissues,c('Other Solid Carcinomas','Other Solid Cancers','Other Sarcomas', "Other Blood Cancers", "Non-Cancerous"))

CMP_annot<-CMP_annot[which(is.element(CMP_annot$cancer_type,tissues)),]

incl_cl_annot<-CMP_annot
save(incl_cl_annot,file=paste(resultPath,'/_incl_cl_annot.RData',sep=''))
write.table(incl_cl_annot,sep="\t",quote=FALSE,file=paste(resultPath,'_incl_cl_annot.tsv',sep=''))

tissues<-sort(tissues)

load(paste(pathdata, "/Robj/basal_exp.RData", sep=""))


#run the main pipeline for each tissue separately
for (ctiss in tissues){

  print(ctiss)
  oldEnv<-ls()
  oldEnv<-ls()

  source(paste(pathscript,'/_VUS2024build_01_mut_CRISPR.R', sep=""))
  newEnv<-ls()
  rm(list=setdiff(newEnv,oldEnv))

  source(paste(pathscript,'/_VUS2024build_02_add_driver_infos.R', sep=""))
  newEnv<-ls()
  rm(list=setdiff(newEnv,oldEnv))

  source(paste(pathscript, '/_VUS2024build_03_add_GDSC_validation.R', sep=""))
  newEnv<-ls()
  rm(list=setdiff(newEnv,oldEnv))

}
