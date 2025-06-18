set.seed(123)
library(CELLector)
library(tidyverse)

####setting paths
pathdata <- "../data/_VUS2024build/"
pathscript <- "pipelines"
resultPath<-'../results/_VUS2024build/20250617paperTest/'

###loading input data
gene_annot <- read_csv(paste(pathdata, "/raw/gene_identifiers_20191101.csv", sep=""))
### gene_identifiers_20191101 downloaded from https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20191101.csv on 20241002

load(paste(pathdata,"/Robj/Sanger_Broad_higQ_scaled_depFC.RData", sep=""))
load(paste(pathdata,'/Robj/Sanger_Broad_higQ_bdep.RData', sep=""))
### Rbjects built from the essentiality matrices downloaded from https://cog.sanger.ac.uk/cmp/download/Project_score_combined_Sanger_v1_Broad_20Q2_20210311.zip on 20241002

load(paste(pathdata,'/Robj/ADaM.RData', sep=""))
load(paste(pathdata,'/Robj/FiPer_outputs.RData', sep=""))
### Rbjects precomputed as in Vinceti et al, BMC Genomics, 2021

cl_variants <- read_csv(paste(pathdata,'/raw/mutations_all_20230202.csv', sep=""))
### mutations_all_20230202 downloaded from  https://cog.sanger.ac.uk/cmp/download/mutations_all_20230202.zip on 20241002

cl_variants <- cbind(cl_variants,gene_annot$hgnc_symbol[match(cl_variants$gene_id,gene_annot$gene_id)])
colnames(cl_variants)[ncol(cl_variants)]<-'gene_symbol_2023'

CMP_annot <- read_csv(paste(pathdata,"/raw/model_list_20240110.csv", sep="")) 
### model_list_20240110.csv downloaded from https://cog.sanger.ac.uk/cmp/download/model_list_20240110.csv on 20241002

###print a few figures
print(paste(nrow(CMP_annot),'annotated models in the Cell Models Passports'))
print(paste('of which',length(which(is.element(CMP_annot$model_id,cl_variants$model_id))),'with mutation data'))

#select only cell lines with data in both CRISPR screens and CMP annotation file
CMP_annot<-CMP_annot[which(is.element(CMP_annot$model_name,colnames(scaled_depFC))),]
print(paste('of which',nrow(CMP_annot),'with high quality CRISPR data'))

#select tissues
tissues<-CMP_annot$cancer_type
st<-summary(as.factor(tissues))

tissues<-sort(setdiff(tissues,names(which(st<5))))
tissues<-setdiff(tissues,c('Other Solid Carcinomas','Other Solid Cancers','Other Sarcomas'))

CMP_annot<-CMP_annot[which(is.element(CMP_annot$cancer_type,tissues)),]

incl_cl_annot<-CMP_annot
save(incl_cl_annot,file=paste(resultPath,'/_incl_cl_annot.RData',sep=''))

tissues<-sort(tissues)

load(paste(pathdata,'Robj/basa_exp_FPKM.RData',sep=''))


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

source(paste(pathscript, '/_VUS2024build_04_First_analysisSummary_andDriverEnrichments.R', sep=""))

