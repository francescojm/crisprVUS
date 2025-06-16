### reverted main to aurora's version
.libPaths(c("/home/aurora.savino/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))
library(CELLector)
library(tidyverse)

####set path
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20220208'

###load input data
gene_annot <- read_csv(paste(pathdata, "/gene_identifiers_20191101.csv", sep=""))

load(paste(pathdata,"/R/Sanger_Broad_higQ_scaled_depFC.RData", sep=""))
load(paste(pathdata,'/R/Sanger_Broad_higQ_bdep.RData', sep=""))

load(paste(pathdata,'/preComputed/ADaM.RData', sep=""))
load(paste(pathdata,'/preComputed/FiPer_outputs.RData', sep=""))

## latest sanger/broad unreleased yet variants hg38 
cl_variants <- read_csv(paste(pathdata,'/mutations_all_latest.csv', sep=""))
cl_variants <- cbind(cl_variants,gene_annot$hgnc_symbol[match(cl_variants$gene_id,gene_annot$gene_id)])
colnames(cl_variants)[ncol(cl_variants)]<-'gene_symbol_2019'

CMP_annot <- read_csv(paste(pathdata,"/model_list_20210611.csv", sep="")) # from https://cog.sanger.ac.uk/cmp/download/model_list_20210611.csv

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

#OPTIONAL: call the code in a bash script and use input_iter to select the tissue to analyze
input_iter<-as.numeric(commandArgs(trailingOnly = TRUE)[2])
print(paste("executing iter", input_iter))

#run the main pipeline for each tissue separately
for (ctiss in tissues[input_iter]){
  
  print(ctiss)
  oldEnv<-ls()
  oldEnv<-ls()

  source(paste(pathscript,'/01_mut_CRISPR.R', sep=""))
  newEnv<-ls()
  rm(list=setdiff(newEnv,oldEnv))

  source(paste(pathscript,'/02_add_driver_infos.R', sep=""))
  newEnv<-ls()
  rm(list=setdiff(newEnv,oldEnv))

  source(paste(pathscript, '/03_add_GDSC_validation.R', sep=""))
   newEnv<-ls()
   rm(list=setdiff(newEnv,oldEnv))

}

