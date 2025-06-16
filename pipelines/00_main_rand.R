#same code as 00_main but calling 01_mut_CRISPR_rand.R, reshuffling cell lines' labels in dependency data
.libPaths(c("/home/aurora.savino/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))
print(.libPaths())
library(CELLector)
library(tidyverse)

pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20220208'

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


print(paste(nrow(CMP_annot),'annotated models in the Cell Models Passports'))

print(paste('of which',length(which(is.element(CMP_annot$model_id,cl_variants$model_id))),'with mutation data'))

CMP_annot<-CMP_annot[which(is.element(CMP_annot$model_name,colnames(scaled_depFC))),]
print(paste('of which',nrow(CMP_annot),'with high quality CRISPR data'))

tissues<-CMP_annot$cancer_type

st<-summary(as.factor(tissues))

tissues<-sort(setdiff(tissues,names(which(st<5))))
tissues<-setdiff(tissues,c('Other Solid Carcinomas','Other Solid Cancers','Other Sarcomas'))

CMP_annot<-CMP_annot[which(is.element(CMP_annot$cancer_type,tissues)),]

incl_cl_annot<-CMP_annot

save(incl_cl_annot,file=paste(resultPath,'/_incl_cl_annot.RData',sep=''))

tissues<-sort(tissues)


library(parallel)
library(sets)
library(binaryLogic)
library(tidyverse)
library(CoRe)

input_iter<-as.numeric(commandArgs(trailingOnly = TRUE)[2])
print(paste("executing iter", input_iter))
for (ctiss in tissues[input_iter]){
  RESTOT_rand<-vector(mode="list", length=1000)
  print(ctiss)
set.seed(62587)
trials<-c(1:1000)
numCores <- detectCores()
fun_iter<-function(x){
  print(x)
  source(paste(pathscript,'/01_mut_CRISPR_rand.R', sep=""), verbose=F)
}
RESTOT_rand <- mclapply(trials, fun_iter, mc.cores = numCores)
save(RESTOT_rand,file=paste(resultPath,'/',ctiss,'_results_rand.RData',sep=''))
}

