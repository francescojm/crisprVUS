.libPaths(c("/home/aurora.savino/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))
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


tissues<-sort(tissues)

variantSpectrum<-function(cl_var,gene){
  
  tmp<-cl_var[cl_var$gene_symbol_2019==gene,c('protein_mutation','gene_symbol_2019')]
  
  aa<-sort(unique(tmp$protein_mutation))
  aa<-setdiff(aa, c("-", "p.?"))
  
  aa<-tmp[match(aa,tmp$protein_mutation),]
  
  return(aa)
}

tot_genes_comb<-0
tot_vars_comb<-0
tot_genes_comb_list<-c()
tot_vars_comb_list<-c()
for (ctiss in tissues){
  
  ts_depFC<-scaled_depFC[,CMP_annot$model_name[CMP_annot$cancer_type==ctiss]]
ts_bdep<-bdep[,CMP_annot$model_name[CMP_annot$cancer_type==ctiss]]

ts_cl_variants<-cl_variants[which(is.element(cl_variants$model_name,colnames(ts_depFC))),]
ts_cl_variants<-ts_cl_variants[which(is.element(ts_cl_variants$gene_symbol_2019,rownames(ts_depFC))),]

genesToTest<-unique(ts_cl_variants$gene_symbol_2019)

vs_spec_cardinality<-unlist(lapply(lapply(1:length(genesToTest),
                                          function(x){
                                            dd<-variantSpectrum(cl_var = ts_cl_variants,gene = genesToTest[x])
                                            dd<-dd$protein_mutation}),'length'))

names(vs_spec_cardinality)<-genesToTest

vs_spec<-lapply(1:length(genesToTest),
                                          function(x){
                                            dd<-variantSpectrum(cl_var = ts_cl_variants,gene = genesToTest[x])
                                            dd<-dd$protein_mutation})
names(vs_spec)<-genesToTest


##### questo esclude gene mutatissimi perche' enormi ##### tipo TTN e anche tanti TSGs
genesToTest<-sort(names(which(vs_spec_cardinality<=10 & vs_spec_cardinality>0)))

##### questo esclude gene pan-cancer cf and common-essential genes
genesToTest<-setdiff(genesToTest,ADaM)
genesToTest<-setdiff(genesToTest,Perc_AUC)

vars<-c()
for(i in 1:length(vs_spec[genesToTest])){
  vars<-c(vars, paste(rep(names(vs_spec[genesToTest])[i], length(unlist(vs_spec[genesToTest][i]))), unlist(vs_spec[genesToTest][i])))
}

tot_genes_comb<-tot_genes_comb+length(genesToTest)
tot_vars_comb<-tot_vars_comb+sum(vs_spec_cardinality[genesToTest])
tot_genes_comb_list<-c(tot_genes_comb_list, genesToTest)
tot_vars_comb_list<-c(tot_vars_comb_list,vars)
}

