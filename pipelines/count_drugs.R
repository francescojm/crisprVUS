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
  
targets<-c()
drugs<-c()
ntests<-0
targets_sel<-c()
drugs_sel<-c()
ntests_sel<-0

for (ctiss in tissues){
drugTargetInfo <- read.table(paste(pathdata,'/raw/drug-target_data_hgvs_clean_Goncalves_et_all.txt', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
gdsc1<-read.csv(paste(pathdata,'/raw/GDSC1_fitted_dose_response_25Feb20.csv', sep=""),header = TRUE,stringsAsFactors = FALSE)
gdsc2<-read.csv(paste(pathdata,'/raw/GDSC2_fitted_dose_response_25Feb20.csv', sep=""),header = TRUE,stringsAsFactors = FALSE)
colnames(gdsc1)[1]<-"DATASET"
colnames(gdsc2)[1]<-"DATASET"

gdscAll<-rbind(gdsc1,gdsc2)

clTiss<-CMP_annot$model_name[CMP_annot$cancer_type==ctiss]

  load(paste(resultPath, "/", ctiss, "_results_ext.RData", sep=""))



  for(x in 1:nrow(RESTOT)){
    
    target<-RESTOT$GENE[x]
    variant<-RESTOT$var[x]
    cellLines<-unlist(str_split(RESTOT$ps_cl[x],', '))
    ids<-which(drugTargetInfo$Gene.Target==target)
     
    drug_ids<-drugTargetInfo$Drug.ID[ids]
    drug_names<-drugTargetInfo$Name[ids]
  
    if (length(ids)>0){
        
      print(x)
      data1<-gdsc1[which(is.element(gdsc1$DRUG_ID,drug_ids) & is.element(gdsc1$CELL_LINE_NAME,cellLines)),]
      data2<-gdsc2[which(is.element(gdsc2$DRUG_ID,drug_ids) & is.element(gdsc2$CELL_LINE_NAME,cellLines)),]
      
      data1<-rbind(data1,data2)
      targets<-c(targets, target)
      drugs<-c(drugs, data1$DRUG_NAME)
      ntests<-ntests+nrow(data1)
      
      if(RESTOT$rank_ratio[x]<1.6 & RESTOT$medFitEff[x]< -0.5 & RESTOT$pval_rand[x]<0.2){
      
      targets_sel<-c(targets_sel, target)
      drugs_sel<-c(drugs_sel, data1$DRUG_NAME)
      ntests_sel<-ntests_sel+nrow(data1)
      }
       
    }
  }
}
