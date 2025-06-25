library(stringr)
library(tidyverse)

load('../results/_VUS2024build/20241002/_totalTestedVariants.RData')
###loading input data
gene_annot <- read_csv("../data/_VUS2024build/raw/gene_identifiers_20191101.csv")

cl_variants <- read_csv('../data/_VUS2024build/raw/mutations_all_20230202.csv')
### mutations_all_20230202 downloaded from  https://cog.sanger.ac.uk/cmp/download/mutations_all_20230202.zip on 20241002
cl_variants <- cbind(cl_variants,gene_annot$hgnc_symbol[match(cl_variants$gene_id,gene_annot$gene_id)])
colnames(cl_variants)[ncol(cl_variants)]<-'gene_symbol_2023'

inTOgen_drivers<-read.table("../data/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv", sep='\t',stringsAsFactors = FALSE,header=TRUE)
load('../results/_VUS2024build/20241002/_allHits.RData')
load('../results/_VUS2024build/20241002/_allDAM_bearing_genes.RData')

inTOgen_ctype_mapping<-read.table('../data/_VUS2024build/raw/intOGen ctype mapping.tsv',sep='\t',header=TRUE,row.names = 1)

DAMs_comutatedActMutations_likelyFunctional<-lapply(1:nrow(allHits),function(x){
  print(x)
  ps_cl<-unlist(str_split(allHits$ps_cl[x],', '))
  mappedCtype<-setdiff(unlist(str_split(inTOgen_ctype_mapping[allHits$ctype[x],1],' | ')),'|')
  
  otherActivatingDrivers<-inTOgen_drivers$SYMBOL[which(is.element(inTOgen_drivers$CANCER_TYPE,mappedCtype) & 
                                                         is.element(inTOgen_drivers$ROLE,c('Act','ambigous')))]
  
  if(length(ps_cl)>1){
    otherMutatedGenes<-Reduce('intersect',lapply(ps_cl,function(cl){
      oMG<-totalTestedVariants$gene_symbol[which(totalTestedVariants$positiveCls==cl)]
    }))  
  }else{
    otherMutatedGenes<-totalTestedVariants$gene_symbol[which(totalTestedVariants$positiveCls==ps_cl)]
    }

  otherMutatedActDrivers<-intersect(otherMutatedGenes,otherActivatingDrivers)
})

DAMs_comutatedActMutations_all<-lapply(1:nrow(allHits),function(x){
  print(x)
  ps_cl<-unlist(str_split(allHits$ps_cl[x],', '))
  mappedCtype<-setdiff(unlist(str_split(inTOgen_ctype_mapping[allHits$ctype[x],1],' | ')),'|')
  
  otherActivatingDrivers<-inTOgen_drivers$SYMBOL[which(is.element(inTOgen_drivers$CANCER_TYPE,mappedCtype) & 
                                                         is.element(inTOgen_drivers$ROLE,c('Act','ambigous')))]
  
  if(length(ps_cl)>1){
    otherMutatedGenes<-Reduce('intersect',lapply(ps_cl,function(cl){
      oMG<-cl_variants$gene_symbol_2023[which(cl_variants$model_name==cl)]
    }))  
  }else{
    otherMutatedGenes<-cl_variants$gene_symbol_2023[which(cl_variants$model_name==ps_cl)]
  }
  
  otherMutatedActDrivers<-intersect(otherMutatedGenes,otherActivatingDrivers)
})

paste(round(100*
        sum(unlist(lapply(DAMs_comutatedActMutations_likelyFunctional,'length'))>0)/length(DAMs_comutatedActMutations_likelyFunctional),2),
      "% of hits are co-occurrent with a likely functional mutation in a cancer-type matching Cancer Driver gene with intOGen role = Activating or Ambigous in that cancer-type")

sum(unlist(lapply(DAMs_comutatedActMutations_likelyFunctional,'length'))>0)

paste(round(100*
              sum(unlist(lapply(DAMs_comutatedActMutations_all,'length'))>0)/length(DAMs_comutatedActMutations_likelyFunctional),2),
      "% of hits are co-occurrent with a mutation (any kind) in a cancer-type matching Cancer Driver gene with intOGen role = Activating or Ambigous in that cancer-type")

sum(unlist(lapply(DAMs_comutatedActMutations_all,'length'))>0)

DAMs_comutatedActMutations_likelyFunctional_also_DAMbearingG<-unlist(
lapply(1:length(DAMs_comutatedActMutations_likelyFunctional),function(x){
  print(x)
  currentCtype<-allHits$ctype[x]
  currentDAM_bearing_genes<-allDAM_bearing_genes$allDAM_bearing_genes[grep(currentCtype,allDAM_bearing_genes$DAMbearing_in)]
  comutatedG<-unlist(DAMs_comutatedActMutations_likelyFunctional[x])
  length(intersect(comutatedG,currentDAM_bearing_genes)>0)
}))

summary(as.factor(DAMs_comutatedActMutations_likelyFunctional_also_DAMbearingG))

DAMs_comutatedActMutations_all_also_DAMbearingG<-unlist(
  lapply(1:length(DAMs_comutatedActMutations_all),function(x){
    print(x)
    currentCtype<-allHits$ctype[x]
    currentDAM_bearing_genes<-allDAM_bearing_genes$allDAM_bearing_genes[grep(currentCtype,allDAM_bearing_genes$DAMbearing_in)]
    comutatedG<-unlist(DAMs_comutatedActMutations_all[x])
    length(intersect(comutatedG,currentDAM_bearing_genes)>0)
  }))

summary(as.factor(DAMs_comutatedActMutations_all_also_DAMbearingG))
