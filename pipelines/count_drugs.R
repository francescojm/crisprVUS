library(tidyverse)

pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250221'

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

tissues<-sort(tissues)

variantSpectrum<-function(cl_var,gene){
  
  tmp<-cl_var[cl_var$gene_symbol_2023==gene,c('protein_mutation','gene_symbol_2023')]
  
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
  ### drug-target_data_hgvs_clean_Goncalves_et_all.txt built from https://www.embopress.org/doi/suppl/10.15252/msb.20199405/suppl_file/msb199405-sup-0003-datasetev2.xlsx on 20241003
  
  
  gdsc1<-read.csv(paste(pathdata,'/raw/GDSC1_fitted_dose_response_27Oct23.csv', sep=""),header = TRUE,stringsAsFactors = FALSE)
  gdsc2<-read.csv(paste(pathdata,'/raw/GDSC2_fitted_dose_response_27Oct23.csv', sep=""),header = TRUE,stringsAsFactors = FALSE)
  ### GDSC1_fitted_dose_response_27Oct23.csv and GDSC2_fitted_dose_response_27Oct23.csv have been downloaded from: 
  ### https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC1_fitted_dose_response_27Oct23.xlsx and https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx
  ### on the 20241003
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
      data1<-gdsc1[which(is.element(gdsc1$DRUG_ID,drug_ids) & is.element(gdsc1$SANGER_MODEL_ID,cellLines)),]
      data2<-gdsc2[which(is.element(gdsc2$DRUG_ID,drug_ids) & is.element(gdsc2$SANGER_MODEL_ID,cellLines)),]
      
      data1<-rbind(data1,data2)
      targets<-c(targets, target)
      drugs<-c(drugs, data1$DRUG_NAME)
      ntests<-ntests+nrow(data1)
      
      if(RESTOT$rank_ratio[x]<1.6 & RESTOT$medFitEff[x]< -0.5 & RESTOT$pval_rand[x]<0.2){
      
        if(nrow(data1)>0){
          targets_sel<-c(targets_sel, target)
          drugs_sel<-c(drugs_sel, data1$DRUG_NAME)
          ntests_sel<-ntests_sel+nrow(data1) 
        }
      }
       
    }
  }
}

print(paste("Number of tests:", ntests_sel))

print(paste("Number of drugs:", length(unique(drugs_sel))))

print(paste("Number of targets:", length(unique(targets_sel))))
