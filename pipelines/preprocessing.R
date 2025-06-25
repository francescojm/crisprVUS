set.seed(123)
library(CELLector)
library(tidyverse)

####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250221/'

# ###RNA-seq data downloaded from https://cog.sanger.ac.uk/cmp/download/rnaseq_merged_20250117.zip on 20250214
# basal_exp<-read.csv(paste(pathdata, "/raw/rnaseq_merged_20250117/rnaseq_merged_rsem_fpkm_20250117.csv", sep=""))
# basal_exp<-basal_exp[-c(1:3),]
# basal_exp<-basal_exp[-39330,]
# rownames(basal_exp)<-basal_exp[,1]
# basal_exp<-basal_exp[,-c(1:3)]
# genes_exp<-rownames(basal_exp)
# basal_exp<-sapply(basal_exp, as.numeric)
# rownames(basal_exp)<-genes_exp
# basal_exp<-log2(basal_exp+1)
# save(basal_exp, file=paste(pathdata, "/Robj/basal_exp.RData", sep=""))

#########Francesco's file
load(file=paste(pathdata, "/Robj/basa_exp_FPKM.RData", sep=""))
CMP_annot <- read_csv(paste(pathdata,"/raw/model_list_20241120.csv", sep="")) 
### model_list_20240110.csv downloaded from https://cog.sanger.ac.uk/cmp/download/model_list_20241120.csv on 20250129
colnames(basal_exp)<-CMP_annot$model_id[match(colnames(basal_exp), CMP_annot$model_name)]
save(basal_exp, file=paste(pathdata, "/Robj/basal_exp.RData", sep=""))

library(openxlsx)
proteome<-read.delim("data/raw/Proteomics_20250211/Protein_matrix_averaged_20250211.tsv")
#downloaded from https://cog.sanger.ac.uk/cmp/download/Proteomics_20250211.zip on the 18/03/2025

rownames(proteome)<-proteome[,2]
proteome<-proteome[,-c(1,2)]
colnames(proteome)<-proteome[1,]
proteome<-proteome[-c(1,2),]
proteome<-t(proteome)
proteome_mat <- matrix(as.numeric(unlist(proteome)),    # Convert to numeric matrix
                       ncol = ncol(proteome))
colnames(proteome_mat)<-colnames(proteome)
rownames(proteome_mat)<-rownames(proteome)
save(proteome_mat, file=paste(pathdata, "/Robj/proteome_mat.RData", sep=""))


###loading input data
gdsc1<-read.csv(paste(pathdata,'/raw/GDSC1_fitted_dose_response_27Oct23.csv', sep=""),header = TRUE,stringsAsFactors = FALSE)
gdsc2<-read.csv(paste(pathdata,'/raw/GDSC2_fitted_dose_response_27Oct23.csv', sep=""),header = TRUE,stringsAsFactors = FALSE)
### GDSC1_fitted_dose_response_27Oct23.csv and GDSC2_fitted_dose_response_27Oct23.csv have been downloaded from: 
### https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC1_fitted_dose_response_27Oct23.xlsx and https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx
### on the 20241003
drugTargetInfo <- read.table(paste(pathdata,'/raw/drug-target_data_hgvs_clean_Goncalves_et_all.txt', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
### drug-target_data_hgvs_clean_Goncalves_et_all.txt built from https://www.embopress.org/doi/suppl/10.15252/msb.20199405/suppl_file/msb199405-sup-0003-datasetev2.xlsx on 20241003

gdscAll<-rbind(gdsc1, gdsc2)

drugs<-as.character(unique(gdscAll$DRUG_ID))
cell_lines<-unique(gdscAll$SANGER_MODEL_ID)
drugs_lnIC50<-matrix(nrow=length(drugs), ncol=length(cell_lines), 
                     dimnames=list(c(drugs), c(cell_lines)))
for(drug in drugs){
  for(cell_line in cell_lines){
    ind<-which(gdscAll$DRUG_ID==drug & gdscAll$SANGER_MODEL_ID==cell_line)
    if(length(ind)==1){
      drugs_lnIC50[drug, cell_line]<-gdscAll$LN_IC50[ind]
    } else if(length(ind)>1){
      print(drug)
      print(cell_line)
      drugs_lnIC50[drug, cell_line]<-log(mean(exp(gdscAll$LN_IC50[ind])))
    }
  }
}
save(drugs_lnIC50, file=paste(pathdata, "/Robj/drugs_lnIC50_gdscAll.RData", sep=""))

