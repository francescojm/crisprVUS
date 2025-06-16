
data_path<-'../data/Data_2025/data/'


## DepMap CRISPR data
CRISPRGeneEffect<-read.csv(paste(data_path,'raw/CRISPRGeneEffect.csv',sep=''))


## CMP variants
cl_variants <- read.csv('../data/_VUS2024build/raw/mutations_all_20230202.csv')


## DepMap Models' annotations
DepMapModels<-read.csv(paste(data_path,'raw/Model.csv',sep=''))

## Cell Model Passports models' annotations
CMP_annot <- read.csv(paste(pathdata,"/raw/model_list_20240110.csv", sep="")) 


## Compute n. of Brad DepMap with valid CMP id:
sum(!is.na(CMP_annot$model_id[match(CRISPRGeneEffect$X,CMP_annot$BROAD_ID)]))