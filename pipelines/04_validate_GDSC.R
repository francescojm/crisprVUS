set.seed(123)
library(CELLector)
library(tidyverse)

####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250221/'

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

tissues<-sort(tissues)

load(paste(pathdata, "/Robj/basal_exp.RData", sep=""))


###################################################
####################################################

for(ctiss in tissues){
load(paste(resultPath,'/',ctiss,'_results.RData',sep=''))

hits<-RESTOT[RESTOT$medFitEff< -0.5 & RESTOT$rank_ratio<1.6,]

allTar<-hits$GENE
allVar<-hits$var

allfnames<-list.files(pattern='RData',path=paste(resultPath,'/_DR_plots/',ctiss,"/", sep=''))
fnames<-unlist(lapply(str_split(allfnames,' _ '),function(x){x[1]}))

allVar<-allVar[which(allTar %in% fnames)]
allTar<-intersect(allTar,fnames)

if(length(allTar)>0){
  

  rRES<-do.call(rbind,lapply(1:length(allTar),function(i){
    
    x<-allTar[i]
    print(x)
    current_fn<-allfnames[match(x,unlist(lapply(str_split(allfnames,' _ '),function(x){x[[1]][1]})))]
    load(paste(resultPath,'/_DR_plots/',ctiss,'/',current_fn,sep=''))
    
    bestRankRatio<-min(SCREENdata$screenInfo$`rank ratio`,na.rm = TRUE)
    
    bestRR_Screen<-
      paste(SCREENdata$screenInfo[which(SCREENdata$screenInfo$`rank ratio`==bestRankRatio),'screen'],collapse=' | ')
    bestRR_drug_id<-
      paste(SCREENdata$screenInfo[which(SCREENdata$screenInfo$`rank ratio`==bestRankRatio),'drug_id'],collapse=' | ')
    bestRR_drug_name<-
      paste(SCREENdata$screenInfo[which(SCREENdata$screenInfo$`rank ratio`==bestRankRatio),'drug_name'],collapse=' | ')
    
    trup5<-sum(rowSums(SCREENdata$perctl< 5,na.rm = TRUE))>0
    trup10<-sum(rowSums(SCREENdata$perctl< 10,na.rm = TRUE))>0
    trup25<-sum(rowSums(SCREENdata$perctl< 25,na.rm = TRUE))>0
    trup50<-sum(rowSums(SCREENdata$perctl< 50,na.rm = TRUE))>0
    maxConcMat<-matrix(rep(log(SCREENdata$screenInfo$max_conc),ncol(SCREENdata$lnIC50)),nrow(SCREENdata$lnIC50),ncol(SCREENdata$lnIC50))
  
    trupMAXc<-sum(rowSums(SCREENdata$lnIC50 < maxConcMat),na.rm = TRUE)>0
    trupMAXc.2<-sum(rowSums(SCREENdata$lnIC50 < log(exp(maxConcMat)/2)),na.rm = TRUE)>0
    trupMAXc.4<-sum(rowSums(SCREENdata$lnIC50 < log(exp(maxConcMat)/4)),na.rm = TRUE)>0
    trupMAXc.10<-sum(rowSums(SCREENdata$lnIC50 < log(exp(maxConcMat)/10)),na.rm = TRUE)>0
    trupMAXc.100<-sum(rowSums(SCREENdata$lnIC50 < log(exp(maxConcMat)/100)),na.rm = TRUE)>0
    trupMAXc.1000<-sum(rowSums(SCREENdata$lnIC50 < log(exp(maxConcMat)/1000)),na.rm = TRUE)>0
    
    RES<-data.frame(t(c(allVar[i],
                        bestRankRatio,
                        bestRR_Screen,
                        bestRR_drug_id,
                        bestRR_drug_name,
                        trup5,trup10,trup25,trup50,trupMAXc,trupMAXc.2,trupMAXc.4,trupMAXc.10,trupMAXc.100,trupMAXc.1000)),stringsAsFactors = FALSE)
    colnames(RES)<-c('Var',
                     'bestRR',
                     'bestRR_Screen',
                     'bestRR_drug_id',
                     'bestRR_drug_name',
                     '5perc','10perc','25perc','50perc','lMAXc','l.5MAXc','l.25MAXc','l.1MAXc','l.01MAXc','l.001MAXc')
    rownames(RES)<-x  
    return(RES)
  }))
  
  save(rRES,file=paste(resultPath,'/_DR_plots/',ctiss,'_DR_validation.RData',sep=''))
  
  rRES<-cbind(rownames(rRES),rRES)
  colnames(rRES)[1]<-'Gene'
  write.table(rRES,quote=FALSE,sep='\t',
              row.names = FALSE,file=paste(resultPath,'/_DR_plots/',ctiss,'_DR_validation.tsv',sep=''))
}

}

