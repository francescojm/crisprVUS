
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

save(incl_cl_annot,file=paste(resultPath,'/_incl_cl_annot.RData',sep=''))

tissues<-sort(tissues)

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

