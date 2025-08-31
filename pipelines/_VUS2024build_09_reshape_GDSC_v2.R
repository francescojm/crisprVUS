set.seed(123)
library(CELLector)
library(tidyverse)

####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250808_bugFixed_and_RR_th.1.71_wr/'

inTOgen_drivers<-read.table(paste(pathdata,"/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", sep=""), sep='\t',stringsAsFactors = FALSE,header=TRUE)
inTOgen_drivers<-unique(sort(inTOgen_drivers$SYMBOL))

decoupleMultipleHits<-function(hitTable){
  vars<-hitTable$Hit
  
  genes<-str_split_fixed(vars,' ',2)
  
  vars<-genes[,2]
  genes<-genes[,1]
  
  iimultiple<-grep(' | ',vars)
  
  ehitTable<-cbind(hitTable[,1],genes,vars,hitTable[,c(1,3:ncol(hitTable))])
  
  finalHit<-ehitTable[setdiff(1:nrow(hitTable),iimultiple),]
  
  finalHit<-rbind(finalHit,
                  
                  do.call('rbind',lapply(iimultiple,function(x){
                    
                    
                    indivar<-setdiff(strsplit(vars[x],' | ')[[1]],'|')
                    tmpHit<-NULL
                    for (i in 1:length(indivar)){
                      tmpHit<-rbind(tmpHit,ehitTable[x,])
                    }
                    
                    tmpHit$vars<-indivar
                    
                    return(tmpHit)
                  })))
  
  colnames(finalHit)[1]<-'cancer_type'
  return(finalHit)
}

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


###################################################
####################################################

# allDRvaliations<-NULL
# 
# for(ctiss in tissues){
#   load(paste(resultPath,'/',ctiss,'_results.RData',sep=''))
#   
#   hits<-RESTOT[RESTOT$medFitEff< -0.5 & RESTOT$rank_ratio<1.71 & RESTOT$hypTest_p<0.20 & RESTOT$empPval<0.20,]
#   
#   allTar<-hits$GENE
#   allVar<-hits$var
#   
#   allfnames<-list.files(pattern='RData',path=paste(resultPath,'/_DR_plots/',ctiss,"/", sep=''))
#   fnames<-unlist(lapply(str_split(allfnames,' _ '),function(x){x[1]}))
#   
#   allVar<-allVar[which(allTar %in% fnames)]
#   allTar<-intersect(allTar,fnames)
#   
#   if(length(allTar)>0){
#     rRES<-do.call(rbind,lapply(1:length(allTar),function(i){
#       
#       x<-allTar[i]
#       print(x)
#       current_fn<-allfnames[match(x,unlist(lapply(str_split(allfnames,' _ '),function(x){x[[1]][1]})))]
#       load(paste(resultPath,'/_DR_plots/',ctiss,'/',current_fn,sep=''))
#       
#       nnd<-nrow(SCREENdata$screenInfo)
#       
#       RES<-cbind(rep(ctiss,nnd),rep(paste(x,allVar[i]),nnd),SCREENdata$screenInfo)
#       
#       colnames(RES)[c(1,2)]<-c('ctype','Hit')
#       rownames(RES)<-NULL  
#       return(RES)
#     }))
#     
#     save(rRES,file=paste(resultPath,'/_DR_plots/',ctiss,'_DR_validation.RData',sep=''))
#     write.table(rRES,quote=FALSE,sep='\t',
#                 row.names = FALSE,file=paste(resultPath,'/_DR_plots/',ctiss,'_DR_validation.tsv',sep=''))
#     allDRvaliations<-rbind(allDRvaliations,rRES)
#   }
# }
# 
# allDRvaliations<-allDRvaliations[!is.na(allDRvaliations$validated),]
# save(allDRvaliations,file=paste(resultPath,'/_all_DR_validations.RData',sep=''))
# write.table(allDRvaliations,quote=FALSE,sep='\t',
#             row.names = FALSE,file=paste(resultPath,'_all_DR_validations.tsv',sep=''))


load(paste(resultPath,'/_all_DR_validations.RData',sep=''))
load(paste(resultPath,'/_allHits.RData',sep=''))

sigs<-paste(allDRvaliations$ctype,allDRvaliations$Hit)

ncts_hits<-length(unique(paste(allHits$ctype,allHits$GENE,allHits$var)))
nvalidable<-length(unique(sigs))

print(paste('of the ',ncts_hits,' cancer-type-specific hits (DAMs or DAMs combinations with optimal RankRatio, fitness effect and pvalues), ',nvalidable,
            ' involve a DAM-bearing genes that is druggable and targeted by a compound with available cancer-type matching drug-response data on GDSC',sep=''))

validated<-unlist(lapply(unique(sigs),function(x){
  ii<-which(sigs==x)
  sum(allDRvaliations[ii,]$validated)
}))

nvalidated <- length(unique(sigs)[which(validated>0)])

print(paste('of these, ',nvalidated,' are validated',sep=''))

allDRvaliations<-allDRvaliations[allDRvaliations$validated==TRUE,]

allSAMs<-decoupleMultipleHits(hitTable = allDRvaliations)
allSAMs<-allSAMs[order(allSAMs$cancer_type),]


print(paste('encompassing ',
            length(unique(paste(allSAMs$cancer_type,allSAMs$genes,allSAMs$vars))),
            ' individual cancer-type-specific drug validated DAMs (SAMs)',sep=''))

ii<-which(is.element(allSAMs$genes,inTOgen_drivers))
kallSAMs<-allSAMs[ii,]

print(paste('encompassing ',
            length(unique(paste(kallSAMs$cancer_type,kallSAMs$genes,kallSAMs$vars))),
            ' individual cancer-type-specific drug validated DAMs (SAMs)',sep=''))

inTOgen_drivers<-read.table(paste(pathdata,"/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", sep=""), sep='\t',stringsAsFactors = FALSE,header=TRUE)
ctypeMapping<-read.csv(paste(pathdata,'/raw/intOGen ctype mapping_AS.csv',sep=''),header = TRUE,row.names = 1, sep=";")


nSAMs<-nrow(allSAMs)

isAcancerDriver<-rep(NA,nSAMs)
isA_ct_specific_Driver<-rep(NA,nSAMs)
for (i in 1:nSAMs){
    isAcancerDriver[i]<-is.element(allSAMs$genes[i],inTOgen_drivers$SYMBOL)
    intTypes<-str_trim(unlist(str_split(ctypeMapping[allSAMs$cancer_type[i],1],'\\|')))
    isA_ct_specific_Driver[i]<-is.element(allSAMs$genes[i],inTOgen_drivers$SYMBOL[which(is.element(inTOgen_drivers$CANCER_TYPE,intTypes))])
}

repurposableDrug<-(!isAcancerDriver | !isA_ct_specific_Driver)

allSAMs<-cbind(allSAMs,isAcancerDriver,isA_ct_specific_Driver,repurposableDrug)

save(allSAMs,file=paste(resultPath,'/_all_SAMs.RData',sep=''))
write.table(allSAMs,quote=FALSE,sep='\t',
            row.names = FALSE,file=paste(resultPath,'_all_SAMs.tsv',sep=''))


