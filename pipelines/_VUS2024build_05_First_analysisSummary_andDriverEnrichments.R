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



decoupleMultipleHits<-function(hitTable){
  vars<-hitTable$var
  iimultiple<-grep(' | ',vars)
  
  finalHit<-hitTable[setdiff(1:nrow(hitTable),iimultiple),]
  
  finalHit<-rbind(finalHit,
                  do.call('rbind',lapply(iimultiple,function(x){
                    indivar<-setdiff(strsplit(vars[x],' | ')[[1]],'|')
                    tmpHit<-NULL
                    for (i in 1:length(indivar)){
                      tmpHit<-rbind(tmpHit,hitTable[x,])
                    }
                    tmpHit$var<-indivar
                    return(tmpHit)
                  })))
  
  }
my.hypTest<-function(x,k,n,N){
  
  PVALS<-phyper(x-1,n,N-n,k,lower.tail=FALSE)
  
  return(PVALS)
}


print(paste(nrow(incl_cl_annot),'cell lines included in the analysis'))
print(paste(length(unique(incl_cl_annot$cancer_type)),'considered cancer types'))

ss<-summary(as.factor(incl_cl_annot$cancer_type))
print(paste('median number of cell lines per cancer type =',median(ss)))
print(paste('min = ', min(ss),'for',names(sort(ss))[1]))
print(paste('max = ', max(ss),'for',names(sort(ss,decreasing=TRUE))[1]))
pdf(paste(resultPath, "exploration/figures/Ncelllines.pdf",sep=""), 7,10)
par(mar=c(4,16,0,2))
barplot(sort(ss),horiz = TRUE,las=2,border=FALSE,col='blue',xlab='n. cell lines')
abline(v=median(ss),lty=2)
dev.off()


 #uncomment this to collect all tested variants
 fc<-dir(resultPath)
 fc<-grep('testedVariants.RData',fc,value=TRUE)
 
 totalTestedVariants<-NULL
 for (i in 1:length(fc)){
   print(i)
   cty<-strsplit(fc[i],'_testedVariants.RData')[[1]]
   load(paste(resultPath,'/',fc[i],sep=''))
 
   variant_UMIs<-paste(ts_cl_variants$gene_symbol,ts_cl_variants$protein_mutation)
   unique_variants<-unique(variant_UMIs)
 
   curRes<-do.call('rbind',lapply(unique_variants,function(x){
       idxs<-which(variant_UMIs==x)
       positiveCls<-paste(ts_cl_variants$model_id[idxs],collapse=', ')
       res<-cbind(cty,x,ts_cl_variants[idxs[1],c(1,2,4,5,6,7,8,9,10,11,12)],positiveCls)
       }))
   totalTestedVariants<-rbind(totalTestedVariants,curRes)
   }
 
 save(totalTestedVariants,file=paste(resultPath,'_totalTestedVariants.RData',sep=''))
 write.table(totalTestedVariants,sep="\t",quote=FALSE,file=paste(resultPath,'_totalTestedVariants.tsv',sep=''))
 #[END] uncomment this to collect all tested variants

load(paste(resultPath,'_totalTestedVariants.RData',sep=''))
ntestedVarCtypeCombos<-nrow(totalTestedVariants)
print(paste(ntestedVarCtypeCombos,'tested variants x cancer type combos'))
ntestedIndividualVariants<-length(unique(totalTestedVariants$x))
print(paste(ntestedIndividualVariants,'tested individual variants'))
ntestedGenes<-length(unique(totalTestedVariants$gene_symbol))
print(paste('involving',ntestedGenes,'genes'))

library(stringr)
var_nlines <- str_count(totalTestedVariants$positiveCls, ",")+1
totalTestedVariants$var_nlines<-var_nlines

pdf(paste(resultPath, "exploration/figures/var_nlines.pdf",sep=""), 7,5)
ggplot(totalTestedVariants, aes(x = var_nlines)) +
  geom_histogram(color = "#000000", fill = "steelblue")+theme_classic()+
  labs(x = "# of cell lines bearing the variant",y = "# of variants")
dev.off()
pdf(paste(resultPath, "exploration/figures/var_nlines_log10.pdf",sep=""), 7,5)
ggplot(totalTestedVariants, aes(x = var_nlines)) +
  geom_histogram(color = "#000000", fill = "steelblue")+theme_classic()+ scale_y_log10()+
  labs(x = "# of cell lines bearing the variant",
    y = "# of variants (log10)"
  )
dev.off()

ct_tested_variants<-unlist(lapply(tissues,function(x){length(which(totalTestedVariants$cty==x))}))
names(ct_tested_variants)<-tissues

 #uncomment this to collect all DAMs
 fc<-dir(resultPath)
 fc<-grep('_results_ext.RData',fc,value=TRUE)
 
 allDAMs<-NULL
 allHits<-NULL
 for (i in 1:length(fc)){
    print(i)
    cty<-strsplit(fc[i],'_results_ext.RData')[[1]]
    load(paste(resultPath,'',fc[i],sep=''))
 
    hitsIdxs<-which(RESTOT$medFitEff< -0.5 & RESTOT$rank_ratio< 1.6 & RESTOT$pval_rand < 0.2)
    currHits<-RESTOT[hitsIdxs,]
    currDAMs<-decoupleMultipleHits(currHits)
    rownames(currHits)<-NULL
    rownames(currDAMs)<-NULL
 
    currHits<-currHits[order(currHits$GENE),]
    currDAMs<-currDAMs[,1:3]
    currDAMs<-currDAMs[order(currDAMs$GENE),]
    allDAMs<-rbind(allDAMs,currDAMs)
    allHits<-rbind(allHits,currHits)
    }
 
 
 
 save(allHits,file=paste(resultPath,'_allHits.RData',sep=''))
 write.table(allHits,sep="\t",quote=FALSE,file=paste(resultPath,'_allHits.tsv',sep=''))
 save(allDAMs,file=paste(resultPath,'_allDAMs.RData',sep=''))
 write.table(allDAMs,sep="\t",quote=FALSE,file=paste(resultPath,'_allDAMs.tsv',sep=''))
 allDAM_bearing_genes<-sort(unique(allHits$GENE))
 id<-match(allDAM_bearing_genes,allHits$GENE)
 
 
 DAMbearing_in<-unlist(lapply(allDAM_bearing_genes,function(x){
     cid<-which(allHits$GENE==x)
     analyses<-paste(sort(allHits$ctype[cid]),collapse=' | ')
   }))
 
 allDAM_bearing_genes<-cbind(allDAM_bearing_genes,allHits[id,c('inTOgen_driver_for','Act','LoF','ambigous')],DAMbearing_in)
 save(allDAM_bearing_genes,file=paste(resultPath,'_allDAM_bearing_genes.RData',sep=''))
 write.table(allDAM_bearing_genes,quote=FALSE,sep='\t',file=paste(resultPath,'_allDAM_bearing_genes.tsv',sep=''))
 #[END] uncomment this to collect all DAMS

load(paste(resultPath,'_allHits.RData',sep=''))
load(paste(resultPath,'_allDAMs.RData',sep=''))
load(paste(resultPath,'_allDAM_bearing_genes.RData',sep=''))
print(paste(nrow(allHits),'hits with rankratio < 1.6, medFitness effect < -0.5  and p < 0.2 across cancer types'))
print(paste(nrow(allDAMs),'individual DAMs across cancer types'))

print(paste('corresponding to',length(unique(paste(allDAMs$GENE,allDAMs$var))),'individual variants'))
print(paste(100*length(paste(allDAMs$GENE,allDAMs$var))/nrow(totalTestedVariants),'% of tested cases (considering variants testes multiple times across different cancer types)'))
print(paste('involving ',length(unique(allDAMs$GENE)),'genes'))
print(paste(100*length(unique(paste(allDAMs$GENE)))/ntestedGenes,'% of tested cases'))

pdf(paste(resultPath, "exploration/figures/DAMs_ctype.pdf",sep=""), 12,10)
par(mar=c(20,5,5,5))
barplot(table(allDAMs$ctype)[order(table(allDAMs$ctype), decreasing=T)], las=2, ylab="# of DAMs") 
dev.off()
pdf(paste(resultPath, "exploration/figures/DAMsbearing_ctype.pdf",sep=""), 12,10)
par(mar=c(20,5,5,5))
allDAMbearing<-allDAMs[!duplicated(allDAMs[,c('ctype','GENE')]),]
barplot(table(allDAMbearing$ctype)[order(table(allDAMbearing$ctype), decreasing=T)], las=2, ylab="# of DAM-bearing genes") 
dev.off()


nDAMbearing<-c()
for(t in tissues){
  allDAMs_sel<-allDAMs[allDAMs$ctype==t,]
  nDAMbearing<-c(nDAMbearing, length(unique(paste(allDAMs_sel$GENE))))
}
names(nDAMbearing)<-tissues
print(paste(median(nDAMbearing),'median number of DAM bearing genes across cancer types'))
print(paste('min = ', min(nDAMbearing),'for',names(sort(nDAMbearing))[1]))
print(paste('max = ', max(nDAMbearing),'for',names(sort(nDAMbearing, decreasing=T))[1]))
# enrichment of cancer driver genes among the DAM-bearing genes

inTOgen_drivers<-read.table(paste(pathdata,"/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", sep=""), sep='\t',stringsAsFactors = FALSE,header=TRUE)
inTOgen_drivers<-unique(sort(inTOgen_drivers$SYMBOL))


N<-ntestedGenes
n<-length(inTOgen_drivers)
k<-length(unique(allDAM_bearing_genes$allDAM_bearing_genes))
x<-length(intersect(inTOgen_drivers,allDAM_bearing_genes$allDAM_bearing_genes))

print(paste('of the',k,'DAM bearing genes,',x,'have been previously reported as cancer drivers (inTOgen) (',100*x/k,'%)'))
print('hypergeometric test p-value:')
my.hypTest(x,k,n,N)

DAMbearing_known_as_cancer_driver<-x

inTOgen_drivers<-read.table(paste(pathdata,"/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", sep=""), sep='\t',stringsAsFactors = FALSE,header=TRUE)

ug<-sort(unique(inTOgen_drivers$SYMBOL))

roles<-do.call('rbind',lapply(ug,function(x){
  ids<-which(inTOgen_drivers$SYMBOL==x)
  inTOgen_Role_Perc<-c(0,0,0)
  names(inTOgen_Role_Perc)<-c('Act','LoF','ambiguous')
  rolSums<-summary(as.factor(inTOgen_drivers$ROLE[ids]))
  inTOgen_Role_Perc[names(rolSums)]<-round(100*rolSums/sum(rolSums),2)
  return(inTOgen_Role_Perc)
}))

rownames(roles)<-ug
putRole<-rep(NA,length(ug))
for (i in 1:nrow(roles)){
  putRole[i]<-colnames(roles)[order(roles[i,],decreasing=TRUE)[1]]  
}
names(putRole)<-ug


oncogenes<-names(which(roles[,'Act']==100))
tsg<-names(which(roles[,'LoF']==100))
ambig<-setdiff(rownames(roles),c(oncogenes,tsg))

n<-length(oncogenes)
x<-length(intersect(oncogenes,allDAM_bearing_genes$allDAM_bearing_genes))
print(paste('of the',k,'DAM bearing genes,',x,'have been previously reported as 100% oncoGenes (inTOgen) (',100*x/k,'%)'))
print('hypergeometric test p-value:')
my.hypTest(x,k,n,N)

nAct<-x
pAct<-my.hypTest(x,k,n,N)


n<-length(tsg)
x<-length(intersect(tsg,allDAM_bearing_genes$allDAM_bearing_genes))
print(paste('of the',k,'DAM bearing genes,',x,'have been previously reported as 100% TSG (inTOgen) (',100*x/k,'%)'))
print('hypergeometric test p-value:')
my.hypTest(x,k,n,N)

nTsg<-x
pTsg<-my.hypTest(x,k,n,N)


n<-length(ambig)
x<-length(intersect(ambig,allDAM_bearing_genes$allDAM_bearing_genes))
print(paste('of the',k,'DAM bearing genes,',x,'have been previously reported as TSGs or oncoGenes (inTOgen) (',100*x/k,'%)'))
print('hypergeometric test p-value:')
my.hypTest(x,k,n,N)

nAmb<-x
pAmb<-my.hypTest(x,k,n,N)


n<-length(ambig)
x<-length(intersect(ambig,allDAM_bearing_genes$allDAM_bearing_genes))
k<-DAMbearing_known_as_cancer_driver
print(paste('of the',k,'DAM bearing genes previously reported to be cancer drivers,',x,'have been previously reported as TSGs or oncoGenes (inTOgen) (',100*x/k,'%)'))
print('hypergeometric test p-value:')
my.hypTest(x,k,n,N)

nAmb<-x
pAmb<-my.hypTest(x,k,n,N)

pdf(paste(resultPath, "exploration/figures/pie_drivers.pdf",sep=""), 7,8)
pie(c(nTsg,nAmb,nAct), col=c('#004add','#d5aff3','red'),border=FALSE,
    main=paste(DAMbearing_known_as_cancer_driver,paste('DAM-bearing genes known as cancer driver\nacross cancer types')),
    labels = c('Frequent LoF','Ambigous','Frequent GoF'),cex=1.2,)
dev.off()

ct_mapping<-read.csv(paste(pathdata,'/raw/intOGen ctype mapping_AS.csv',sep=''),header = TRUE,row.names = 1, sep=";")


COMPOSITION<-NULL
COMPOSITIONp<-NULL

CT_DAM_Bearing<-NULL

for (i in 1:length(tissues)){
  ctiss<-tissues[i]
  current_Dam_bearing<-allHits$GENE[which(allHits$ctype==ctiss)]
  
  CT_DAM_Bearing[[i]]<-current_Dam_bearing
  
  intoTypes<-setdiff(unlist(strsplit(ct_mapping[ctiss,1],' | ')),'|')#intogen cancer type identifier
  currentOG<-unique(inTOgen_drivers[which(is.element(inTOgen_drivers$CANCER_TYPE,intoTypes) & inTOgen_drivers$ROLE=='Act'),'SYMBOL'])#oncogenes from Intogen in the specific cancer type
  currentTSG<-unique(inTOgen_drivers[which(is.element(inTOgen_drivers$CANCER_TYPE,intoTypes) & inTOgen_drivers$ROLE=='LoF'),'SYMBOL'])#tumor suppressor genes from Intogen in the specific cancer type
  both<-intersect(currentOG,currentTSG)#ambiguous OG or TSG from Intogen in the specific cancer type
  currentOG<-setdiff(currentOG,both)#oncogenes (non ambiguous) from Intogen in the specific cancer type
  currentTSG<-setdiff(currentTSG,both)#TSG (non ambiguous) from Intogen in the specific cancer type
  
  classOG<-intersect(current_Dam_bearing,currentOG)
  classAMB<-intersect(current_Dam_bearing,both)
  classTSG<-intersect(current_Dam_bearing,currentTSG)
  novel<-setdiff(current_Dam_bearing,c(classOG,classAMB,classTSG))
  
  RES<-c(length(classTSG),length(classAMB),length(classOG),length(novel))
  bg<-length(unique(totalTestedVariants$gene_symbol[totalTestedVariants$cty==ctiss]))#number of tested genes in the specific tissue
  
  OGp<-my.hypTest(length(classOG),length(current_Dam_bearing),length(currentOG),bg)
  AMBp<-my.hypTest(length(classAMB),length(current_Dam_bearing),length(both),bg)
  TSGp<-my.hypTest(length(classTSG),length(current_Dam_bearing),length(currentTSG),bg)
  
  RESp<-c(TSGp,AMBp,OGp,NA)
  
  COMPOSITION<-cbind(COMPOSITION,RES)
  COMPOSITIONp<-cbind(COMPOSITIONp,RESp)
}

names(CT_DAM_Bearing)<-tissues

oo<-order(colSums(COMPOSITION))
colnames(COMPOSITION)<-tissues
colnames(COMPOSITIONp)<-tissues
rownames(COMPOSITION)<-c('TSG','Amb','OG','Novel')
rownames(COMPOSITIONp)<-c('TSG','Amb','OG','Novel')

pdf(paste(resultPath, "exploration/figures/ActLoFenrichment.pdf", sep=""), 12,9)
par(mfrow=c(1,3))
par(mar=c(4,12,1,2))
barplot(colSums(COMPOSITION)[oo],horiz=TRUE,las=2,xlab='n. DAM-bearing genes',xlim=c(0,200),cex.names=0.6,border=NA,
        col='#75b4d9')
abline(v= median(colSums(COMPOSITION)),lty=2,col='darkgray')

barplot(100*COMPOSITION[1:3,oo]/t(matrix(rep(colSums(COMPOSITION[,oo]),3),ncol(COMPOSITION),3)),
        horiz=TRUE,las=2,xlab='% DAM-bearing genes',xlim=c(0,25),cex.names=0.6,border=NA,
        col=c('#004add','#d5aff3','red'))

COMPOSITIONp[COMPOSITION==0]<-NA
par(mar=c(4,10,1.5,2))
plot(-log10(COMPOSITIONp[1,oo]),1:ncol(COMPOSITIONp),col='#004add',pch=16,xlim=c(0,6), xlab='enrichment -log10 pval',yaxt='n',frame.plot=FALSE,ylab='')
points(-log10(COMPOSITIONp[2,oo]),1:ncol(COMPOSITIONp),col='#d5aff3',pch=16)
points(-log10(COMPOSITIONp[3,oo]),1:ncol(COMPOSITIONp),col='red',pch=16)
abline(v= -log10(0.05),lty=2)
dev.off()

print(paste('median n. of DAM bearing genes across cancer types = ', median(colSums(COMPOSITION)), sep=''))
print(paste('min = ', sort(colSums(COMPOSITION))[1], ' for ', names(sort(colSums(COMPOSITION)))[1], sep=''))
print(paste('max = ', sort(colSums(COMPOSITION),decreasing = TRUE)[1],
            ' for ', names(sort(colSums(COMPOSITION),decreasing = TRUE))[1], sep=''))

print(paste('DAM-bearing genes enriched for cancer type specific TSGs for ',length(which(COMPOSITIONp[1,]<0.05)),' cancer types (',round(100*length(which(COMPOSITIONp[1,]<0.05))/36,2),'%)',sep=''))
print(paste('DAM-bearing genes enriched for cancer type specific abmgigous drivers for ',length(which(COMPOSITIONp[2,]<0.05)),' cancer types (',round(100*length(which(COMPOSITIONp[2,]<0.05))/36,2),'%)',sep=''))
print(paste('DAM-bearing genes enriched for cancer type specific OG drivers for ',length(which(COMPOSITIONp[3,]<0.05)),' cancer types (',round(100*length(which(COMPOSITIONp[3,]<0.05))/36,2),'%)',sep=''))




