
inTOgen_drivers<-read.table(paste(pathdata,"/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv", sep=""), sep='\t',stringsAsFactors = FALSE,header=TRUE)
load(paste(resultPath,'/',ctiss,'_results.RData',sep=''))

inTOgen_annot<-do.call(rbind,lapply(RESTOT$GENE,function(x){
  ids<-which(inTOgen_drivers$SYMBOL==x)

  inTOgen_driver_for<-paste(sort(unique(inTOgen_drivers$CANCER_TYPE[ids])), collapse = ' | ')
  
  if(inTOgen_driver_for==""){
    inTOgen_driver_for<-'none'
  }
  
  inTOgen_Role_Perc<-c(0,0,0)
  names(inTOgen_Role_Perc)<-c('Act','LoF','ambiguous')
  rolSums<-summary(as.factor(inTOgen_drivers$ROLE[ids]))
  inTOgen_Role_Perc[names(rolSums)]<-round(100*rolSums/sum(rolSums),2)
  
  data.frame(inTOgen_driver_for=inTOgen_driver_for,Act=inTOgen_Role_Perc[1],LoF=inTOgen_Role_Perc[2],ambigous=inTOgen_Role_Perc[3],row.names = x)
}))

RESTOT<-cbind(RESTOT,inTOgen_annot)
save(RESTOT,file=paste(resultPath,'/',ctiss,'_results.RData',sep=''))
write.table(RESTOT,quote=FALSE,sep='\t',row.names = FALSE,file=paste(resultPath,'/',ctiss,'_results.tsv',sep=''))
