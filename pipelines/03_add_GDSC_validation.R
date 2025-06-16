#load GDSC data
drugTargetInfo <- read.table(paste(pathdata,'/raw/drug-target_data_hgvs_clean_Goncalves_et_all.txt', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
gdsc1<-read.csv(paste(pathdata,'/raw/GDSC1_fitted_dose_response_25Feb20.csv', sep=""),header = TRUE,stringsAsFactors = FALSE)
gdsc2<-read.csv(paste(pathdata,'/raw/GDSC2_fitted_dose_response_25Feb20.csv', sep=""),header = TRUE,stringsAsFactors = FALSE)
colnames(gdsc1)[1]<-"DATASET"
colnames(gdsc2)[1]<-"DATASET"

gdscAll<-rbind(gdsc1,gdsc2)

clTiss<-CMP_annot$model_name[CMP_annot$cancer_type==ctiss]

#Load DAM info
load(paste(resultPath,'/',ctiss,'_results.RData',sep=''))

#create results folder
if(!dir.exists(paste(resultPath,'/_DR_plots',sep=''))){
  dir.create(paste(resultPath,'/_DR_plots',sep=''))
}

if(!dir.exists(paste(resultPath,'/_DR_plots/',ctiss,sep=''))){
  dir.create(paste(resultPath,'/_DR_plots/',ctiss,sep=''))
}

#computing the rank ratio for each drug targeting a DAM-bearing-gene
RES<-lapply(1:nrow(RESTOT),function(x){
    
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
      
      ud<-drug_names
      udid<-drug_ids
      
      zscores<-matrix(NA,nrow = length(udid),ncol = length(cellLines),dimnames = list(udid,cellLines))
      
      rankRatio<-matrix(NA,nrow = length(udid),ncol = length(cellLines),dimnames = list(udid,cellLines))
      
      lnIC50<-matrix(NA,nrow = length(udid),ncol = length(cellLines),dimnames = list(udid,cellLines))
      perctl<-matrix(NA,nrow = length(udid),ncol = length(cellLines),dimnames = list(udid,cellLines))
      concRatio<-matrix(NA,nrow = length(udid),ncol = length(cellLines),dimnames = list(udid,cellLines))
      additionalInfos<-as.data.frame(matrix(NA,nrow=length(udid),ncol=7))
      colnames(additionalInfos)<-c('screen','drug_id','drug_name','put_target','min_conc','max_conc','rank ratio')
      rownames(additionalInfos)<-udid
      
      if(nrow(data1)>0){
        for (i in 1:nrow(data1)){
          zscores[as.character(data1$DRUG_ID)[i],data1$CELL_LINE_NAME[i]]<-data1$Z_SCORE[i]
          lnIC50[as.character(data1$DRUG_ID)[i],data1$CELL_LINE_NAME[i]]<-data1$LN_IC50[i]
          
          additionalInfos[as.character(data1$DRUG_ID)[i],'screen']<-data1$DATASET[i]
          additionalInfos[as.character(data1$DRUG_ID)[i],'drug_id']<-data1$DRUG_ID[i]
          additionalInfos[as.character(data1$DRUG_ID)[i],'drug_name']<-data1$DRUG_NAME[i]
          additionalInfos[as.character(data1$DRUG_ID)[i],'put_target']<-data1$PUTATIVE_TARGET[i]
          additionalInfos[as.character(data1$DRUG_ID)[i],'min_conc']<-data1$MIN_CONC[i]
          additionalInfos[as.character(data1$DRUG_ID)[i],'max_conc']<-data1$MAX_CONC[i]
          
          allPattern<-gdscAll$LN_IC50[gdscAll$DATASET==data1$DATASET[i] & gdscAll$DRUG_ID==data1$DRUG_ID[i] & is.element(gdscAll$CELL_LINE_NAME,clTiss)]
          names(allPattern)<-gdscAll$CELL_LINE_NAME[gdscAll$DATASET==data1$DATASET[i] & gdscAll$DRUG_ID==data1$DRUG_ID[i] & is.element(gdscAll$CELL_LINE_NAME,clTiss)]
          
          hits<-match(cellLines,names(allPattern)[order(allPattern)])
          n<-length(which(!is.na(allPattern[hits])))
          
          additionalInfos[as.character(data1$DRUG_ID)[i],'rank ratio']<-sum(hits,na.rm = TRUE)/((n*(n+1))/2)
          
          perctl[as.character(data1$DRUG_ID)[i],data1$CELL_LINE_NAME[i]]<-round(100*match(data1$LN_IC50[i],sort(allPattern))/length(allPattern),2)
          concRatio[as.character(data1$DRUG_ID)[i],data1$CELL_LINE_NAME[i]]<- data1$LN_IC50[i]/log(data1$MAX_CONC[i])
          
          bg<-rep(rgb(0,0,255,alpha = 110,maxColorValue = 255),length(allPattern))
          names(bg)<-names(allPattern)
          bg[cellLines]<-'red'
          
          
          dname<-data1$DRUG_NAME[i]
          dname<-str_replace_all(dname,'/','|')
          pdf(gsub("\\|", "",paste( resultPath,'/_DR_plots/',ctiss,'/',target,' - ',
                    data1$DATASET[i],'_',
                    data1$DRUG_ID[i],'_',
                    dname,'.pdf',sep='')),8.85,4.20)
          
          par(mfrow=c(1,2))
          plot(sort(allPattern),
               col=bg[order(allPattern)],xlab='cell lines',ylab='ln IC50',pch=16,
               main = paste(data1$DRUG_NAME[i],' [',target,'] - rankRatio = ',round(sum(hits,na.rm = TRUE)/((n*(n+1))/2),2),sep=''))
          abline(h=log(data1$MAX_CONC[i]),col='gray',lty=2)
          abline(h=log(data1$MAX_CONC[i]),col='gray',lty=2)
          
          plot(0,0,frame.plot = FALSE,axes = FALSE,col=NA,xlab='',ylab='')
          legend('left',pch=16,col=c('red',"#0000FF6E"),legend=c(variant,'others'),bg = 'white',title=paste(target,'status'))
          dev.off()
          
        }
        
        SCREENdata<-list(screenInfo=additionalInfos,
                         lnIC50=lnIC50,
                         Zscores=zscores,
                         perctl=perctl,
                         concRatio=concRatio)
        
        if(RESTOT$rank_ratio[x]<1.6 & RESTOT$medFitEff[x]< -0.5 & min(SCREENdata[["screenInfo"]]$'rank ratio', na.rm=T)<1.6){bb<-'match'}else{bb<-''}
        
        save(SCREENdata,file=paste(resultPath,'/_DR_plots/',ctiss,'/',target,' _ ',bb,paste(gsub("\\?|\\*|!|\\|", "", variant),collapse='AND'),'_screenRes.RData',sep=''))
        
      
        
        }
          
      }
    
  })


