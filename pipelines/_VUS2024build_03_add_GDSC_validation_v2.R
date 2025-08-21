my.hypTest<-function(x,k,n,N){
  
  PVALS<-phyper(x-1,n,N-n,k,lower.tail=FALSE)
  
  return(PVALS)
}

#load GDSC data
drugTargetInfo <- read.table(paste(pathdata,'/raw/drug-target_data_hgvs_clean_Goncalves_et_all.txt', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
### drug-target_data_hgvs_clean_Goncalves_et_all.txt built from https://www.embopress.org/doi/suppl/10.15252/msb.20199405/suppl_file/msb199405-sup-0003-datasetev2.xlsx on 20241003

gdsc1<-read.csv(paste(pathdata,'/raw/GDSC1_fitted_dose_response_27Oct23.csv', sep=""),header = TRUE,stringsAsFactors = FALSE)
gdsc2<-read.csv(paste(pathdata,'/raw/GDSC2_fitted_dose_response_27Oct23.csv', sep=""),header = TRUE,stringsAsFactors = FALSE)
### GDSC1_fitted_dose_response_27Oct23.csv and GDSC2_fitted_dose_response_27Oct23.csv have been downloaded from: 
### https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC1_fitted_dose_response_27Oct23.xlsx and https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx
### on the 20241003

drugTargetInfo<-drugTargetInfo[which(is.element(as.character(drugTargetInfo$Drug.ID),union(as.character(gdsc1$DRUG_ID),as.character(gdsc2$DRUG_ID)))),]

gdsc1$uDRUG_ID<-paste(gdsc1$DRUG_ID,'_gdsc1',sep='')
gdsc2$uDRUG_ID<-paste(gdsc2$DRUG_ID,'_gdsc2',sep='')

gdscAll<-rbind(gdsc1,gdsc2)

clTiss<-CMP_annot$model_id[CMP_annot$cancer_type==ctiss]

#Load DAM info
load(paste(resultPath,'/',ctiss,'_results.RData',sep=''))

#create results folder
if(!dir.exists(paste(resultPath,'/_DR_plots',sep=''))){
  dir.create(paste(resultPath,'/_DR_plots',sep=''))
}

if(!dir.exists(paste(resultPath,'/_DR_plots/',ctiss,sep=''))){
  dir.create(paste(resultPath,'/_DR_plots/',ctiss,sep=''))
}

# populating data matrices, drug annotations object and mutant cell lines
RES<-lapply(1:nrow(RESTOT),function(x){
    
    target<-RESTOT$GENE[x]
    variant<-RESTOT$var[x]
    cellLines<-unlist(str_split(RESTOT$ps_cl[x],', '))
    ids<-which(drugTargetInfo$Gene.Target==target)
     
    drug_ids<-drugTargetInfo$Drug.ID[ids]
    drug_names<-drugTargetInfo$Name[ids]
  
    drug_ids<-sort(unique(gdscAll$uDRUG_ID[is.element(gdscAll$DRUG_ID,drug_ids)]))
    drug_names<-
      gdscAll$DRUG_NAME[match(drug_ids,gdscAll$uDRUG_ID)]
    
    if (length(ids)>0){
      
      print(x)
     
      drug_response_data<-gdscAll[which(is.element(gdscAll$uDRUG_ID,drug_ids) & is.element(gdscAll$SANGER_MODEL_ID,clTiss)),]
      udid<-unique(drug_response_data$uDRUG_ID)
      ucl<-unique(drug_response_data$SANGER_MODEL_ID)
      
      zscores<-matrix(NA,nrow = length(udid),
                         ncol = length(ucl),
                      dimnames = list(udid,ucl))
      
      lnIC50<-matrix(NA,nrow = length(udid),
                     ncol = length(ucl),
                     dimnames = list(udid,ucl))
                     
      if(nrow(drug_response_data)>0){
        for (i in 1:nrow(drug_response_data)){
          zscores[as.character(drug_response_data$uDRUG_ID)[i],drug_response_data$SANGER_MODEL_ID[i]]<-drug_response_data$Z_SCORE[i]
          lnIC50[as.character(drug_response_data$uDRUG_ID)[i],drug_response_data$SANGER_MODEL_ID[i]]<-drug_response_data$LN_IC50[i]
        }
      }
      
      additionalInfos<-as.data.frame(matrix(NA,nrow=length(udid),ncol=6))
      colnames(additionalInfos)<-c('screen','drug_id','drug_name','put_target','min_conc','max_conc')
      rownames(additionalInfos)<-sort(udid)
      
      ii<-match(rownames(additionalInfos),drug_response_data$uDRUG_ID)
      
      additionalInfos[,'screen']<-drug_response_data$DATASET[ii]
      additionalInfos[,'drug_id']<-as.character(drug_response_data$DRUG_ID[ii])
      additionalInfos[,'drug_name']<-drug_response_data$DRUG_NAME[ii]
      additionalInfos[,'put_target']<-drug_response_data$PUTATIVE_TARGET[ii]
      additionalInfos[,'min_conc']<-drug_response_data$MIN_CONC[ii]
      additionalInfos[,'max_conc']<-drug_response_data$MAX_CONC[ii]
      
    
     MutantCellLines<-intersect(cellLines,colnames(lnIC50))
     
     if(length(MutantCellLines)>0){
       
       # Computing SAMs scores       
       ndrugs<-nrow(additionalInfos)
       
       SAMscores<-do.call('rbind',lapply(1:ndrugs,function(y){
         
         drpattern<-sort(lnIC50[rownames(additionalInfos)[y],])
         drpattern<-drpattern[!is.na(drpattern)]
         
         medianLnIC50ofMutant<-median(drpattern[MutantCellLines],na.rm = TRUE)
         
         hits<-match(MutantCellLines,names(drpattern))
         hits<-hits[!is.na(hits)]
         nhits<-length(hits)
         rankRatio<-sum(hits)/(nhits*(nhits+1)/2)
         
         k<-max(hits)
         x<-length(hits)
         n<-length(hits)
         N<-length(which(!is.na(drpattern)))
         
         HG_pval<-my.hypTest(x,k,n,N)
         
         NT<-1000
         
         randRankRatio<-rep(NA,NT)
         for (rt in 1:NT){
           hits<-match(MutantCellLines,names(drpattern[sample(length(drpattern))]))
           hits<-hits[!is.na(hits)]
           nhits<-length(hits)
           randRankRatio[rt]<-sum(hits)/(nhits*(nhits+1)/2)
         }
         
         EMP_pval<-(1+(length(which(randRankRatio<=rankRatio))))/(1+NT)
         
         
         res<-c(medianLnIC50ofMutant,rankRatio,HG_pval,EMP_pval)
       }))
       
       colnames(SAMscores)<-c('medLn50_mutCLs','rankRatio','HG_pval','EMP_pval')
       
       additionalInfos<-cbind(additionalInfos,SAMscores)  
       
       validated<-(additionalInfos$medLn50_mutCLs < log(additionalInfos$max_conc)) & 
         (additionalInfos$rankRatio <= 1.5) &
         (additionalInfos$HG_pval < 0.20) &
         (additionalInfos$EMP_pval < 0.20)
       
       additionalInfos<-cbind(additionalInfos,validated)
       
       SAMs<-additionalInfos[which(additionalInfos$validated),]
       
       if(nrow(SAMs)>0){
        lapply(1:nrow(SAMs),function(p){
         
         allPattern<-lnIC50[rownames(SAMs)[p],]
         bg<-rep(rgb(0,0,255,alpha = 110,maxColorValue = 255),length(allPattern))
         names(bg)<-names(allPattern)
         bg[MutantCellLines]<-'red'
         
         dname<-SAMs$drug_name[p]
         dname<-str_replace_all(dname,'/','|')
         
         pdf(gsub("\\|", "",paste( resultPath,'/_DR_plots/',ctiss,'/',target,' - ',
                                   SAMs$screen[p],'_',
                                   SAMs$drug_id[p],'_',
                                   dname,'.pdf',sep='')),8.85,4.20)
         
         par(mfrow=c(1,2))
         plot(sort(allPattern),
              col=bg[order(allPattern)],xlab='cell lines',ylab='ln IC50',pch=16,
              main = paste(dname,' [',target,']\nrankRatio = ',round(SAMs$rankRatio[p],2),
                           '\nHGp = ',formatC(SAMs$HG_pval[p], format = "e", digits = 2),
                           '\nEMPp =',formatC(SAMs$EMP_pval[p], format = "e", digits = 2),sep=''),
              cex.main = 0.8)
         
         abline(h=log(SAMs$max_conc[p]),col='gray',lty=2)
         abline(h=SAMs$medLn50_mutCLs[p],col='pink',lty=1)
         
         plot(0,0,frame.plot = FALSE,axes = FALSE,col=NA,xlab='',ylab='')
         legend('left',pch=16,col=c('red',"#0000FF6E"),legend=c(variant,'others'),
                bg = 'white',title=paste(target,'status'),cex=0.8)
         dev.off()
         
       })
       }
       
       SCREENdata<-list(screenInfo=additionalInfos,
                        lnIC50=lnIC50,
                        Zscores=zscores,
                        MutantCellLines=MutantCellLines)
       
       save(SCREENdata,
            file=paste(resultPath,'/_DR_plots/',ctiss,'/',target,' _ ',paste(gsub("\\?|\\*|!|>|\\|", "", variant),collapse='AND'),'_screenRes.RData',sep=''))
     }
    
    
     }
    })











