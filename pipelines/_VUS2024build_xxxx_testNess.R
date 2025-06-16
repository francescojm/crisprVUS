load(paste(pathdata,"/Robj/Sanger_Broad_higQ_scaled_depFC.RData", sep=""))

ness_genes<-colSums(scaled_depFC < -0.5)
ndep_cls<-rowSums(scaled_depFC < -0.5)

ndep_cls<-ndep_cls[ndep_cls>0]



growthProperties<-CMP_annot[match(colnames(scaled_depFC),CMP_annot$model_name),'growth_properties']
growthProperties<-growthProperties$growth_properties

growthProperties<-CMP_annot[match(colnames(scaled_depFC),CMP_annot$model_name),'growth_properties']
growthProperties<-growthProperties$growth_properties

modelTreatment<-CMP_annot[match(colnames(scaled_depFC),CMP_annot$model_name),'model_treatment']
modelTreatment<-modelTreatment$model_treatment


growthProperties<-growthProperties$growth_properties

meth_values<-read.table('/Volumes/Macintosh HD/Users/francesco.iorio/Downloads/F2_METH_CELL_DATA.txt',header=TRUE)
idMatch<-read.csv('/Volumes/Macintosh HD/Users/francesco.iorio/Downloads/methSampleId_2_cosmicIds.csv',header=TRUE)

idtomatch<-paste('X',idMatch$Sentrix_ID,'_',idMatch$Sentrix_Position,sep='')
colnames(meth_values)<-idMatch$cosmic_id[match(colnames(meth_values),idtomatch)]
colnames(meth_values)<-CMP_annot$model_name[match(colnames(meth_values),CMP_annot$COSMIC_ID)]

subdep<-scaled_depFC[,ii]
submeth<-meth_values[,ii]

subness_genes<-colSums(subdep < -0.5)
subndep_cls<-rowSums(subdep < -0.5)
submet_burd<-colMeans(submeth)

gtotest<-setdiff(names(which(subndep_cls>2 & subndep_cls<800)),ADaM)

res<-do.call('rbind',lapply(gtotest,function(x){
  print(x)
  p<-t.test(submet_burd~(subdep[x,]< -0.5))$p.value
  ES<-
    (mean(submet_burd[which(subdep[x,]< -0.5)])-
       mean(submet_burd[which(subdep[x,]>= -0.5)]))/
    sd(submet_burd)
  ret<-c(ES,p)
}))

rownames(res)<-gtotest
colnames(res)<-c('es','p')

plot(res[,1],-log10(res[,2]))
par(xpd=TRUE)
identify(res[,1],-log10(res[,2]),rownames(res))




tissues<-CMP_annot[match(colnames(scaled_depFC),CMP_annot$model_name),'tissue']

boxplot(mutBurden$mutational_burden~as.factor(tissues$tissue),horizontal = TRUE,las=2)
boxplot(ness_genes~as.factor(growthProperties))



gtotest<-setdiff(names(which(ndep_cls>2 & ndep_cls<800)),ADaM)



res<-do.call('rbind',
  lapply(gtotest,function(x){
    p<-t.test(ness_genes~(scaled_depFC[x,]< -0.5))$p.value
    ES<-
      (mean(ness_genes[which(scaled_depFC[x,]< -0.5)])-
         mean(ness_genes[which(scaled_depFC[x,]>= -0.5)]))/
      sd(ness_genes)
    ret<-c(ES,p)
    }))
rownames(res)<-gtotest
colnames(res)<-c('es','p')

plot(res[,1],-log10(res[,2]))
identify(res[,1],-log10(res[,2]),labels=rownames(res))



res<-do.call('rbind',
             lapply(gtotest,function(x){
               p<-t.test(scaled_depFC[x,setdiff(1:ncol(scaled_depFC),c(68,683,692))],
                         scaled_depFC[x,c(68,683,692)])$p.value
               ES<-
                 (mean(scaled_depFC[x,setdiff(1:ncol(scaled_depFC),c(68,683,692))])-
                    mean(scaled_depFC[x,c(68,683,692)]))/
                 sd(scaled_depFC[x,])
               ret<-c(ES,p)
             }))
rownames(res)<-gtotest
colnames(res)<-c('es','p')



68 683 692

