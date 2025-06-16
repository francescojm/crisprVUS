library(readr)
library(stringr)

broadVar<-read.csv(paste(pathdata,'/raw/vcf_muts_broad_wes_prod.csv',sep=""),
                   stringsAsFactors = FALSE,header=TRUE)
sangerVar<-read.csv(paste(pathdata,'/raw/vcf_muts_sanger_wes_prod.csv',sep=""),
                   stringsAsFactors = FALSE,header=TRUE)

broadVar<-broadVar[,colnames(sangerVar)]

allVar<-rbind(broadVar,sangerVar)


gene_annot<-read_csv(paste(pathdata,'/gene_identifiers_20191101.csv', sep=""))

load(file = paste(resultPath,'/_incl_cl_annot.RData', sep=""))

ctypes<-sort(unique(incl_cl_annot$cancer_type))

all_HITS<-NULL

for (ct in ctypes){
  print(ct)
  load(paste(resultPath,"/",ct,'_results.RData',sep=''))
  ids<-which(RESTOT$rank_ratio<1.6 & RESTOT$medFitEff < -.5)
  RESTOT<-RESTOT[ids,]
  
  current_ct_res<-do.call(rbind,lapply(1:length(ids),function(ii){
    currentLine<-RESTOT[ii,]
    GENE<-currentLine$GENE
    variantsSet<-str_trim(unlist(str_split(currentLine$var,'[|]')))
    GENEinfo<-gene_annot[gene_annot$hgnc_symbol==GENE,c(1,3:6,8)]
    
    RES<-do.call(rbind,lapply(variantsSet,function(cv){
      allVar[which(allVar$gene_id == as.character(GENEinfo[1]) & allVar$protein==cv)[1],c(1,2,3,6,12,13)]
    }))
    
    RES<-cbind(rep(as.character(GENEinfo[2]),length(variantsSet)),RES)
    RES<-cbind(rep(as.character(GENEinfo[3]),length(variantsSet)),RES)
    RES<-cbind(rep(as.character(GENEinfo[4]),length(variantsSet)),RES)
    RES<-cbind(rep(as.character(GENEinfo[5]),length(variantsSet)),RES)
    
    colnames(RES)[1:4]<-colnames(GENEinfo)[c(5,4,3,2)]
    abbr<-paste(RES[,1],str_sub(RES[,9],3,str_length(RES[,9])),sep='-')
    RES<-cbind(abbr,RES)
  }))
  
  current_ct_res<-cbind(rep(ct,nrow(current_ct_res)),current_ct_res)
  colnames(current_ct_res)[1]<-'analysis'
  all_HITS<-rbind(all_HITS,current_ct_res)
}

vs_to_assess<-all_HITS[match(unique(all_HITS$abbr),all_HITS$abbr),]
vs_to_assess<-vs_to_assess[,-1]

vs_to_assess<-vs_to_assess[order(vs_to_assess$abbr),]

vs_to_assess<-vs_to_assess[vs_to_assess$protein!='p.?',]

save(vs_to_assess,file='../results/20210718/_vs_to_asses.RData')
write.table(vs_to_assess,quote=FALSE,sep='\t',row.names = FALSE,file='../results/20210718/_vs_to_asses.tsv')




