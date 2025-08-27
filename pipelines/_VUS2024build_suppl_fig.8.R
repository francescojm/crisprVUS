####setting paths
home<-'/Users/francesco.iorio/Dropbox/CODING/vus/VUS 2025/'
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250808_bugFixed_and_RR_th.1.71_wr/'

load(paste(resultPath,'_allDAMs.RData',sep=''))
load(paste(resultPath,'_incl_cl_annot.RData',sep=''))

ncell<-nrow(incl_cl_annot)

ctiss<-unique(allDAMs$ctype)

ucls<-unique(unlist(str_split(allDAMs$ps_cl,', ')))

nDAMsPerCCL<-unlist(lapply(ucls,function(U){
  length(grep(U,allDAMs$ps_cl))
}))

names(nDAMsPerCCL)<-ucls

print(paste('DAMs are observed in ',length(ucls),' cell lines (',100*length(ucls)/ncell,'%)',sep=''))

pdf(paste(resultPath,'_figures_source/n.DAMs across n. CCLs harbouring them.pdf',sep=''),6,6)
hist(RES$NacrossCLs,100,main = 'n. DAMs across n. of CCLs in which they are observed',xlab='n. DAMs',ylab='n. CCLs')
dev.off()

print(paste(length(which(nDAMsPerCCL<4))/length(nDAMsPerCCL),' of CCLs with at least one DAMs no more than 3 DAMs'))

utiss<-unique(allDAMs$ctype)

nDAMsAcrossCLS<-lapply(utiss,function(ut){
  currentCLS<-incl_cl_annot$model_id[which(incl_cl_annot$cancer_type==ut)]
  nocc<-rep(0,length(currentCLS))
  names(nocc)<-currentCLS
  ii<-intersect(currentCLS,RES$uniqueCLs)
  nocc[ii]<-RES$NacrossCLs[match(ii,RES$uniqueCLs)]
  nocc<-sort(nocc,decreasing=TRUE)
  return(nocc)
})

names(nDAMsAcrossCLS)<-utiss

oo<-order(unlist(lapply(nDAMsAcrossCLS,'sum')),decreasing=TRUE)

nDAMsAcrossCLS<-nDAMsAcrossCLS[oo]


# Pad with zeros to same length
my_list<-nDAMsAcrossCLS
max_len <- max(lengths(my_list))
mat <- sapply(my_list, function(x) {
  c(x, rep(0, max_len - length(x)))
})

# Convert to matrix
mat <- as.matrix(mat)

# Stacked barplot

pdf(paste(resultPath,'_figures_source/n.DAMs_across_cell_lines_in_each_ctype.pdf',sep=''),15,8)
par(mar=c(17,4,2,2))
tt<-barplot(mat,
        beside = FALSE,
        col = c("skyblue", "salmon", "lightgreen", "orange"),
        border=FALSE,las=2,ylab='n. DAMs across cell lines',ylim=c(0,280))

text(tt,unlist(lapply(nDAMsAcrossCLS,sum))+5,nCellLinesWithDAMs)
dev.off()

nCellLinesWithDAMs<-unlist(lapply(nDAMsAcrossCLS,function(x){length(which(x>0))}))

nclsAcrossCtype<-unlist(lapply(names(nDAMsAcrossCLS),function(dd){
  length(which(incl_cl_annot$cancer_type==dd))
}))

pdf(paste(resultPath,'_figures_source/perc_cell_lines_with_DAMs_in_each_ctype.pdf',sep=''),15,8)
par(mar=c(17,4,2,2))
tt<-barplot(100*nCellLinesWithDAMs/nclsAcrossCtype,las=2,ylim=c(0,100),border=FALSE,ylab='% of cell lines with DAMs')
  
text(tt,100*nCellLinesWithDAMs/nclsAcrossCtype+2,nclsAcrossCtype)

dev.off()

