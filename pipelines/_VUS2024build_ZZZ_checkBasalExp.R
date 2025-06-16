
tract<-read.table(paste(pathdata,'/raw/Tractability_pacini_et_al_2024.txt',sep=''),sep="\t",
                  header=TRUE,row.names = 1)

load(paste(resultPath,'/_allHits.RData',sep=''))


NwithExpdata<-length(which(!is.na(allHits$avgBasalEXP_fpkm_in_ps_cl)))

print(paste(NwithExpdata,' hits have all positive cell lines with availble gene expression data',sep=''))

cexRanges<-5*(allHits$percBasalEXP_of_ps_cl/100)

ii<-which(allHits$avgBasalEXP_fpkm_in_ps_cl>5 & allHits$percBasalEXP_of_ps_cl>50 &
          allHits$medFitEff< -0.8 & allHits$inTOgen_driver_for=='none' & 
            (tract[allHits$GENE,1]<=4 | tract[allHits$GENE,2]<=4 | tract[allHits$GENE,3]<=4))

color_pal <- colorRampPalette(c("blue","red"))
colors <- color_pal(100)

plot(allHits$medFitEff[ii],log10(allHits$avgBasalEXP_fpkm_in_ps_cl[ii]+1),frame.plot=FALSE,
     xlab='DAM-bearing gene med-fitness effect',
     ylab='DAM-bearing gene log10 [Avg Expression(FPKM) + 1]',pch=16,cex=1.2,
     col=colors[round(allHits$percBasalEXP_of_ps_cl[ii])])

par(xpd=TRUE)
identify(allHits$medFitEff[ii],log10(allHits$avgBasalEXP_fpkm_in_ps_cl[ii]+1),
         paste(allHits$ctype[ii],'\n',allHits$GENE[ii],allHits$var[ii]),cex=0.5)




# allIndex<-do.call('rbind',lapply(1:nrow(allHits),function(x){
#   pscl<-unlist(str_split(allHits$ps_cl[x],', '))
#   g<-rep(allHits$GENE[x],length(pscl))
#   cbind(g,pscl)
# }))
# 
# allIndex<-allIndex[is.element(allIndex[,1],rownames(basal_exp)) & 
#                      is.element(allIndex[,2],colnames(basal_exp)),]
# 
# DAMSexpression<-basal_exp[allIndex]
# OtherExpression<-basal_exp[setdiff(as.vector(matrix(1:length(basal_exp),
#                                           nrow = nrow(basal_exp),
#                                           ncol = ncol(basal_exp))),as.vector(allIndex))]
# 
# boxplot(log10(OtherExpression+1),log10(DAMSexpression+1))
# 
# hist(log10(OtherExpression+1),100,xlim=c(0,6))
# hist(log10(DAMSexpression+1),100,xlim=c(0,6))


