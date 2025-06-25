pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250221/'

load(paste(resultPath, "misclass_expr_rank.RData", sep=""))

misclass_result_all<-misclass_result_all[-1,]
all_expr<-matrix(nrow=1,ncol=4)
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  cor_other<-misclass_result_all[i,othertype]
  p_other<-misclass_result_all[i,paste(othertype, "pval")]
  line<-misclass_result_all[i,"line"]
  all_expr<-rbind(all_expr, c(line, cor_tum, cor_other, p_other))
}
all_expr<-as.data.frame(all_expr)
colnames(all_expr)<-c("line", "cor_tum", "cor_other", "p_other")
all_expr$p_other[all_expr$p_other==0]<-0.01
all_expr$line[all_expr$line!="TYK-nu"]<-NA
library(ggrepel)
p1<-ggplot(all_expr, aes(x=as.numeric(cor_tum), y=as.numeric(cor_other), label=line, colour=-log10(as.numeric(p_other))))+geom_point()+theme_classic()+
  geom_abline(slope = 1)+geom_label_repel(max.overlaps = 10)+xlab("Correlation with tumor of origin")+
  ylab("Correlation with alternate cancer type")+ guides(colour=guide_legend(title="-log10(p-value)"))

p1_hist<-ggplot(all_expr, aes(x=as.numeric(cor_other)-as.numeric(cor_tum), label=line, colour=-log10(as.numeric(p_other))))+
  geom_histogram()+theme_classic()+geom_vline(xintercept = 0, colour="red")+xlab("Correlation difference")

print(paste("Number of tested combinations", nrow(misclass_result_all)))
print(paste("Number of positive tests", sum((as.numeric(all_expr$cor_other)-as.numeric(all_expr$cor_tum))>0, na.rm=T)))

load(paste(resultPath, "misclass_ess_rank.RData", sep=""))

misclass_result_all<-misclass_result_all[-1,]
all_expr<-matrix(nrow=1,ncol=4)
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  cor_other<-misclass_result_all[i,othertype]
  p_other<-misclass_result_all[i,paste(othertype, "pval")]
  line<-misclass_result_all[i,"line"]
  all_expr<-rbind(all_expr, c(line, cor_tum, cor_other, p_other))
}
all_expr<-as.data.frame(all_expr)
colnames(all_expr)<-c("line", "cor_tum", "cor_other", "p_other")
all_expr$p_other[all_expr$p_other==0]<-0.01
#all_expr$line[all_expr$line!="TYK-nu"]<-NA
library(ggrepel)
p2<-ggplot(all_expr, aes(x=as.numeric(cor_tum), y=as.numeric(cor_other), label=line, colour=-log10(as.numeric(p_other))))+geom_point()+theme_classic()+
  geom_abline(slope = 1)+geom_label_repel(max.overlaps = 10)+xlab("Correlation with tumor of origin")+
  ylab("Correlation with alternate cancer type")+ guides(colour=guide_legend(title="-log10(p-value)"))

p2_hist<-ggplot(all_expr, aes(x=as.numeric(cor_other)-as.numeric(cor_tum), label=line, colour=-log10(as.numeric(p_other))))+
  geom_histogram()+theme_classic()+geom_vline(xintercept = 0, colour="red")+xlab("Correlation difference")

print(paste("Number of tested combinations", nrow(misclass_result_all)))
print(paste("Number of positive tests", sum((as.numeric(all_expr$cor_other)-as.numeric(all_expr$cor_tum))>0, na.rm=T)))

load(paste(resultPath, "misclass_drugs_rank.RData", sep=""))

misclass_result_all<-misclass_result_all[-1,]
all_expr<-matrix(nrow=1,ncol=4)
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  cor_other<-misclass_result_all[i,othertype]
  p_other<-misclass_result_all[i,paste(othertype, "pval")]
  line<-misclass_result_all[i,"line"]
  all_expr<-rbind(all_expr, c(line, cor_tum, cor_other, p_other))
}
all_expr<-as.data.frame(all_expr)
colnames(all_expr)<-c("line", "cor_tum", "cor_other", "p_other")
all_expr$p_other[all_expr$p_other==0]<-0.01
all_expr$line[all_expr$line!="TYK-nu"]<-NA

library(ggrepel)
p3<-ggplot(all_expr, aes(x=as.numeric(cor_tum), y=as.numeric(cor_other), label=line, colour=-log10(as.numeric(p_other))))+geom_point()+theme_classic()+
  geom_abline(slope = 1)+geom_label_repel(max.overlaps = 10)+xlab("Correlation with tumor of origin")+
  ylab("Correlation with alternate cancer type")+ guides(colour=guide_legend(title="-log10(p-value)"))

p3_hist<-ggplot(all_expr, aes(x=as.numeric(cor_other)-as.numeric(cor_tum), label=line, colour=-log10(as.numeric(p_other))))+
  geom_histogram()+theme_classic()+geom_vline(xintercept = 0, colour="red")+xlab("Correlation difference")

print(paste("Number of tested combinations", nrow(misclass_result_all)))
print(paste("Number of positive tests", sum((as.numeric(all_expr$cor_other)-as.numeric(all_expr$cor_tum))>0, na.rm=T)))

load(paste(resultPath, "misclass_prot_rank.RData", sep=""))

misclass_result_all<-misclass_result_all[-1,]
all_expr<-matrix(nrow=1,ncol=4)
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  cor_other<-misclass_result_all[i,othertype]
  p_other<-misclass_result_all[i,paste(othertype, "pval")]
  line<-misclass_result_all[i,"line"]
  all_expr<-rbind(all_expr, c(line, cor_tum, cor_other, p_other))
}
all_expr<-as.data.frame(all_expr)
colnames(all_expr)<-c("line", "cor_tum", "cor_other", "p_other")
all_expr$p_other[all_expr$p_other==0]<-0.01
#all_expr$line[all_expr$line!="TYK-nu"]<-NA
library(ggrepel)
p4<-ggplot(all_expr, aes(x=as.numeric(cor_tum), y=as.numeric(cor_other), label=line, colour=-log10(as.numeric(p_other))))+geom_point()+theme_classic()+
  geom_abline(slope = 1)+geom_label_repel(max.overlaps = 10)+xlab("Correlation with tumor of origin")+
  ylab("Correlation with alternate cancer type")+ guides(colour=guide_legend(title="-log10(p-value)"))

p4_hist<-ggplot(all_expr, aes(x=as.numeric(cor_other)-as.numeric(cor_tum), label=line, colour=-log10(as.numeric(p_other))))+
  geom_histogram()+theme_classic()+geom_vline(xintercept = 0, colour="red")+xlab("Correlation difference")

print(paste("Number of tested combinations", nrow(misclass_result_all)))
print(paste("Number of positive tests", sum((as.numeric(all_expr$cor_other)-as.numeric(all_expr$cor_tum))>0, na.rm=T)))

pdf(paste(resultPath, "misclass_expr_scatter.pdf", sep=""),5,5)
print(p1)
dev.off()

pdf(paste(resultPath, "misclass_ess_scatter.pdf", sep=""),5,5)
print(p2)
dev.off()

pdf(paste(resultPath, "misclass_drugs_scatter.pdf", sep=""),5,5)
print(p3)
dev.off()

pdf(paste(resultPath, "misclass_prot_scatter.pdf", sep=""),5,5)
print(p4)
dev.off()

pdf(paste(resultPath, "misclass_expr_hist.pdf", sep=""),5,5)
print(p1_hist)
dev.off()

pdf(paste(resultPath, "misclass_ess_hist.pdf", sep=""),5,5)
print(p2_hist)
dev.off()

pdf(paste(resultPath, "misclass_drugs_hist.pdf", sep=""),5,5)
print(p3_hist)
dev.off()

pdf(paste(resultPath, "misclass_prot_hist.pdf", sep=""),5,5)
print(p4_hist)
dev.off()