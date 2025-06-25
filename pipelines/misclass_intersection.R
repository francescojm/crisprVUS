library(pheatmap)
####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250221/'

load(paste(resultPath, "misclass_expr_rank.RData", sep=""))

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  max_cor<-which.max(misclass_result_all[i,c(7:42)])
  if(colnames(misclass_result_all[i,c(7:42)])[max_cor]==othertype & misclass_result_all[i,-c(1:42)][max_cor]<0.05){
  ind_sel<-c(ind_sel, i)
  }
}

misclass_expr_sel<-misclass_result_all[ind_sel,]

load(paste(resultPath, "misclass_ess_rank.RData", sep=""))

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  max_cor<-which.max(misclass_result_all[i,c(7:42)])
  if(colnames(misclass_result_all[i,c(7:42)])[max_cor]==othertype & misclass_result_all[i,-c(7:42)][max_cor]<0.05){
    ind_sel<-c(ind_sel, i)
  }
}

misclass_ess_sel<-misclass_result_all[ind_sel,]

load(paste(resultPath, "misclass_drugs_rank.RData", sep=""))

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  max_cor<-which.max(misclass_result_all[i,c(7:42)])
  if(colnames(misclass_result_all[i,c(7:42)])[max_cor]==othertype & misclass_result_all[i,-c(1:42)][max_cor]<0.05){
    ind_sel<-c(ind_sel, i)
  }
}

misclass_drugs_sel<-misclass_result_all[ind_sel,]

load(paste(resultPath, "misclass_prot_rank.RData", sep=""))

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  max_cor<-which.max(misclass_result_all[i,c(7:42)])
  if(colnames(misclass_result_all[i,c(7:42)])[max_cor]==othertype & misclass_result_all[i,-c(1:42)][max_cor]<0.05){
    ind_sel<-c(ind_sel, i)
  }
}

misclass_prot_sel<-misclass_result_all[ind_sel,]



Reduce(intersect, list(misclass_expr_sel$line, misclass_drugs_sel$line, misclass_ess_sel$line, misclass_prot_sel$line))

table(c(unique(misclass_expr_sel$line), unique(misclass_drugs_sel$line), unique(misclass_ess_sel$line), unique(misclass_prot_sel$line)))

#SIDM00460 in 3

################################
#######plots
#################################

load(paste(resultPath, "misclass_expr_rank.RData", sep=""))
misclass_expr<-misclass_result_all


load(paste(resultPath, "misclass_ess_rank.RData", sep=""))
misclass_ess<-misclass_result_all


load(paste(resultPath, "misclass_drugs_rank.RData", sep=""))
misclass_drugs<-misclass_result_all

load(paste(resultPath, "misclass_prot_rank.RData", sep=""))
misclass_prot<-misclass_result_all


###########
## Look for cell lines consistently not at the top
##########

load(paste(resultPath, "misclass_expr_rank.RData", sep=""))

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  p_othertype<-misclass_result_all[i,paste(othertype, "pval")]
  cor_othertype<-misclass_result_all[i,othertype]
  if(!is.na(cor_othertype)){
    if(cor_othertype>cor_tum & !is.na(p_othertype) & p_othertype<0.25){
    ind_sel<-c(ind_sel, i)
    }
  }
}

misclass_expr_sel<-misclass_result_all[ind_sel,]

load(paste(resultPath, "misclass_ess_rank.RData", sep=""))

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  p_othertype<-misclass_result_all[i,paste(othertype, "pval")]
  cor_othertype<-misclass_result_all[i,othertype]
  if(cor_othertype>cor_tum & !is.na(p_othertype) & p_othertype<0.25){
    ind_sel<-c(ind_sel, i)
  }
}

misclass_ess_sel<-misclass_result_all[ind_sel,]

load(paste(resultPath, "misclass_drugs_rank.RData", sep=""))

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  p_othertype<-misclass_result_all[i,paste(othertype, "pval")]
  cor_othertype<-misclass_result_all[i,othertype]
  if(cor_othertype>cor_tum & !is.na(p_othertype) & p_othertype<0.25){
    ind_sel<-c(ind_sel, i)
  }
}

misclass_drugs_sel<-misclass_result_all[ind_sel,]

load(paste(resultPath, "misclass_prot_rank.RData", sep=""))

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  p_othertype<-misclass_result_all[i,paste(othertype, "pval")]
  cor_othertype<-misclass_result_all[i,othertype]
  if(cor_othertype>cor_tum & !is.na(p_othertype) & p_othertype<0.25){
    ind_sel<-c(ind_sel, i)
  }
}

misclass_prot_sel<-misclass_result_all[ind_sel,]



misclass_check_all<-matrix(nrow=length(unique(c(c(misclass_expr$line), c(misclass_drugs$line), 
                                             c(misclass_ess$line), c(misclass_prot$line)))), ncol=4)
rownames(misclass_check_all)<-unique(c(c(misclass_expr$line), c(misclass_drugs$line), 
                                       c(misclass_ess$line), c(misclass_prot$line)))
all_lines<-unique(c(c(misclass_expr$line), c(misclass_drugs$line), 
                    c(misclass_ess$line), c(misclass_prot$line)))
misclass_check_all<-list()
misclass_check_all[[1]]<-NA
misclass_check_all[[2]]<-NA
misclass_check_all[[3]]<-NA
misclass_check_all[[4]]<-NA

for(line in all_lines){
  if(!is.na(match(line, misclass_expr_sel$line))){
  misclass_check_all[[1]]<-c(misclass_check_all[[1]], paste(line, misclass_expr_sel$other_type[which(misclass_expr_sel$line==line)]))
  }
  if(!is.na(match(line, misclass_drugs_sel$line))){
    misclass_check_all[[2]]<-c(misclass_check_all[[2]], paste(line, misclass_drugs_sel$other_type[which(misclass_drugs_sel$line==line)]))
  }
  if(!is.na(match(line, misclass_ess_sel$line))){
    misclass_check_all[[3]]<-c(misclass_check_all[[3]], paste(line, misclass_ess_sel$other_type[which(misclass_ess_sel$line==line)]))
  }
  if(!is.na(match(line, misclass_prot_sel$line))){
    misclass_check_all[[4]]<-c(misclass_check_all[[4]], paste(line, misclass_prot_sel$other_type[which(misclass_prot_sel$line==line)]))
  }
}

which(table(c(unique(misclass_check_all[[1]]), unique(misclass_check_all[[2]]), unique(misclass_check_all[[3]]), unique(misclass_check_all[[4]])))>3)

corrs_SIDM00227<-rbind(misclass_expr[which(misclass_expr$line=="SIDM00227" & misclass_expr$other_type=="Pancreatic Carcinoma") ,c(7:42)],
                 misclass_drugs[which(misclass_drugs$line=="SIDM00227" & misclass_drugs$other_type=="Pancreatic Carcinoma"),c(7:42)],
                 misclass_ess[which(misclass_ess$line=="SIDM00227" & misclass_ess$other_type=="Pancreatic Carcinoma") ,c(7:42)],
                 misclass_prot[which(misclass_prot$line=="SIDM00227" & misclass_prot$other_type=="Pancreatic Carcinoma") ,c(7:42)])

rownames(corrs_SIDM00227)<-c("Transcriptome", "Drug screens", "Essentiality", "Proteome")
order_SIDM00227<-order(rowMeans(apply(-corrs_SIDM00227,1, rank)))

pdf(paste(resultPath, "SIDM00227_misclass.pdf", sep="/"),10,5)
pheatmap(cellwidth = 15, cellheight = 15, corrs_SIDM00227[,order_SIDM00227],  cluster_rows = F, cluster_cols=F, scale="row")
dev.off()

cl<-"SIDM00460"
tum<-"B-Lymphoblastic Leukemia"
corrs_SIDM00460<-rbind(misclass_expr[which(misclass_expr$line==cl & misclass_expr$other_type==tum) ,c(7:42)],
                  misclass_drugs[which(misclass_drugs$line==cl & misclass_drugs$other_type==tum),c(7:42)],
                  misclass_ess[which(misclass_ess$line==cl & misclass_ess$other_type==tum) ,c(7:42)],
                  misclass_prot[which(misclass_prot$line==cl & misclass_prot$other_type==tum) ,c(7:42)])

rownames(corrs_SIDM00460)<-c("Transcriptome", "Drug screens", "Essentiality", "Proteome")
order_SIDM00460<-order(rowMeans(apply(-corrs_SIDM00460,1, rank)))

pdf(paste(resultPath, "SIDM00460_misclass.pdf", sep="/"),10,5)
pheatmap(cellwidth = 15, cellheight = 15, corrs_SIDM00460[,order_SIDM00460],  cluster_rows = F, cluster_cols=F, scale="row")
dev.off()


