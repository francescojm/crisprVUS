load("~/VUS/VUS/results/20220208/misclass_expr_rank2.RData")

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  max_cor<-which.max(misclass_result_all[i,c(7:39)])
  if(colnames(misclass_result_all[i,c(7:39)])[max_cor]==othertype & misclass_result_all[i,-c(1:39)][max_cor]<0.05){
  ind_sel<-c(ind_sel, i)
  }
}

misclass_expr_sel<-misclass_result_all[ind_sel,]


load("~/VUS/VUS/results/20220208/misclass_ess_rank.RData")

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  max_cor<-which.max(misclass_result_all[i,c(7:39)])
  if(colnames(misclass_result_all[i,c(7:39)])[max_cor]==othertype & misclass_result_all[i,-c(1:39)][max_cor]<0.05){
    ind_sel<-c(ind_sel, i)
  }
}

misclass_ess_sel<-misclass_result_all[ind_sel,]


load("~/VUS/VUS/results/20220208/misclass_drugs_rank.RData")

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  max_cor<-which.max(misclass_result_all[i,c(7:39)])
  if(colnames(misclass_result_all[i,c(7:39)])[max_cor]==othertype & misclass_result_all[i,-c(1:39)][max_cor]<0.05){
    ind_sel<-c(ind_sel, i)
  }
}

misclass_drugs_sel<-misclass_result_all[ind_sel,]

load("~/VUS/VUS/results/20220208/misclass_prot_rank.RData")

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  max_cor<-which.max(misclass_result_all[i,c(7:39)])
  if(colnames(misclass_result_all[i,c(7:39)])[max_cor]==othertype & misclass_result_all[i,-c(1:39)][max_cor]<0.05){
    ind_sel<-c(ind_sel, i)
  }
}

misclass_prot_sel<-misclass_result_all[ind_sel,]



Reduce(intersect, list(misclass_expr_sel$line, misclass_drugs_sel$line, misclass_ess_sel$line, misclass_prot_sel$line))

table(c(unique(misclass_expr_sel$line), unique(misclass_drugs_sel$line), unique(misclass_ess_sel$line), unique(misclass_prot_sel$line)))

#OV-7 in 3

################################
#######plots
#################################

load("~/VUS/VUS/results/20220208/misclass_expr_rank2.RData")
misclass_expr<-misclass_result_all


load("~/VUS/VUS/results/20220208/misclass_ess_rank.RData")
misclass_ess<-misclass_result_all


load("~/VUS/VUS/results/20220208/misclass_drugs_rank.RData")
misclass_drugs<-misclass_result_all

load("~/VUS/VUS/results/20220208/misclass_prot_rank.RData")
misclass_prot<-misclass_result_all


corrs_ov7<-rbind(misclass_expr[which(misclass_expr$line=="OV-7" & misclass_expr$other_type=="Thyroid Gland Carcinoma") ,c(7:39)],
                 misclass_drugs[which(misclass_drugs$line=="OV-7" & misclass_drugs$other_type=="Thyroid Gland Carcinoma"),c(7:39)],
                 misclass_ess[which(misclass_ess$line=="OV-7" & misclass_ess$other_type=="Thyroid Gland Carcinoma") ,c(7:39)],
                 misclass_prot[which(misclass_prot$line=="OV-7" & misclass_prot$other_type=="Thyroid Gland Carcinoma") ,c(7:39)])

rownames(corrs_ov7)<-c("Transcriptome", "Drug screens", "Essentiality", "Proteome")
order_ov7<-order(rowMeans(apply(-corrs_ov7,1, rank)))

pdf("OV7_misclass.pdf",10,5)
pheatmap(cellwidth = 15, cellheight = 15, corrs_ov7[,order_ov7],  cluster_rows = F, cluster_cols=F, scale="row")
dev.off()

###########
## Look for cell lines consistently not at the top
##########

load("~/VUS/VUS/results/20220208/misclass_expr_rank2.RData")

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  p_othertype<-misclass_result_all[i,paste(othertype, "pval")]
  cor_othertype<-misclass_result_all[i,othertype]
  if(cor_othertype>cor_tum & !is.na(p_othertype) & p_othertype<0.05){
    ind_sel<-c(ind_sel, i)
  }
}

misclass_expr_sel<-misclass_result_all[ind_sel,]


load("~/VUS/VUS/results/20220208/misclass_ess_rank.RData")

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  p_othertype<-misclass_result_all[i,paste(othertype, "pval")]
  cor_othertype<-misclass_result_all[i,othertype]
  if(cor_othertype>cor_tum & !is.na(p_othertype) & p_othertype<0.05){
    ind_sel<-c(ind_sel, i)
  }
}

misclass_ess_sel<-misclass_result_all[ind_sel,]


load("~/VUS/VUS/results/20220208/misclass_drugs_rank.RData")

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  p_othertype<-misclass_result_all[i,paste(othertype, "pval")]
  cor_othertype<-misclass_result_all[i,othertype]
  if(cor_othertype>cor_tum & !is.na(p_othertype) & p_othertype<0.05){
    ind_sel<-c(ind_sel, i)
  }
}

misclass_drugs_sel<-misclass_result_all[ind_sel,]

load("~/VUS/VUS/results/20220208/misclass_prot_rank.RData")

misclass_result_all<-misclass_result_all[-1,]
ind_sel<-c()
for(i in 1:nrow(misclass_result_all)){
  tum<-misclass_result_all$tumor[i]
  othertype<-misclass_result_all$other_type[i]
  cor_tum<-misclass_result_all[i,tum]
  p_othertype<-misclass_result_all[i,paste(othertype, "pval")]
  cor_othertype<-misclass_result_all[i,othertype]
  if(cor_othertype>cor_tum & !is.na(p_othertype) & p_othertype<0.05){
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

which(table(c(unique(misclass_check_all[[1]]), unique(misclass_check_all[[2]]), unique(misclass_check_all[[3]]), unique(misclass_check_all[[4]])))==4)

corrs_tknu<-rbind(misclass_expr[which(misclass_expr$line=="TYK-nu" & misclass_expr$other_type=="Thyroid Gland Carcinoma") ,c(7:39)],
                 misclass_drugs[which(misclass_drugs$line=="TYK-nu" & misclass_drugs$other_type=="Thyroid Gland Carcinoma"),c(7:39)],
                 misclass_ess[which(misclass_ess$line=="TYK-nu" & misclass_ess$other_type=="Thyroid Gland Carcinoma") ,c(7:39)],
                 misclass_prot[which(misclass_prot$line=="TYK-nu" & misclass_prot$other_type=="Thyroid Gland Carcinoma") ,c(7:39)])

rownames(corrs_tknu)<-c("Transcriptome", "Drug screens", "Essentiality", "Proteome")
order_tknu<-order(rowMeans(apply(-corrs_tknu,1, rank)))

pdf("TYK-nu_misclass.pdf",10,5)
pheatmap(cellwidth = 15, cellheight = 15, corrs_tknu[,order_tknu],  cluster_rows = F, cluster_cols=F, scale="row")
dev.off()


corrs_HCC1806<-rbind(unique(misclass_expr[which(misclass_expr$line=="HCC1806" & misclass_expr$other_type=="Esophageal Squamous Cell Carcinoma") ,c(7:39)]),
                     unique(misclass_drugs[which(misclass_drugs$line=="HCC1806" & misclass_drugs$other_type=="Esophageal Squamous Cell Carcinoma"),c(7:39)]),
                            unique(misclass_ess[which(misclass_ess$line=="HCC1806" & misclass_ess$other_type=="Esophageal Squamous Cell Carcinoma") ,c(7:39)]),
                                   unique(misclass_prot[which(misclass_prot$line=="HCC1806" & misclass_prot$other_type=="Esophageal Squamous Cell Carcinoma") ,c(7:39)]))

rownames(corrs_HCC1806)<-c("Transcriptome", "Drug screens", "Essentiality", "Proteome")
order_HCC1806<-order(rowMeans(apply(-corrs_HCC1806,1, rank)))

pdf("HCC1806_misclass.pdf",10,5)
pheatmap(cellwidth = 15, cellheight = 15, corrs_HCC1806[,order_HCC1806],  cluster_rows = F, cluster_cols=F, scale="row")
dev.off()

