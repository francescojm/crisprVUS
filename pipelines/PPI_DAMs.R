load("/Volumes/home/VUS/VUS/data/preComputed/PPI_mat.RData")
load("/Volumes/home/VUS/VUS/data/preComputed/PPI.RData")
write.csv(PPI, file="PPI_cytoscape.csv", quote=F)


library(pheatmap)
library(ggplot2)
library(openxlsx)
library(ReactomePA)
library(clusterProfiler)
library(tidyverse)

path_data<-"/data"
path_results<-"/results/20220208"
home<-"/Volumes/home/VUS/VUS"

setwd(paste(home, "/", path_results, sep=""))

##load dei dati
tissues<-gsub("_results_ext.RData", "", list.files(pattern="results_ext.RData"))

results<-list()
ind_iter<-0
for(ctiss in tissues){
  ind_iter<-ind_iter+1
  load(paste(ctiss, "_results_ext.RData", sep=""))
  results[[ind_iter]] <- RESTOT
}
names(results)<-tissues

##driver genes
inTOgen_drivers<-read.table(paste(home, '/data/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
driver_genes<-unique(inTOgen_drivers$SYMBOL)

gene_annot <- read_csv(paste(home, "/",path_data, "/gene_identifiers_20191101.csv", sep=""))
load(paste(home, "/", path_data,'/R/Sanger_Broad_higQ_bdep.RData', sep=""))
load(paste(home, "/", path_data,'/preComputed/ADaM.RData', sep=""))
load(paste(home, "/", path_data,'/preComputed/FiPer_outputs.RData', sep=""))

load(paste(home, "/", path_data, '/R/Sanger_Broad_higQ_scaled_depFC.RData', sep=""))

CMP_annot <- read_csv(paste(home, "/",path_data,"/model_list_20210611.csv", sep="")) # from https://cog.sanger.ac.uk/cmp/download/model_list_20210611.csv

## latest sanger/broad unreleased yet variants hg38
cl_variants <- read_csv(paste(home, "/", path_data,'/mutations_all_latest.csv', sep=""))
cl_variants <- cbind(cl_variants,gene_annot$hgnc_symbol[match(cl_variants$gene_id,gene_annot$gene_id)])
colnames(cl_variants)[ncol(cl_variants)]<-'gene_symbol_2019'

cl_variants <- cbind(cl_variants,CMP_annot$model_name[match(cl_variants$model_id,CMP_annot$model_id)])
colnames(cl_variants)[ncol(cl_variants)] <- 'model_name'


cl_variants <- cbind(cl_variants,CMP_annot$cancer_type[match(cl_variants$model_id,CMP_annot$model_id)])
colnames(cl_variants)[ncol(cl_variants)] <- 'cancer_type'

colnames(cl_variants)[12] <- 'protein_mutation'
CMP_annot<-CMP_annot[which(is.element(CMP_annot$model_name,colnames(scaled_depFC))),]


load("/Volumes/home/VUS/VUS/results/20220208/hits.RData")




probs1<-c()
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)
  gene<-results[[ctiss]][ind,"GENE"]
  cl<-results[[ctiss]][ind,"ps_cl"]
  
  for(i in 1:length(ind)){
    g<-gene[i]
    line<-cl[i]
    
    if(g %in% rownames(PPI_mat)){
      if(length(grep(",",line))==1){
        lines<-unlist(strsplit(line, ", "))
        for(line in lines){
          PPI_sel<-PPI_mat[g, ]
          first_neigh<-names(which(PPI_sel==1))
          first_neigh<-intersect(first_neigh, rownames(scaled_depFC))
          rand_neigh<-setdiff(colnames(PPI_mat)[sample(1:ncol(PPI_mat),length(first_neigh), replace=F)],first_neigh)
          rand_neigh<-intersect(rownames(scaled_depFC), rand_neigh)
          rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
          rand_lines<-setdiff(rand_lines, lines)
          if(length(rand_neigh)>0){
          rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
          for(neigh in 1:length(rand_neigh)){
            probs1<-c(probs1,sum(scaled_depFC[rand_neigh[neigh], rand_lines]<scaled_depFC[intersect(rownames(scaled_depFC), rand_neigh[neigh]), line])/length(intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])))
            
          }
          }
        }
      } else {
        PPI_sel<-PPI_mat[g, ]
        first_neigh<-names(which(PPI_sel==1))
        first_neigh<-intersect(first_neigh, rownames(scaled_depFC))
        rand_neigh<-setdiff(colnames(PPI_mat)[sample(1:ncol(PPI_mat),length(first_neigh), replace=F)],first_neigh)
        rand_neigh<-intersect(rownames(scaled_depFC), rand_neigh)
        rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
        rand_lines<-setdiff(rand_lines, line)
        if(length(rand_neigh)>0){
        for(neigh in 1:length(rand_neigh)){
          probs1<-c(probs1,sum(scaled_depFC[rand_neigh[neigh], rand_lines]<scaled_depFC[intersect(rownames(scaled_depFC), rand_neigh[neigh]), line])/length(intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])))
                   
        }
        }
      }
      
    }
  }
}

hist(probs1, breaks=100)
sum(probs1<0.05)/length(probs1)

probs2<-c()
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)
  gene<-results[[ctiss]][ind,"GENE"]
  cl<-results[[ctiss]][ind,"ps_cl"]
  
  for(i in 1:length(ind)){
    g<-gene[i]
    line<-cl[i]
    
    if(g %in% rownames(PPI_mat)){
      if(length(grep(",",line))==1){
        lines<-unlist(strsplit(line, ", "))
        for(line in lines){
          PPI_sel<-PPI_mat[g, ]
          first_neigh<-names(which(PPI_sel==1))
          first_neigh<-intersect(first_neigh, rownames(scaled_depFC))
          rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
          rand_lines<-setdiff(rand_lines, lines)
          if(length(first_neigh)>0){
            rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
            for(neigh in 1:length(first_neigh)){
              probs2<-c(probs2,sum(scaled_depFC[first_neigh[neigh], rand_lines]<scaled_depFC[first_neigh[neigh], line])/length(intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])))
              
            }
          }
        }
      } else {
        PPI_sel<-PPI_mat[g, ]
        first_neigh<-names(which(PPI_sel==1))
        first_neigh<-intersect(first_neigh, rownames(scaled_depFC))
        rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
        rand_lines<-setdiff(rand_lines, line)
        if(length(first_neigh)>0){
          for(neigh in 1:length(first_neigh)){
            probs2<-c(probs2,sum(scaled_depFC[first_neigh[neigh], rand_lines]<scaled_depFC[first_neigh[neigh], line])/length(intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])))
            
          }
        }
      }
      
    }
  }
}

sum(probs2<0.05)/length(probs2)

df<-data.frame(ratio=c(probs2, probs1), type=c(rep("DAMs", length(probs2)), rep("random genes", length(probs1))))
ggplot(df, aes(x=ratio, colour=type, group=type))+geom_density()
ks.test(probs1, probs2)

probs<-c()
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)
  gene<-results[[ctiss]][ind,"GENE"]
  cl<-results[[ctiss]][ind,"ps_cl"]
  
  for(i in 1:length(ind)){
    g<-gene[i]
    line<-cl[i]
    
    if(g %in% rownames(PPI_mat)){
      if(length(grep(",",line))==1){
        lines<-unlist(strsplit(line, ", "))
        for(line in lines){
          PPI_sel<-PPI_mat[g, ]
          first_neigh<-names(which(PPI_sel==1))
          rand_neigh<-setdiff(colnames(PPI_mat),first_neigh)
          first_neigh<-intersect(first_neigh, rownames(bdep))
          for(neigh in 1:length(first_neigh)){
            probs<-c(probs,sum(bdep[intersect(rownames(bdep), rand_neigh), line]>bdep[first_neigh, line][neigh])/length(intersect(rownames(bdep), rand_neigh)))
            
          }
        }
      } else {
        PPI_sel<-PPI_mat[g, ]
        first_neigh<-names(which(PPI_sel==1))
        rand_neigh<-setdiff(colnames(PPI_mat),first_neigh)
        first_neigh<-intersect(first_neigh, rownames(bdep))
        for(neigh in 1:length(first_neigh)){
          probs<-c(probs,sum(bdep[intersect(rownames(bdep), rand_neigh), line]>bdep[first_neigh, line][neigh])/length(intersect(rownames(bdep), rand_neigh)))
          
        }
      }
      
    }
  }
}

hist(probs, breaks=100)


#######Francesco
######################

actual<-c()
others<-c()
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)
  gene<-results[[ctiss]][ind,"GENE"]
  cl<-results[[ctiss]][ind,"ps_cl"]
  
  for(i in 1:length(ind)){
    g<-gene[i]
    line<-cl[i]
    
    if(g %in% rownames(PPI_mat)){
      if(length(grep(",",line))==1){
        lines<-unlist(strsplit(line, ", "))
        
        PPI_sel<-PPI_mat[g, ]
        first_neigh<-names(which(PPI_sel==1))
        first_neigh<-intersect(first_neigh, rownames(scaled_depFC))
        rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
        rand_lines<-setdiff(rand_lines, lines)
        if(length(first_neigh)>0){
          rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
          actual<-c(actual,unlist(scaled_depFC[first_neigh, lines]))
          others<-c(others,unlist(scaled_depFC[first_neigh, rand_lines]))
        }
        
      } else {
        PPI_sel<-PPI_mat[g, ]
        first_neigh<-names(which(PPI_sel==1))
        first_neigh<-intersect(first_neigh, rownames(scaled_depFC))
        rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
        rand_lines<-setdiff(rand_lines, line)
        if(length(first_neigh)>0){
          actual<-c(actual,unlist(scaled_depFC[first_neigh, line]))
          others<-c(others,unlist(scaled_depFC[first_neigh, rand_lines]))
          
        }
      }
      
    }
  }
}

df<-data.frame(logFC=c(actual, others), type=c(rep("DAM-bearing lines", length(actual)), rep("other lines", length(others))))
ggplot(df, aes(x=logFC, colour=type, group=type))+geom_density()
ks.test(actual, others, alternative = "greater")
ks.test(actual, others, alternative = "less")

###second neighbours

actual2<-c()
others2<-c()
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)
  gene<-results[[ctiss]][ind,"GENE"]
  cl<-results[[ctiss]][ind,"ps_cl"]
  
  for(i in 1:length(ind)){
    g<-gene[i]
    line<-cl[i]
    
    if(g %in% rownames(PPI_mat)){
      if(length(grep(",",line))==1){
        lines<-unlist(strsplit(line, ", "))
        
        PPI_sel<-PPI_mat[g, ]
        first_neigh<-names(which(PPI_sel==1))
        PPI_sel<-PPI_mat[first_neigh, ]
        first_neigh<-names(which(PPI_sel==1))
        first_neigh<-setdiff(intersect(first_neigh, rownames(scaled_depFC)),g)
        rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
        rand_lines<-setdiff(rand_lines, lines)
        if(length(first_neigh)>0){
          rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
          actual2<-c(actual2,unlist(scaled_depFC[first_neigh, lines]))
          others2<-c(others2,unlist(scaled_depFC[first_neigh, rand_lines]))
        }
        
      } else {
        PPI_sel<-PPI_mat[g, ]
        first_neigh<-names(which(PPI_sel==1))
        PPI_sel<-PPI_mat[first_neigh, ]
        first_neigh<-names(which(PPI_sel==1))
        first_neigh<-setdiff(intersect(first_neigh, rownames(scaled_depFC)),g)
        rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
        rand_lines<-setdiff(rand_lines, line)
        if(length(first_neigh)>0){
          actual2<-c(actual2,unlist(scaled_depFC[first_neigh, line]))
          others2<-c(others2,unlist(scaled_depFC[first_neigh, rand_lines]))
          
        }
      }
      
    }
  }
}

df<-data.frame(logFC=c(actual2, others2), type=c(rep("DAM-bearing lines", length(actual2)), rep("other lines", length(others2))))
ggplot(df, aes(x=logFC, colour=type, group=type))+geom_density()
ks.test(actual2, others2, alternative = "greater")
ks.test(actual2, others2, alternative = "less")
median(actual2)-median(others2)

#######
##final
#########

actual<-list()
others<-list()
for(neigh_order in 1:5){
  actual[[neigh_order]]<-c(NA)
  others[[neigh_order]]<-c(NA)
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)
  gene<-results[[ctiss]][ind,"GENE"]
  cl<-results[[ctiss]][ind,"ps_cl"]
  
  for(i in 1:length(ind)){
    g<-gene[i]
    line<-cl[i]
    line<-unlist(strsplit(line, ", "))
    
    if(g %in% rownames(PPI_mat)){
      PPI_sel<-PPI_mat[g, ]
      first_neigh<-names(which(PPI_sel==1))
      if(neigh_order>1){
      for(iter in 2:neigh_order){
        PPI_sel<-PPI_mat[first_neigh, ]
        if(is.null(nrow(PPI_sel))){
          first_neigh<-names(which(PPI_sel==1))
        } else {
        first_neigh<-colnames(PPI_sel)[unique(which(PPI_sel==1, arr.ind=T)[,2])]
        }
      }
      }
      first_neigh<-setdiff(intersect(first_neigh, rownames(scaled_depFC)),g)
      
        rand_lines<-intersect(colnames(scaled_depFC), CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
        rand_lines<-setdiff(rand_lines, line)
        
        if(length(first_neigh)>0){
          actual[[neigh_order]]<-c(actual[[neigh_order]],unlist(scaled_depFC[first_neigh, line]))
          others[[neigh_order]]<-c(others[[neigh_order]],unlist(scaled_depFC[first_neigh, rand_lines]))
          
        }
      }
      
    }
  }
}

md<-c()
for(neigh_order in 1:5){
  md<-c(md, median(actual[[neigh_order]], na.rm=T)-median(others[[neigh_order]], na.rm=T))
}
pdf("Neighbor_order.pdf",7,7)
plot(1:5, md, pch=19, xlab="neighbour order", ylab="median difference")
dev.off()

df<-data.frame(logFC=c(actual[[1]], others[[1]]), type=c(rep("DAM-bearing lines", length(actual[[1]])), rep("other lines", length(others[[1]]))))
pdf("Distribution_logFC_neighbours.pdf",8,8)
ggplot(df, aes(x=logFC, colour=type, group=type))+geom_density()+theme_classic()
dev.off()
ks.test(actual2, others2, alternative = "greater")
ks.test(actual2, others2, alternative = "less")
median(actual2)-median(others2)

