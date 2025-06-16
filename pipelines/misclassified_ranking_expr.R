.libPaths(c("/home/aurora.savino/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))
######################
### Cell lines with mutations in driver genes of other tissues
### Are they mis-classified cell lines?
######################
##import expression data
library(openxlsx)
setwd("~/VUS/VUS/results/20220208/")
load("../../data/rnaseq_fpkm_20210329.RData")
expression<-log2(expression+1)

library(tidyverse)
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
inTOgen_drivers<-read.table('../../data/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv',sep='\t',stringsAsFactors = FALSE,header=TRUE)
driver_genes<-unique(inTOgen_drivers$SYMBOL)


##############################
##taking tumor burden into account
##############################
gene_annot <- read_csv("../../data/gene_identifiers_20191101.csv")
load('../../data/R/Sanger_Broad_higQ_scaled_depFC.RData')

CMP_annot <- read_csv("../../data/model_list_20210611.csv") # from https://cog.sanger.ac.uk/cmp/download/model_list_20210611.csv
CMP_annot<-CMP_annot[which(is.element(CMP_annot$model_name,colnames(scaled_depFC))),]
     
## latest sanger/broad unreleased yet variants hg38
cl_variants <- read_csv('../../data/mutations_all_latest.csv')
cl_variants <- cbind(cl_variants,gene_annot$hgnc_symbol[match(cl_variants$gene_id,gene_annot$gene_id)])
colnames(cl_variants)[ncol(cl_variants)]<-'gene_symbol_2019'

cl_variants <- cbind(cl_variants,CMP_annot$model_name[match(cl_variants$model_id,CMP_annot$model_id)])
colnames(cl_variants)[ncol(cl_variants)] <- 'model_name'


cl_variants <- cbind(cl_variants,CMP_annot$cancer_type[match(cl_variants$model_id,CMP_annot$model_id)])
colnames(cl_variants)[ncol(cl_variants)] <- 'cancer_type'

colnames(cl_variants)[12] <- 'protein_mutation'

cancer_match<-read.xlsx("cancer_match_ext.xlsx")

###testare se le linee cellulari dello stesso tessuto sono più simili tra loro o meno
all_lines<-unique(CMP_annot$model_name[CMP_annot$cancer_type %in% tissues])
all_lines_sel<-intersect(all_lines, colnames(expression))
extracor<-cor(expression[,all_lines_sel])

##�due tipi di randomizzazione
#1) int_lines vs random ext_lines scelte tra quelle dei tumori senza il driver
#2) int_lines vs random ext_lines

set.seed(75973759)
misclass_result_all <-as.data.frame(matrix(nrow=1,ncol=33*2+6))
colnames(misclass_result_all)<-c("line", "gene", "inexpr", "tumor", "other_type", "average_int_corr", tissues, paste(tissues, "pval"))

for(ctiss in tissues){
  print(ctiss)
  known<-inTOgen_drivers$SYMBOL[inTOgen_drivers$CANCER_TYPE %in% (cancer_match[which(cancer_match[,1]==ctiss),2])]
  total_found<-intersect(results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2], driver_genes)
  novel_all<-setdiff(total_found, known)
  
  for(novel in novel_all){
  ##cell lines with mutation in the novel genes
  all_lines_ctiss<-unique(CMP_annot$model_name[CMP_annot$cancer_type==ctiss])
  lines_ctiss_sel<-intersect(CMP_annot$model_name, unique(cl_variants$model_name[cl_variants$gene_symbol==novel&cl_variants$cancer_type==ctiss]))
  
  all_lines_ctiss<-intersect(all_lines_ctiss, colnames(expression))
  lines_ctiss_sel<-intersect(lines_ctiss_sel, colnames(expression))
  
  type<-na.omit(unique(cancer_match[which(cancer_match[,2] %in% inTOgen_drivers$CANCER_TYPE[inTOgen_drivers$SYMBOL==novel]), 1]))
  
  if(!is.null(type) & length(type)>0 & length(lines_ctiss_sel)>0 & length(all_lines_ctiss)>1){
  #pi� di un tipo
  for(othertype in type){
    misclass_result <-as.data.frame(matrix(nrow=length(lines_ctiss_sel),ncol=33+6))
    colnames(misclass_result)<-c("line", "gene", "inexpr", "tumor", "other_type", "average_int_corr", tissues)
    rownames(misclass_result)<-lines_ctiss_sel
    misclass_result$line<-lines_ctiss_sel 
    misclass_result$gene<-rep(novel, length(lines_ctiss_sel))
    
    misclass_result$inexpr<-ifelse(misclass_result$line %in% lines_ctiss_sel, "Yes", "No")
    misclass_result$other_type<-othertype
    misclass_result$tumor<-ctiss
    
    cor_int<-cor(expression[,all_lines_ctiss])
    
    if(length(lines_ctiss_sel)>1){
      for(line in lines_ctiss_sel){
        cor_int_line<-cor_int[line,-which(colnames(cor_int)==line)]
        
        misclass_result$average_int_corr[misclass_result$line==line]<-mean(cor_int_line)
        corrs<-c()
        for(other_tiss in tissues){
          all_lines_othertype<-unique(CMP_annot$model_name[CMP_annot$cancer_type==other_tiss])
          all_lines_othertype<-intersect(all_lines_othertype, colnames(expression))
          cor_ext_line<-extracor[line,all_lines_othertype]
          corrs<-c(corrs, mean(cor_ext_line))
        }
        misclass_result[line,tissues]<-corrs
        }
      
      
    } else if(length(lines_ctiss_sel)==1){
      cor_int_line<-cor_int[lines_ctiss_sel,-which(colnames(cor_int)==lines_ctiss_sel)]
      misclass_result$average_int_corr<-mean(cor_int_line)
      
      corrs<-c()
      for(other_tiss in tissues){
        all_lines_othertype<-unique(CMP_annot$model_name[CMP_annot$cancer_type==other_tiss])
        all_lines_othertype<-intersect(all_lines_othertype, colnames(expression))
        cor_ext_line<-extracor[lines_ctiss_sel,all_lines_othertype]
        corrs<-c(corrs, mean(cor_ext_line))
      }
      misclass_result[lines_ctiss_sel,tissues]<-corrs
    }
    
    misclass_result<-misclass_result[which(misclass_result$inexpr=="Yes"),]
    
    ####randomization
    
    misclass_rand<-list()
    
    for(iter in 1:100){
      ind<-0
      #reshuffle labels of other tissues
      toreshuff<-setdiff(all_lines_sel, all_lines_ctiss)
      reshuffled<-toreshuff[sample(1:length(toreshuff), length(toreshuff), replace = F)]
      extracor_reshuff<-extracor
      colnames(extracor_reshuff)[match(toreshuff, colnames(extracor_reshuff))]<-reshuffled
      
      if(length(lines_ctiss_sel)>1){
        for(line in lines_ctiss_sel){
          
          ind<-ind+1
          
          if(iter==1){
            misclass_rand[[ind]]<-matrix(nrow=100, ncol=33)
            colnames(misclass_rand[[ind]])<-tissues
          }
          
          corrs<-c()
          for(other_tiss in tissues){
            all_lines_othertype<-unique(CMP_annot$model_name[CMP_annot$cancer_type==other_tiss])
            all_lines_othertype<-intersect(all_lines_othertype, colnames(expression))
            cor_ext_line<-extracor_reshuff[line,all_lines_othertype]
            corrs<-c(corrs, mean(cor_ext_line))
          }
          misclass_rand[[ind]][iter,tissues]<-corrs
        }
        
        
      } else if(length(lines_ctiss_sel)==1){
        
        ind<-ind+1
        
        if(iter==1){
          misclass_rand[[ind]]<-matrix(nrow=100, ncol=33)
          colnames(misclass_rand[[ind]])<-tissues
        }
      
        corrs<-c()
        for(other_tiss in tissues){
          all_lines_othertype<-unique(CMP_annot$model_name[CMP_annot$cancer_type==other_tiss])
          all_lines_othertype<-intersect(all_lines_othertype, colnames(expression))
          cor_ext_line<-extracor_reshuff[lines_ctiss_sel,all_lines_othertype]
          corrs<-c(corrs, mean(cor_ext_line))
        }
        misclass_rand[[ind]][iter,tissues]<-corrs
      }
    }
    if(nrow(misclass_result)!=length(misclass_rand)){
      print("Not matching")
      print(ctiss)
      print(novel)
    }
    
    
    misclass_result_p<-matrix(nrow=nrow(misclass_result), ncol=length(tissues))
    colnames(misclass_result_p)<-tissues
    for(ind in 1:nrow(misclass_result)){
    empirical_p<-colSums(t(t(misclass_rand[[ind]][,tissues])>=unlist(misclass_result[ind,tissues])))/100
    misclass_result_p[ind,tissues]<-empirical_p
    }
    colnames(misclass_result_p)<-paste(colnames(misclass_result_p), "pval")
    misclass_result<-cbind.data.frame(misclass_result, misclass_result_p)
    
  }
    
    misclass_result_all<-rbind(misclass_result_all, misclass_result)
    
  }
  }
}
   
save(misclass_result_all, file="misclass_expr_rank2.RData")
