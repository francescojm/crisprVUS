set.seed(123)
library(CELLector)
library(tidyverse)

####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250221/'

######################
### Cell lines with mutations in driver genes of other tissues
### Are they mis-classified cell lines?
######################
##import expression data
load(paste(pathdata, "/Robj/basal_exp.RData", sep=""))
load(paste(pathdata, "/Robj/proteome_mat.RData", sep=""))
load(paste(pathdata, "/Robj/drugs_lnIC50_gdscAll.RData", sep=""))

###loading input data
gene_annot <- read_csv(paste(pathdata, "/raw/gene_identifiers_20241212.csv", sep=""))
### gene_identifiers_20191101 downloaded from https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20241212.csv on 20250221

CMP_annot <- read_csv(paste(pathdata,"/raw/model_list_20241120.csv", sep="")) 
### model_list_20240110.csv downloaded from https://cog.sanger.ac.uk/cmp/download/model_list_20241120.csv on 20250129

scaled_depFC<-read.csv(paste(pathdata,'/raw/CRISPRGeneEffect.csv', sep=""), row.names = 1)
colnames(scaled_depFC)<-gsub("\\..*","",colnames(scaled_depFC))
###scaled essentiality matrices downloaded from https://depmap.org/portal/data_page/?tab=allData on 20250129 (24Q4)

toremove<-which(is.na(CMP_annot$model_id[match(rownames(scaled_depFC),CMP_annot$BROAD_ID)]))
scaled_depFC<-scaled_depFC[-toremove,]
rownames(scaled_depFC)<-CMP_annot$model_id[match(rownames(scaled_depFC),CMP_annot$BROAD_ID)]

scaled_depFC<-t(scaled_depFC)

#remove genes with missing values
scaled_depFC<-scaled_depFC[-which(rowSums(is.na(scaled_depFC))>0),]
#create a binarized matrix (essential vs non-essential)
bdep<-scaled_depFC
bdep[scaled_depFC<=(-0.5)]<-1
bdep[scaled_depFC>(-0.5)]<-0

load(paste(pathdata,'/Robj/ADaM.RData', sep=""))
load(paste(pathdata,'/Robj/FiPer_outputs.RData', sep=""))
### Rbjects precomputed as in Vinceti et al, BMC Genomics, 2021

cl_variants <- read_csv(paste(pathdata,'/raw/mutations_all_20241212.csv', sep=""))
### mutations_all_20230202 downloaded from  https://cog.sanger.ac.uk/cmp/download/mutations_all_20241212.zip on 20250221

cl_variants <- cbind(cl_variants,gene_annot$hgnc_symbol[match(cl_variants$gene_id,gene_annot$gene_id)])
colnames(cl_variants)[ncol(cl_variants)]<-'gene_symbol_2023'


#select only cell lines with data in both CRISPR screens and CMP annotation file
CMP_annot<-CMP_annot[which(is.element(CMP_annot$model_id,colnames(scaled_depFC))),]

#select tissues
tissues<-CMP_annot$cancer_type
st<-summary(as.factor(tissues))

tissues<-sort(setdiff(tissues,names(which(st<5))))
tissues<-setdiff(tissues,c('Other Solid Carcinomas','Other Solid Cancers','Other Sarcomas', "Other Blood Cancers", "Non-Cancerous"))

CMP_annot<-CMP_annot[which(is.element(CMP_annot$cancer_type,tissues)),]

tissues<-sort(tissues)


##driver genes
inTOgen_drivers<-read.table(paste(pathdata,"/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", sep=""), sep='\t',stringsAsFactors = FALSE,header=TRUE)
driver_genes<-unique(inTOgen_drivers$SYMBOL)
  
cancer_match<-read.csv(paste(pathdata,'/raw/intOGen ctype mapping_AS.csv',sep=''),header = TRUE,row.names = 1, sep=";")

setwd(resultPath)

tissues_res<-gsub("_results.RData", "", list.files(pattern="results.RData"))

results<-list()
ind_iter<-0
for(ctiss in tissues_res){
  ind_iter<-ind_iter+1
  load(paste(ctiss, "_results_ext.RData", sep=""))
  results[[ind_iter]] <- RESTOT
}
names(results)<-tissues_res

cancer_match_long_CMP<-c()
cancer_match_long_into<-c()
for(i in 1:nrow(cancer_match)){
  cancer_match_long_into<-c(cancer_match_long_into, unlist(strsplit(cancer_match[i,1], " \\| ")))
  cancer_match_long_CMP<-c(cancer_match_long_CMP, rep(rownames(cancer_match)[i],length(unlist(strsplit(cancer_match[i,1], " \\| ")))))
}
cancer_match_long<-cbind(cancer_match_long_CMP, cancer_match_long_into)

##################################
### count the number of gene-cell line to be tested
###################################

comb_tests<-c()
  for(ctiss in tissues_res){
    print(ctiss)
    #all genes known to be driver in the specific cancer type
    known<-unique(inTOgen_drivers$SYMBOL[ inTOgen_drivers$CANCER_TYPE %in% setdiff(unlist(strsplit( cancer_match[ctiss,1], " | ")), "|")])
    #all genes known to be driver in any cancer type that are found as DAM-bearing
    total_found<-intersect(results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand < 0.2], driver_genes)
    novel_all<-setdiff(total_found, known)
    
    for(novel in novel_all){
      print(novel)
      type<-na.omit(unique(cancer_match_long[which(cancer_match_long[,2] %in% inTOgen_drivers$CANCER_TYPE[inTOgen_drivers$SYMBOL==novel]), 1]))
      comb<-paste(ctiss, unlist(strsplit(results[[ctiss]]$ps_cl[results[[ctiss]]$GENE==novel], ", ")))
      
    ##cell lines with mutation in the novel genes found through our analysis
      comb_tests<-c(comb_tests, paste(rep(comb, each = length(type)), type, sep = "."))
      
    }
  }
print(paste("combinations of gene-cell line- alternate cancer type to be tested for misclassification:", length(comb_tests)))
print(paste("number of cell lines-alternate cancer type to be tested for misclassification:", length(unique(comb_tests))))

#################################################
###test whether cell lines of the same tissue are more similar to the selected cl than cell lines of other tissues 
################################################

cl_similarity<-function(omic, omic_name){
  inn<-0
  all_lines<-unique(CMP_annot$model_id[CMP_annot$cancer_type %in% tissues])
all_lines_sel<-intersect(all_lines, colnames(omic))
extracor<-cor(omic[,all_lines_sel], use = "pairwise.complete.obs")

set.seed(75973759)
misclass_result_all <-as.data.frame(matrix(nrow=1,ncol=length(tissues)*2+6))
colnames(misclass_result_all)<-c("line", "gene", "inexpr", "tumor", "other_type", "average_int_corr", tissues, paste(tissues, "pval"))

for(ctiss in tissues_res){
  print(ctiss)
  #all genes known to be driver in the specific cancer type
  known<-unique(inTOgen_drivers$SYMBOL[ inTOgen_drivers$CANCER_TYPE %in% setdiff(unlist(strsplit( cancer_match[ctiss,1], " | ")), "|")])
  #all genes known to be driver in any cancer type that are found as DAM-bearing
  total_found<-intersect(results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand < 0.2], driver_genes)
  novel_all<-setdiff(total_found, known)
  
  for(novel in novel_all){
  
  all_lines_ctiss<-unique(CMP_annot$model_id[CMP_annot$cancer_type==ctiss])
  ##cell lines with mutation in the novel genes found through our analysis
  lines_ctiss_sel<-unlist(strsplit(results[[ctiss]]$ps_cl[results[[ctiss]]$GENE==novel], ", "))
  
  all_lines_ctiss<-intersect(all_lines_ctiss, colnames(omic))
  lines_ctiss_sel<-intersect(lines_ctiss_sel, colnames(omic))
  
  type<-na.omit(unique(cancer_match_long[which(cancer_match_long[,2] %in% inTOgen_drivers$CANCER_TYPE[inTOgen_drivers$SYMBOL==novel]), 1]))
  
  if(!is.null(type) & length(type)>0 & length(lines_ctiss_sel)>0 & length(all_lines_ctiss)>1){
  
  for(othertype in type){
    inn<-inn+1
    print(inn)
    misclass_result <-as.data.frame(matrix(nrow=length(lines_ctiss_sel),ncol=length(tissues)+6))
    colnames(misclass_result)<-c("line", "gene", "inexpr", "tumor", "other_type", "average_int_corr", tissues)
    rownames(misclass_result)<-lines_ctiss_sel
    misclass_result$line<-lines_ctiss_sel 
    misclass_result$gene<-rep(novel, length(lines_ctiss_sel))
    
    misclass_result$inexpr<-ifelse(misclass_result$line %in% lines_ctiss_sel, "Yes", "No")
    misclass_result$other_type<-othertype
    misclass_result$tumor<-ctiss
    
    cor_int<-cor(omic[,all_lines_ctiss], use = "pairwise.complete.obs")
    
    if(length(lines_ctiss_sel)>1){
      for(line in lines_ctiss_sel){
        cor_int_line<-cor_int[line,-which(colnames(cor_int)==line)]
        
        misclass_result$average_int_corr[misclass_result$line==line]<-mean(cor_int_line)
        corrs<-c()
        for(other_tiss in tissues){
          all_lines_othertype<-unique(CMP_annot$model_id[CMP_annot$cancer_type==other_tiss])
          all_lines_othertype<-intersect(all_lines_othertype, colnames(omic))
          cor_ext_line<-extracor[line,all_lines_othertype]
          corrs<-c(corrs, mean(cor_ext_line, na.rm=T))
        }
        misclass_result[line,tissues]<-corrs
        }
      
      
    } else if(length(lines_ctiss_sel)==1){
      cor_int_line<-cor_int[lines_ctiss_sel,-which(colnames(cor_int)==lines_ctiss_sel)]
      misclass_result$average_int_corr<-mean(cor_int_line)
      
      corrs<-c()
      for(other_tiss in tissues){
        all_lines_othertype<-unique(CMP_annot$model_id[CMP_annot$cancer_type==other_tiss])
        all_lines_othertype<-intersect(all_lines_othertype, colnames(omic))
        cor_ext_line<-extracor[lines_ctiss_sel,all_lines_othertype]
        corrs<-c(corrs, mean(cor_ext_line, na.rm=T))
      }
      misclass_result[lines_ctiss_sel,tissues]<-corrs
    }
    
    misclass_result<-misclass_result[which(misclass_result$inexpr=="Yes"),]
    
    ####randomization
    
    misclass_rand<-list()
    
    for(iter in 1:1000){
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
            misclass_rand[[ind]]<-matrix(nrow=1000, ncol=length(tissues))
            colnames(misclass_rand[[ind]])<-tissues
          }
          
          corrs<-c()
          for(other_tiss in tissues){
            all_lines_othertype<-unique(CMP_annot$model_id[CMP_annot$cancer_type==other_tiss])
            all_lines_othertype<-intersect(all_lines_othertype, colnames(omic))
            cor_ext_line<-extracor_reshuff[line,all_lines_othertype]
            corrs<-c(corrs, mean(cor_ext_line, na.rm=T))
          }
          misclass_rand[[ind]][iter,tissues]<-corrs
        }
        
        
      } else if(length(lines_ctiss_sel)==1){
        
        ind<-ind+1
        
        if(iter==1){
          misclass_rand[[ind]]<-matrix(nrow=1000, ncol=length(tissues))
          colnames(misclass_rand[[ind]])<-tissues
        }
      
        corrs<-c()
        for(other_tiss in tissues){
          all_lines_othertype<-unique(CMP_annot$model_id[CMP_annot$cancer_type==other_tiss])
          all_lines_othertype<-intersect(all_lines_othertype, colnames(omic))
          cor_ext_line<-extracor_reshuff[lines_ctiss_sel,all_lines_othertype]
          corrs<-c(corrs, mean(cor_ext_line, na.rm=T))
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
    empirical_p<-colSums(t(t(misclass_rand[[ind]][,tissues])>=unlist(misclass_result[ind,tissues])))/1000
    misclass_result_p[ind,tissues]<-empirical_p
    }
    colnames(misclass_result_p)<-paste(colnames(misclass_result_p), "pval")
    misclass_result<-cbind.data.frame(misclass_result, misclass_result_p)
    
    misclass_result_all<-rbind(misclass_result_all, misclass_result)
  }
    
    
    
  }
  }
}
   
save(misclass_result_all, file=paste("misclass_", omic_name, "_rank.RData", sep=""))
}

#######################
## run with different omics
######################

cl_similarity(omic=basal_exp, omic_name = "expr")
cl_similarity(omic=proteome_mat, omic_name = "prot")
cl_similarity(omic=drugs_lnIC50, omic_name = "drugs")
cl_similarity(omic=scaled_depFC, omic_name = "ess")
