########################################
##Generate tables for VUS database
#########################################
library(tidyverse)

path_data<-"/data"
path_results<-"/results/20250221"
home<-"E:/VUS_2024build"

setwd(paste(home, "/", path_results, sep=""))
dir.create("DB_tables")

##load essentiality data
CMP_annot <- read_csv(paste(home, path_data,"/raw/model_list_20241120.csv", sep="")) 
### model_list_20240110.csv downloaded from https://cog.sanger.ac.uk/cmp/download/model_list_20241120.csv on 20250129

scaled_depFC<-read.csv(paste(home, path_data,'/raw/CRISPRGeneEffect.csv', sep=""), row.names = 1)
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

##############################
##cell lines annotations
#############################
#select only cell lines with data in both CRISPR screens and CMP annotation file
CMP_annot<-CMP_annot[which(is.element(CMP_annot$model_id,colnames(scaled_depFC))),]

#select tissues
tissues<-CMP_annot$cancer_type
st<-summary(as.factor(tissues))

tissues<-sort(setdiff(tissues,names(which(st<5))))
tissues<-setdiff(tissues,c('Other Solid Carcinomas','Other Solid Cancers','Other Sarcomas', "Other Blood Cancers", "Non-Cancerous"))

CMP_annot<-CMP_annot[which(is.element(CMP_annot$cancer_type,tissues)),]

tissues<-sort(tissues)

tissue_names<-data.frame(tissue_name=tissues)
write.csv(tissue_names, "DB_tables/tissues.csv", quote=F, row.names = F)


CellLine<-CMP_annot[, c("model_id", "model_name", "cancer_type")]
colnames(CellLine)<-c("cell_line_id", "cell_line_name", "tissue_name")
#check that the id is unique
sum(duplicated(CellLine$CellLineID))==0

write.csv(CellLine, "DB_tables/CellLine.csv", quote=F, row.names = F)


#######################
### Variants
####################
###loading input data
gene_annot <- read_csv(paste(home, path_data, "/raw/gene_identifiers_20241212.csv", sep=""))
### gene_identifiers_20191101 downloaded from https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20241212.csv on 20250221

cl_variants <- read_csv(paste(home, path_data,'/raw/mutations_all_20241212.csv', sep=""))
### mutations_all_20230202 downloaded from  https://cog.sanger.ac.uk/cmp/download/mutations_all_20241212.zip on 20250221

cl_variants <- cbind(cl_variants,gene_annot$hgnc_symbol[match(cl_variants$gene_id,gene_annot$gene_id)])
colnames(cl_variants)[ncol(cl_variants)]<-'gene_symbol_2023'

##cancer type
cl_variants <- cbind(cl_variants,CMP_annot$cancer_type[match(cl_variants$model_id,CMP_annot$model_id)])
colnames(cl_variants)[ncol(cl_variants)] <- 'cancer_type'

##used for gene filtering
load(paste(home, path_data,'/RObj/ADaM.RData', sep=""))
load(paste(home, path_data,'/RObj/FiPer_outputs.RData', sep=""))
### Rbjects precomputed as in Vinceti et al, BMC Genomics, 2021

##########
### Load DAMs
##############
tissues<-gsub("_results_ext.RData", "", list.files(pattern="results_ext.RData"))

results<-list()
ind_iter<-0
for(ctiss in tissues){
  ind_iter<-ind_iter+1
  load(paste(ctiss, "_results_ext.RData", sep=""))
  results[[ind_iter]] <- RESTOT
}
names(results)<-tissues

####define all selected vars

vars_tot<-c()
genes_tot<-c()
tiss_tot<-c()
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)
  vars<-strsplit(results[[ctiss]]$var[ind], " \\| ")
  genes<-rep(results[[ctiss]]$GENE[ind], lapply(vars, length))
  vars_tot<-c(vars_tot, paste(genes, unlist(vars), sep="-"))
  genes_tot<-c(genes_tot, genes)
  tiss_tot<-c(tiss_tot, rep(ctiss, length(genes)))
}

DAMs<-cbind(gsub("p\\.", "", vars_tot), genes_tot, tiss_tot)
colnames(DAMs)<-c("variant_id", "gene_id", "tissue")
write.csv(DAMs, file="DB_tables/DAM.csv", quote=F, row.names = F)

allvars<-unique(vars_tot)
allgenes<-unique(genes_tot)

scaled_depFC_sel<-scaled_depFC[allgenes,unique(CMP_annot$model_id)]

library(reshape2)
scaled_depFC_sel_melt<-melt(data = scaled_depFC_sel, varnames = c("gene", "line"))
Essentiality<-scaled_depFC_sel_melt
colnames(Essentiality)<-c("gene_id", "cell_line_id", "essentiality_value")

#check that the id is unique
sum(duplicated(paste(Essentiality$gene_id, Essentiality$cell_line_id)))==0
write.csv(Essentiality, file="DB_tables/Essentiality.csv", quote=F, row.names = F)

#####################

CellLineHasVariant<-cl_variants[(cl_variants$model_id %in% CellLine$cell_line_id & paste(cl_variants$gene_symbol_2023, cl_variants$protein_mutation, sep="-") %in% allvars),]
CellLineHasVariant$VariantID<-paste(CellLineHasVariant$gene_symbol_2023, CellLineHasVariant$protein_mutation, sep="-")
CellLineHasVariant<-CellLineHasVariant[,c("model_id", "VariantID")]
colnames(CellLineHasVariant)<-c("CellLineID", "VariantID")

#check that the id is unique
CellLineHasVariant<-unique(CellLineHasVariant)
sum(duplicated(paste(CellLineHasVariant$CellLineID, CellLineHasVariant$VariantID)))==0

#togliere le lettere maiuscole subito dopo "del" in variants
CellLineHasVariant$VariantID[grep("del[[:upper:]]+$",CellLineHasVariant$VariantID)]<-gsub("del[[:upper:]]+$","del",CellLineHasVariant$VariantID[grep("del[[:upper:]]+$",CellLineHasVariant$VariantID)])

CellLineHasVariant$VariantID<-gsub("p.", "", CellLineHasVariant$VariantID)

CellLineHasVariant$DAM<-0
CellLineHasVariant$DAM[which(paste(CellLineHasVariant$VariantID, CMP_annot$cancer_type[match(CellLineHasVariant$CellLineID, CMP_annot$model_id)]) %in% gsub("p\\.", "", paste(vars_tot, tiss_tot)))]<-1
CellLineHasVariant$tissue_name<-CMP_annot$cancer_type[match(CellLineHasVariant$CellLineID, CMP_annot$model_id)]

colnames(CellLineHasVariant)<-c("cell_line_id", "variant_id", "dam", "tissue_name")
CellLineHasVariant<-CellLineHasVariant[,c("cell_line_id", "variant_id",  "tissue_name", "dam")]
write.csv(CellLineHasVariant, file="DB_tables/CellLineHasVariant.csv", quote=F, row.names = F)

#####################################################################
################################### BY VAR
#########################################################################
inTOgen_drivers<-read.table(paste(home, '/data/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
### inTOgene drivers downloaded from https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip on 20241002
driver_genes<-unique(inTOgen_drivers$SYMBOL)

genes<-cbind(gene_annot$hgnc_symbol, ifelse(gene_annot$hgnc_symbol %in% driver_genes, 1, 0))
colnames(genes)<-c("gene", "is_driver")
genes<-genes[-which(duplicated(genes[,"gene"])),]
write.csv(genes, file="DB_tables/genes.csv", quote=F, row.names = F)

#################
#### each cancer type separately
###################

mapping<-read.csv(paste(home, path_data,'/raw/intOGen ctype mapping_AS.csv',sep=''),header = TRUE,row.names = 1, sep=";")

cancer_match_long_CMP<-c()
cancer_match_long_into<-c()
for(i in 1:nrow(mapping)){
  cancer_match_long_into<-c(cancer_match_long_into, unlist(strsplit(mapping[i,1], " \\| ")))
  cancer_match_long_CMP<-c(cancer_match_long_CMP, rep(rownames(mapping)[i],length(unlist(strsplit(mapping[i,1], " \\| ")))))
}
mapping_intogen<-cbind(cancer_match_long_CMP, cancer_match_long_into)

cancer_match_long_CMP<-c()
cancer_match_long_into<-c()
for(i in 1:nrow(mapping)){
  cancer_match_long_into<-c(cancer_match_long_into, unlist(strsplit(mapping[i,2], " \\| ")))
  cancer_match_long_CMP<-c(cancer_match_long_CMP, rep(rownames(mapping)[i],length(unlist(strsplit(mapping[i,2], " \\| ")))))
}
mapping_cosmic<-cbind(cancer_match_long_CMP, cancer_match_long_into)

load("summary_byvar.RData")
summary_vars$ID<-paste(summary_vars$Gene, summary_vars$Var, sep="-")
summary_vars$Num_tiss<-summary_vars$Perc_tiss*33
summary_vars$Num_Intogen[is.na(summary_vars$Num_Intogen)]<-0

summary_vars$Num_tiss_Intogen<-NA
summary_vars$Num_tiss_COSMIC<-NA
summary_vars$Match_tiss_Intogen<-NA
summary_vars$Match_tiss_COSMIC<-NA
summary_vars$Matching_inboth<-NA
summary_vars$Matching_union<-NA
summary_vars$Match_tiss_Intogen_name<-NA
summary_vars$Match_tiss_COSMIC_name<-NA

load("all_tiss_comsic_byvar.RData")
load("all_tiss_intogen_byvar.RData")

for(i in 1:nrow(summary_vars)){
  type_tiss<-unique(unlist(strsplit(summary_vars[i, "Tiss"], " \\| ")))
  type_intogen<-na.omit(unique(unlist(strsplit(summary_vars[i, "Type_Intogen"], " \\| "))))
  type_cosmic<-na.omit(unique(unlist(strsplit(summary_vars[i, "Type_COSMIC"], " \\| "))))
  
  summary_vars$Num_tiss_Intogen[i]<-length(type_intogen)
  summary_vars$Num_tiss_COSMIC[i]<-length(type_cosmic)
  
  matching_intogen<-unique(mapping_intogen[which(mapping_intogen[,2] %in% type_intogen),1])
  
  matching_cosmic<-unique(mapping_cosmic[which(mapping_cosmic[,2] %in% type_cosmic),1])
  
  summary_vars$Match_tiss_Intogen[i]<-length(intersect(type_tiss, matching_intogen))
  summary_vars$Match_tiss_COSMIC[i]<-length(intersect(type_tiss, matching_cosmic))
  summary_vars$Matching_inboth[i]<-length(Reduce(intersect, list(matching_intogen, type_tiss, matching_cosmic)))
  summary_vars$Matching_union[i]<-length(union(intersect(matching_intogen, type_tiss), intersect(matching_cosmic, type_tiss)))
  summary_vars$Match_tiss_Intogen_name[i]<-paste(intersect(type_tiss, matching_intogen), collapse="|")
  summary_vars$Match_tiss_COSMIC_name[i]<-paste(intersect(type_tiss, matching_cosmic), collapse="|")
}

summary_vars_bytiss<-rbind(summary_vars,summary_vars)
summary_vars_bytiss$tissue<-"ALL"
summary_vars_bytiss$num_patients<-c(summary_vars[,c("Num_COSMIC")], summary_vars[,c( "Num_Intogen")])
summary_vars_bytiss$database<-c(rep("COSMIC", nrow(summary_vars)), rep("Intogen", nrow(summary_vars)))

for(i in 1:nrow(summary_vars)){
  type_tiss<-unique(unlist(strsplit(summary_vars[i, "Tiss"], " \\| ")))
  type_intogen<-na.omit(unique(unlist(strsplit(summary_vars[i, "Type_Intogen"], " \\| "))))
  type_cosmic<-na.omit(unique(unlist(strsplit(summary_vars[i, "Type_COSMIC"], " \\| "))))
  
  matching_intogen<-unique(mapping_intogen[which(mapping_intogen[,2] %in% type_intogen),1])
  matching_cosmic<-unique(mapping_cosmic[which(mapping_cosmic[,2] %in% type_cosmic),1])
  
  ###matched tissues cosmic
  if(length(intersect(type_tiss, matching_cosmic))>0){
    num_cosmic_bytiss<-c()
    for(tiss in intersect(type_tiss, matching_cosmic)){
      tiss_cosmic<-mapping_cosmic[mapping_cosmic[,1]==tiss,2]
      num_cosmic_bytiss<-c(num_cosmic_bytiss, sum(unlist(all_tiss_cosmic[[summary_vars$ID[i]]][tiss_cosmic]), na.rm=T))
    }
    num_cosmic_bytiss<-c(num_cosmic_bytiss, sum(num_cosmic_bytiss))
    names(num_cosmic_bytiss)<-c(intersect(type_tiss, matching_cosmic), "ALL_MATCHING")
    toadd<-data.frame(summary_vars[i,], num_patients=num_cosmic_bytiss)
    toadd$tissue<-names(num_cosmic_bytiss)
    toadd$database<-"COSMIC"
    summary_vars_bytiss<-rbind(summary_vars_bytiss, toadd)
  } else {
    toadd<-data.frame(summary_vars[i,], num_patients=0)
    toadd$tissue<-"ALL_MATCHING"
    toadd$database<-"COSMIC"
    summary_vars_bytiss<-rbind(summary_vars_bytiss, toadd)
  }
  
  ###matched tissues intogen
  if(length(intersect(type_tiss, matching_intogen))>0){
    num_intogen_bytiss<-c()
    for(tiss in intersect(type_tiss, matching_intogen)){
      tiss_intogen<-mapping_intogen[mapping_intogen[,1]==tiss,2]
      num_intogen_bytiss<-c(num_intogen_bytiss, sum(unlist(all_tiss_intogen[[summary_vars$ID[i]]][tiss_intogen]), na.rm=T))
    }
    num_intogen_bytiss<-c(num_intogen_bytiss, sum(num_intogen_bytiss))
    names(num_intogen_bytiss)<-c(intersect(type_tiss, matching_intogen), "ALL_MATCHING")
    toadd<-data.frame(summary_vars[i,], num_patients=num_intogen_bytiss)
    toadd$tissue<-names(num_intogen_bytiss)
    toadd$database<-"Intogen"
    summary_vars_bytiss<-rbind(summary_vars_bytiss, toadd)
  } else {
    toadd<-data.frame(summary_vars[i,], num_patients=0)
    toadd$tissue<-"ALL_MATCHING"
    toadd$database<-"Intogen"
    summary_vars_bytiss<-rbind(summary_vars_bytiss, toadd)
  }
  
}

###########################
#### add VEP
##############################

vep <- readRDS(paste(home, '/data/raw/vep_annotations.rds', sep=""))
vars_vep<-vep$GRCh38[,c("gene_name", "protein_mutation", "sift_prediction", "polyphen_prediction")]
vars_vep$var<-paste(vars_vep$gene_name, gsub("p.", "", vars_vep$protein_mutation), sep="-")
vars_vep$polyphen_prediction[is.na(vars_vep$polyphen_prediction)]<-"unknown"
vars_vep$sift_prediction[is.na(vars_vep$sift_prediction)]<-"unknown"

#remove uppercase letters after "del" in DAMs, but also in intogen (no del mapped)
vars_vep$var[grep("del[[:upper:]]+$",vars_vep$var)]<-gsub("del[[:upper:]]+$","del",vars_vep$var[grep("del[[:upper:]]+$",vars_vep$var)])

#taking the most frequent prediction per AA variant
vars_vep_summ<-data.frame(var=unique(vars_vep$var))
sift_max<-c()
poly_max<-c()
for(v in unique(vars_vep$var)){
  sift<-vars_vep$sift_prediction[which(vars_vep$var==v)]
  sm<-setdiff(names(sort(table(sift))), "unknown")
  sift_max<-c(sift_max, ifelse(length(sm)==0, "unknown", sm))
  poly<-vars_vep$polyphen_prediction[which(vars_vep$var==v)]
  pm<-setdiff(names(sort(table(poly))), "unknown")
  poly_max<-c(poly_max, ifelse(length(pm)==0, "unknown", pm))
}

vars_vep_summ$polyphen_prediction<-poly_max
vars_vep_summ$sift_prediction<-sift_max

############create a summary data.frame
summary_vars_bytiss$driver<-ifelse(summary_vars_bytiss$Gene %in% driver_genes, "Driver", "Non-driver")
summary_vars_bytiss$SIFT<-vars_vep_summ$sift_prediction[match(summary_vars_bytiss$ID, vars_vep_summ$var)]
summary_vars_bytiss$PolyPhen<-vars_vep_summ$polyphen_prediction[match(summary_vars_bytiss$ID, vars_vep_summ$var)]

summary_vars_bytiss$SIFT<-gsub("\\(.*\\)", "", summary_vars_bytiss$SIFT)
summary_vars_bytiss$PolyPhen<-gsub("\\(.*\\)", "", summary_vars_bytiss$PolyPhen)

#######also the initial summary_vars
summary_vars$driver<-ifelse(summary_vars$Gene %in% driver_genes, "Driver", "Non-driver")
summary_vars$SIFT<-vars_vep_summ$sift_prediction[match(summary_vars$ID, vars_vep_summ$var)]
summary_vars$PolyPhen<-vars_vep_summ$polyphen_prediction[match(summary_vars$ID, vars_vep_summ$var)]

summary_vars$SIFT<-gsub("\\(.*\\)", "", summary_vars$SIFT)
summary_vars$PolyPhen<-gsub("\\(.*\\)", "", summary_vars$PolyPhen)




summary_sel<-summary_vars_bytiss[which(summary_vars_bytiss$Perc_tiss>(0)),]
summary_sel<-summary_sel[which(summary_sel$tissue=="ALL_MATCHING"),]
summary_sel$Num_COSMIC<-rep(summary_sel$num_patients[summary_sel$database=="COSMIC"], each=2)
summary_sel$Num_Intogen<-rep(summary_sel$num_patients[summary_sel$database=="Intogen"], each=2)
summary_sel$tot_patients<-rowSums(summary_sel[,c("Num_COSMIC", "Num_Intogen")])
summary_sel<-summary_sel[which(summary_sel$Matching_union>0),]

summary_sel<-summary_sel[summary_sel$tissue=="ALL_MATCHING"&summary_sel$database=="COSMIC",]
summary_sel<-summary_sel[,c("ID", "Gene","Var","driver","Perc_tiss", "Tiss", "tot_patients", "Num_Intogen", "Num_COSMIC",
                            
                            "Num_tiss_Intogen", "Num_tiss_COSMIC",
                            "Type_Intogen", "Type_COSMIC",
                            "Match_tiss_Intogen", "Match_tiss_COSMIC",
                            "Match_tiss_Intogen_name", "Match_tiss_COSMIC_name", "Matching_inboth", "Matching_union",
                            "SIFT", "PolyPhen")]
summary_sel<-summary_sel[order(summary_sel$tot_patients, decreasing = T),]


summary_vars_bytiss$driver<-ifelse(summary_vars_bytiss$Gene %in% driver_genes, "Driver", "Non-driver")
summary_vars_bytiss$SIFT<-vars_vep_summ$sift_prediction[match(summary_vars_bytiss$ID, vars_vep_summ$var)]
summary_vars_bytiss$PolyPhen<-vars_vep_summ$polyphen_prediction[match(summary_vars_bytiss$ID, vars_vep_summ$var)]

summary_vars_bytiss$SIFT<-gsub("\\(.*\\)", "", summary_vars_bytiss$SIFT)
summary_vars_bytiss$PolyPhen<-gsub("\\(.*\\)", "", summary_vars_bytiss$PolyPhen)


Variant<-summary_vars_bytiss[,c("ID", "Gene", "SIFT", "PolyPhen")]
colnames(Variant)<-c("variant_id", "gene_id", "sift", "polyphen")

Variant<-unique(Variant)

#check that the id is unique
sum(duplicated(Variant$variant_id))==0

Variant$sift[is.na(Variant$sift)]<-""
Variant$polyphen[is.na(Variant$polyphen)]<-""

write.csv(Variant, file="DB_tables/Variant.csv", quote=F, row.names = F)

#############
load("all_tiss_comsic_byvar.RData")
load("all_tiss_intogen_byvar.RData")

matching_intogen<-unique(mapping_intogen[which(mapping_intogen[,2] %in% type_intogen),1])
matching_cosmic<-unique(mapping_cosmic[which(mapping_cosmic[,2] %in% type_cosmic),1])

all_cosmic<-matrix(0, nrow=length(all_tiss_cosmic), ncol=length(tissue_names$tissue_name))
rownames(all_cosmic)<-names(all_tiss_cosmic)
colnames(all_cosmic)<-tissue_names$tissue_name

for(i in 1:length(all_tiss_cosmic)){
  for(j in 1:length(names(all_tiss_cosmic[[i]]))){
    tiss<-mapping_cosmic[mapping_cosmic[,2]==names(all_tiss_cosmic[[i]])[j],1]
    if(length(tiss)>0 ){

    all_cosmic[names(all_tiss_cosmic[i]),tiss]<-all_cosmic[names(all_tiss_cosmic[i]),tiss]+unlist(all_tiss_cosmic[[i]][j])
      
    }
  }
}


all_intogen<-matrix(0, nrow=length(all_tiss_intogen), ncol=length(tissue_names$tissue_name))
rownames(all_intogen)<-names(all_tiss_intogen)
colnames(all_intogen)<-tissue_names$tissue_name

for(i in 1:length(all_tiss_intogen)){
  for(j in 1:length(names(all_tiss_intogen[[i]]))){
    tiss<-mapping_intogen[mapping_intogen[,2]==names(all_tiss_intogen[[i]])[j],1]
      if(length(tiss)>0 ){
    
        all_intogen[names(all_tiss_intogen[i]),tiss]<-all_intogen[names(all_tiss_intogen[i]),tiss]+unlist(all_tiss_intogen[[i]][j])
      }
    
  }
}

all_tot<-all_cosmic+all_intogen

TissueHasVariant<-data.frame(variant_id=rep(rownames(all_tot), ncol(all_tot)),
                                    tissue_name=rep(colnames(all_tot), each=nrow(all_tot)),
                                    n_patients=c(all_tot))


#remove uppercase letters after "del" in DAMs, but also in intogen (no del mapped)
vars_tot_adj<-vars_tot
vars_tot_adj[grep("del[[:upper:]]+$",vars_tot_adj)]<-gsub("del[[:upper:]]+$","del",vars_tot_adj[grep("del[[:upper:]]+$",vars_tot_adj)])

TissueHasVariant$dam<-0
TissueHasVariant$dam[which(paste(TissueHasVariant$variant_id, TissueHasVariant$tissue_name) %in% gsub("p\\.", "", paste(vars_tot_adj, tiss_tot)))]<-1

write.csv(TissueHasVariant, file="DB_tables/TissueHasVariant.csv", quote=F, row.names = F)


#############################
##### drugs??
############################

##load dei dati

setwd(paste(home, "/results/20250221/_DR_plots", sep=""))
tissues_drugs<-gsub("_DR_validation.RData", "", list.files(pattern="validation.RData"))

results_drugs<-list()
ind_iter<-0
for(ctiss in tissues_drugs){
  ind_iter<-ind_iter+1
  load(paste(ctiss, "_DR_validation.RData", sep=""))
  results_drugs[[ind_iter]] <- rRES
}
names(results_drugs)<-tissues_drugs


genes_tot<-c()
vars_tot<-c()
tiss_tot<-c()
drugs_tot<-c()
for(i in 1:length(results_drugs)){
  vars_list<-strsplit(results_drugs[[i]]$Var[results_drugs[[i]]$bestRR<1.6], " \\| ")
  genes<-rownames(results_drugs[[i]])[results_drugs[[i]]$bestRR<1.6]
  drugs<-results_drugs[[i]]$bestRR_drug_name[results_drugs[[i]]$bestRR<1.6]
  names(vars_list)<-genes

  genes_tot<-c(genes_tot, unlist(lapply(genes, function(x){rep(x,length(unlist(vars_list[x])))})))
  vars_tot<-c(vars_tot, gsub("p.", "", unlist(vars_list)))
  tiss_tot<-c(tiss_tot, rep(tissues_drugs[i], length(unlist(vars_list))))
  names(vars_list)<-drugs
  drugs_tot<-c(drugs_tot, unlist(lapply(drugs, function(x){rep(x,length(unlist(vars_list[x])))})))

}

summary_drugs<-data.frame(genes=genes_tot, vars=vars_tot, tiss=tiss_tot, drugs_tot)
summary_drugs$var_id<-paste(summary_drugs$genes, summary_drugs$vars, sep="-")


View(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id),])
View(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id & summary_vars_bytiss$tissue=="ALL" & summary_vars_bytiss$database=="COSMIC"),])

print(paste("Number of SAMs that are also reliable DAMs", nrow(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id & summary_vars_bytiss$tissue=="ALL" & summary_vars_bytiss$database=="COSMIC"),])))

print(paste("Number of SAMs that are also reliable DAMs in at least one patient", nrow(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id & summary_vars_bytiss$tissue=="ALL" & summary_vars_bytiss$database=="COSMIC" & summary_vars_bytiss$num_patients>0 ),])))

summary_vars_bytiss$driver<-ifelse(summary_vars_bytiss$Gene %in% driver_genes, "Driver", "Non-driver")
summary_vars_bytiss$SIFT<-vars_vep$sift_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]
summary_vars_bytiss$PolyPhen<-vars_vep$polyphen_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]

summary_vars_bytiss$SIFT<-gsub("\\(.*\\)", "", summary_vars_bytiss$SIFT)
summary_vars_bytiss$PolyPhen<-gsub("\\(.*\\)", "", summary_vars_bytiss$PolyPhen)
summary_vars_bytiss$tot_patients<-summary_vars_bytiss$Num_Intogen+summary_vars_bytiss$Num_COSMIC
summary_vars_bytiss<-summary_vars_bytiss[,c("ID","Gene","Var","driver","Perc_tiss", "Tiss", "tot_patients", "Num_Intogen", "Num_COSMIC",
                                            
                                            "Num_tiss_Intogen", "Num_tiss_COSMIC",
                                            "Type_Intogen", "Type_COSMIC",
                                            "Match_tiss_Intogen", "Match_tiss_COSMIC",
                                            "Match_tiss_Intogen_name", "Match_tiss_COSMIC_name", "Matching_inboth", "Matching_union",
                                            "SIFT", "PolyPhen", "tissue")]
summary_vars_bytiss<-summary_vars_bytiss[order(summary_vars_bytiss$tot_patients, decreasing = T),]

##############
###DRUGS######
#############
genes<-cbind(gene_annot$hgnc_symbol, ifelse(gene_annot$hgnc_symbol %in% driver_genes, 1, 0))
colnames(genes)<-c("gene", "is_driver")
genes<-genes[-which(duplicated(genes[,"gene"])),]

setwd(home)
library(openxlsx)
#load GDSC data
drugTargetInfo <- read.table('data/raw/drug-target_data_hgvs_clean_Goncalves_et_all.txt',sep='\t',stringsAsFactors = FALSE,header=TRUE)
### drug-target_data_hgvs_clean_Goncalves_et_all.txt built from https://www.embopress.org/doi/suppl/10.15252/msb.20199405/suppl_file/msb199405-sup-0003-datasetev2.xlsx on 20241003


drug_targets<-drugTargetInfo[, c("Drug.ID", "Gene.Target")]
colnames(drug_targets)<-c("drug_gdsc_id", "gene_id")
drug_targets<-drug_targets[drug_targets[,2]!="nan",]
drug_targets<-drug_targets[which(drug_targets$gene_id %in% genes[,"gene"]),]
write.csv(drug_targets, "drug_targets_gene.csv", row.names = F, quote = F)

gdsc1<-read.csv('data/raw/GDSC1_fitted_dose_response_27Oct23.csv',header = TRUE,stringsAsFactors = FALSE)
gdsc2<-read.csv('data/raw/GDSC2_fitted_dose_response_27Oct23.csv',header = TRUE,stringsAsFactors = FALSE)
### GDSC1_fitted_dose_response_27Oct23.csv and GDSC2_fitted_dose_response_27Oct23.csv have been downloaded from: 
### https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC1_fitted_dose_response_27Oct23.xlsx and https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx
### on the 20241003

colnames(gdsc1)[1]<-"DATASET"
colnames(gdsc2)[1]<-"DATASET"

gdscAll<-rbind(gdsc1,gdsc2)

drugs<-data.frame(drugTargetInfo$Drug.ID, drugTargetInfo$Name)
drugs<-unique(drugs)
colnames(drugs)<-c("drug_gdsc_id", "drug_name")
drugs<-drugs[drugs[,1] %in% drug_targets[,1],]
write.csv(drugs, "drug.csv", row.names = F, quote = F)

sensitivity<-data.frame(drug_gdsc_id=gdscAll$DRUG_ID, cell_line_id=gdscAll$SANGER_MODEL_ID,
                        gdsc=gdscAll$DATASET, value=gdscAll$LN_IC50, conc_min=gdscAll$MIN_CONC, conc_max=gdscAll$MAX_CONC)
sensitivity<-sensitivity[sensitivity[,1] %in% c(drugs[,2]),]

colnames(sensitivity)<-c("drug_gdsc_id", "cell_line_id", "gdsc_version", "IC50", "conc_min", "conc_max")
sensitivity[, "gdsc_version"]<-gsub("GDSC", "", sensitivity[, "gdsc_version"])
sensitivity<-sensitivity[which(sensitivity$drug_gdsc_id %in% drugs$drug_gdsc_id),]
sensitivity<-sensitivity[which(sensitivity$cell_line_id %in% CellLine$cell_line_id),]

write.csv(sensitivity, "sensitivity_profiles.csv", row.names=F, quote = F)


path_results<-"/results/20250221/_DR_plots"

setwd(paste(home, "/", path_results, sep=""))

##load dei dati
tissues_drugs<-gsub("_DR_validation.RData", "", list.files(pattern="validation.RData"))

results_drugs<-list()
ind_iter<-0
for(ctiss in tissues_drugs){
  ind_iter<-ind_iter+1
  load(paste(ctiss, "_DR_validation.RData", sep=""))
  results_drugs[[ind_iter]] <- rRES
}
names(results_drugs)<-tissues_drugs

mat_tot<-matrix(nrow=1, ncol=5)
colnames(mat_tot)<-c("drug_id", "var_id", "gene", "tissue", "gdsc")
for(i in 1:length(results_drugs)){
  vars_list<-strsplit(results_drugs[[i]]$Var[results_drugs[[i]]$bestRR<1.6], " \\| ")
  genes<-rownames(results_drugs[[i]])[results_drugs[[i]]$bestRR<1.6]
  drugs<-results_drugs[[i]]$bestRR_drug_id[results_drugs[[i]]$bestRR<1.6]
  gd<-results_drugs[[i]]$bestRR_Screen[results_drugs[[i]]$bestRR<1.6]
  names(vars_list)<-genes
  drugs_list<-strsplit(drugs, " \\| ")
  names(drugs_list)<-genes
  gd_list<-strsplit(gd, " \\| ")
  names(gd_list)<-genes

  for(g in genes){
    mat_tmp<-expand.grid(drugs_list[[g]], vars_list[[g]])
    colnames(mat_tmp)<-c("drug_id", "var_id")
    mat_tmp$gene<-g
    mat_tmp$tissue<-tissues_drugs[i]
    mat_tmp$gdsc<-gd_list[[g]][match(mat_tmp$drug_id, drugs_list[[g]])]
    mat_tot<-rbind(mat_tot, mat_tmp)
    }

}

mat_tot$var_id<-gsub("p.", "", unlist(mat_tot$var_id))
mat_tot<-mat_tot[-1,]
write.xlsx(mat_tot, "SAM.xlsx", rowNames=T)

