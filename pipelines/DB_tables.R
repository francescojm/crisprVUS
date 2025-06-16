########################################
##Generate tables for VUS database
#########################################
library(tidyverse)

path_data<-"/data"
path_results<-"/results/20220208"
home<-"/Volumes/home/VUS/VUS"

setwd(paste(home, "/", path_results, sep=""))
dir.create("DB_tables")

##load essentiality data
load(paste(home, "/", path_data,'/R/Sanger_Broad_higQ_bdep.RData', sep=""))
load(paste(home, "/", path_data, '/R/Sanger_Broad_higQ_scaled_depFC.RData', sep=""))

##############################
##cell lines annotations
#############################
CMP_annot <- read_csv(paste(home, "/",path_data,"/model_list_20210611.csv", sep="")) # from https://cog.sanger.ac.uk/cmp/download/model_list_20210611.csv
##select cell lines with essentiality data
CMP_annot<-CMP_annot[which(is.element(CMP_annot$model_name,colnames(scaled_depFC))),]
tissues<-CMP_annot$cancer_type
st<-summary(as.factor(tissues))
tissues<-sort(setdiff(tissues,names(which(st<5))))
tissues<-setdiff(tissues,c('Other Solid Carcinomas','Other Solid Cancers','Other Sarcomas'))
tissues<-sort(tissues)
CMP_annot<-CMP_annot[which(is.element(CMP_annot$cancer_type,tissues)),]
CellLine<-CMP_annot[, c("model_id", "model_name", "cancer_type")]
colnames(CellLine)<-c("cell_line_id", "cell_line_name", "tissue_name")
#check that the id is unique
sum(duplicated(CellLine$CellLineID))==0

write.csv(CellLine, "DB_tables/CellLine.csv", quote=F, row.names = F)


#######################
### Variants
####################

## latest sanger/broad unreleased yet variants hg38
cl_variants <- read_csv(paste(home, "/", path_data,'/mutations_all_latest.csv', sep=""))

##gene annotation
gene_annot <- read_csv(paste(home, "/",path_data, "/gene_identifiers_20191101.csv", sep=""))
cl_variants <- cbind(cl_variants,gene_annot$hgnc_symbol[match(cl_variants$gene_id,gene_annot$gene_id)])
colnames(cl_variants)[ncol(cl_variants)]<-'gene_symbol_2019'
##model names
cl_variants <- cbind(cl_variants,CMP_annot$model_name[match(cl_variants$model_id,CMP_annot$model_id)])
colnames(cl_variants)[ncol(cl_variants)] <- 'model_name'
##cancer type
cl_variants <- cbind(cl_variants,CMP_annot$cancer_type[match(cl_variants$model_id,CMP_annot$model_id)])
colnames(cl_variants)[ncol(cl_variants)] <- 'cancer_type'
colnames(cl_variants)[12] <- 'protein_mutation'


##used for gene filtering
load(paste(home, "/", path_data,'/preComputed/ADaM.RData', sep=""))
load(paste(home, "/", path_data,'/preComputed/FiPer_outputs.RData', sep=""))

##########
### Load DAMs
##############
setwd(paste(home, "/", path_results, sep=""))

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

scaled_depFC_sel<-scaled_depFC[allgenes,unique(CMP_annot$model_name)]
colnames(scaled_depFC_sel)<-CMP_annot$model_id[match(colnames(scaled_depFC_sel), CMP_annot$model_name)]

library(reshape2)
scaled_depFC_sel_melt<-melt(data = scaled_depFC_sel, varnames = c("gene", "line"))
Essentiality<-scaled_depFC_sel_melt
colnames(Essentiality)<-c("gene_id", "cell_line_id", "essentiality_value")

#check that the id is unique
sum(duplicated(paste(Essentiality$gene_id, Essentiality$CellLineID)))==0
write.csv(Essentiality, file="DB_tables/Essentiality.csv", quote=F, row.names = F)

#####################

CellLineHasVariant<-cl_variants[(cl_variants$model_id %in% CellLine$cell_line_id & paste(cl_variants$gene_symbol_2019, cl_variants$protein_mutation, sep="-") %in% allvars),]
CellLineHasVariant$VariantID<-paste(CellLineHasVariant$gene_symbol_2019, CellLineHasVariant$protein_mutation, sep="-")
CellLineHasVariant<-CellLineHasVariant[,c("model_id", "VariantID")]
colnames(CellLineHasVariant)<-c("CellLineID", "VariantID")

#check that the id is unique
sum(duplicated(paste(CellLineHasVariant$CellLineID, CellLineHasVariant$VariantID)))==0
CellLineHasVariant<-unique(CellLineHasVariant)

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
###########################
#### VEP
##############################

vars_vep<-read.csv(paste(home, "/", path_data, "/allvars_id_20220215_hg38.annotated_variants.txt", sep=""), sep="\t", header=T)


#################
#### each cancer type separately
###################
load("summary_byvar.RData")
summary_vars$ID<-paste(summary_vars$Gene, summary_vars$Var, sep="-")
summary_vars$Perc_Intogen<-summary_vars$Perc_Intogen/100
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
library(openxlsx)
cancer_match_ext<-read.xlsx("cancer_match_ext.xlsx",1)
cosmic_match<-read.xlsx("cosmic_match.xlsx",1)

for(i in 1:nrow(summary_vars)){
  type_tiss<-unique(unlist(strsplit(summary_vars[i, "Tiss"], " \\| ")))
  type_intogen<-na.omit(unique(unlist(strsplit(summary_vars[i, "Type_Intogen"], " \\| "))))
  type_cosmic<-na.omit(unique(unlist(strsplit(summary_vars[i, "Type_COSMIC"], " \\| "))))

  summary_vars$Num_tiss_Intogen[i]<-length(type_intogen)
  summary_vars$Num_tiss_COSMIC[i]<-length(type_cosmic)

  matching_intogen<-unique(cancer_match_ext[which(cancer_match_ext[,2] %in% type_intogen),1])

  matching_cosmic<-unique(cosmic_match[which(cosmic_match[,3] %in% type_cosmic),1])

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

#colonna num_patients
#colonna tissue (anche ALL e ALL_MATCHING)
#colonna db
for(i in 1:nrow(summary_vars)){
  type_tiss<-unique(unlist(strsplit(summary_vars[i, "Tiss"], " \\| ")))
  type_intogen<-na.omit(unique(unlist(strsplit(summary_vars[i, "Type_Intogen"], " \\| "))))
  type_cosmic<-na.omit(unique(unlist(strsplit(summary_vars[i, "Type_COSMIC"], " \\| "))))

  matching_intogen<-unique(cancer_match_ext[which(cancer_match_ext[,2] %in% type_intogen),1])

  matching_cosmic<-unique(cosmic_match[which(cosmic_match[,3] %in% type_cosmic),1])

  ###matched tissues cosmic
  if(length(intersect(type_tiss, matching_cosmic))>0){
    num_cosmic_bytiss<-c()
    for(tiss in intersect(type_tiss, matching_cosmic)){
      tiss_cosmic<-cosmic_match[cosmic_match[,1]==tiss,3]
      num_cosmic_bytiss<-c(num_cosmic_bytiss, sum(all_tiss_cosmic[[summary_vars$ID[i]]][tiss_cosmic], na.rm=T))
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
      tiss_intogen<-cancer_match_ext[cancer_match_ext[,1]==tiss,2]
      num_intogen_bytiss<-c(num_intogen_bytiss, sum(all_tiss_intogen[[summary_vars$ID[i]]][tiss_intogen], na.rm=T))
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

#####save list of DAMs  with at least one matching patient

summary_vars_bytiss$SIFT<-vars_vep$transcript_consequences_sift_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]
summary_vars_bytiss$PolyPhen<-vars_vep$transcript_consequences_polyphen_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]

summary_vars_bytiss$SIFT<-gsub("\\(.*\\)", "", summary_vars_bytiss$SIFT)
summary_vars_bytiss$PolyPhen<-gsub("\\(.*\\)", "", summary_vars_bytiss$PolyPhen)

Variant<-summary_vars_bytiss[,c("ID", "Gene", "SIFT", "PolyPhen")]
colnames(Variant)<-c("variant_id", "gene_id", "sift", "polyphen")

###add some missing variants that were not in any patient, so not included in the Variant table
Variant_missing<-setdiff(gsub(".p", "", vars_tot), rownames(summary_vars_bytiss))
gene_missing<-gsub("\\..*","",Variant_missing)
Variant_missing<-gsub("[.]", "-", Variant_missing)

Variant_missing_SIFT<-vars_vep$transcript_consequences_sift_prediction[match(Variant_missing, vars_vep$variant_id)]
Variant_missing_PolyPhen<-vars_vep$transcript_consequences_polyphen_prediction[match(Variant_missing, vars_vep$variant_id)]

Variant_missing_SIFT<-gsub("\\(.*\\)", "", Variant_missing_SIFT)
Variant_missing_PolyPhen<-gsub("\\(.*\\)", "", Variant_missing_PolyPhen)

Variant_missing_df<-data.frame(variant_id=Variant_missing, gene_id=gene_missing,
                               sift=Variant_missing_SIFT, polyphen=Variant_missing_PolyPhen)

Variant<-rbind.data.frame(Variant, Variant_missing)
colnames(Variant)<-c("variant_id", "gene_id", "sift", "polyphen")

Variant<-unique(Variant)

#check that the id is unique
sum(duplicated(Variant$VariantID))==0

Variant$sift[is.na(Variant$sift)]<-"missing"
Variant$polyphen[is.na(Variant$polyphen)]<-"missing"

#togliere le lettere maiuscole subito dopo "del" in variants
Variant$variant_id[grep("del[[:upper:]]+$",Variant$variant_id)]<-gsub("del[[:upper:]]+$","del",Variant$variant_id[grep("del[[:upper:]]+$",Variant$variant_id)])

Variant$is_driver<-ifelse(Variant$gene_id %in% driver_genes, 1, 0)

write.csv(Variant, file="DB_tables/Variant.csv", quote=F, row.names = F)

#############
load("all_tiss_comsic_byvar.RData")
load("all_tiss_intogen_byvar.RData")

cosmic_match<-read.xlsx("cosmic_match.xlsx")
cancer_match_ext<-read.xlsx("cancer_match_ext.xlsx")

matching_intogen<-unique(cancer_match_ext[which(cancer_match_ext[,2] %in% type_intogen),1])
matching_cosmic<-unique(cosmic_match[which(cosmic_match[,3] %in% type_cosmic),1])

all_cosmic<-matrix(0, nrow=length(all_tiss_cosmic), ncol=length(tissues))
rownames(all_cosmic)<-names(all_tiss_cosmic)
colnames(all_cosmic)<-tissues

for(i in 1:length(all_tiss_cosmic)){
  for(j in 1:length(names(all_tiss_cosmic[[i]]))){
    tiss<-cosmic_match[match(names(all_tiss_cosmic[[i]])[j], cosmic_match[,3]),1]
    if(length(tiss)>0 ){
      if(!is.na(tiss)){
    all_cosmic[names(all_tiss_cosmic[i]),tiss]<-all_tiss_cosmic[[i]][j]
      }
    }
  }
}


cancer_match_ext[cancer_match_ext[,1]=="Gastric carcinoma",1]<-"Gastric Carcinoma"
all_intogen<-matrix(0, nrow=length(all_tiss_intogen), ncol=length(tissues))
rownames(all_intogen)<-names(all_tiss_intogen)
colnames(all_intogen)<-tissues

for(i in 1:length(all_tiss_intogen)){
  for(j in 1:length(names(all_tiss_intogen[[i]]))){
    tiss<-cancer_match_ext[match(names(all_tiss_intogen[[i]])[j], cancer_match_ext[,2]),1]
    if(length(tiss)>0 ){
      if(!is.na(tiss)){
        all_intogen[names(all_tiss_intogen[i]),tiss]<-all_tiss_intogen[[i]][j]
      }
    }
  }
}

all_tot<-all_cosmic+all_intogen

TissueHasVariant<-data.frame(variant_id=rep(rownames(all_tot), ncol(all_tot)),
                                    tissue_name=rep(colnames(all_tot), each=nrow(all_tot)),
                                    n_patients=c(all_tot))


TissueHasVariant$dam<-0
TissueHasVariant$dam[which(paste(TissueHasVariant$variant_id, TissueHasVariant$tissue_name) %in% gsub("p\\.", "", paste(vars_tot, tiss_tot)))]<-1

write.csv(TissueHasVariant, file="DB_tables/TissueHasVariant.csv", quote=F, row.names = F)


#############################
##### drugs??
############################

##load dei dati
##driver genes
inTOgen_drivers<-read.table(paste(home, '/data/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
driver_genes<-unique(inTOgen_drivers$SYMBOL)

setwd(paste(home, "/results/20220208/_DR_plots", sep=""))
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

setwd("/Volumes/home/VUS/VUS/results/20220208")

summary_vars_bytiss$driver<-ifelse(summary_vars_bytiss$Gene %in% driver_genes, "Driver", "Non-driver")
summary_vars_bytiss$SIFT<-vars_vep$transcript_consequences_sift_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]
summary_vars_bytiss$PolyPhen<-vars_vep$transcript_consequences_polyphen_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]

summary_vars_bytiss$SIFT<-gsub("\\(.*\\)", "", summary_vars_bytiss$SIFT)
summary_vars_bytiss$PolyPhen<-gsub("\\(.*\\)", "", summary_vars_bytiss$PolyPhen)
summary_vars_bytiss$tot_patients<-summary_vars_bytiss$Num_Intogen+summary_vars_bytiss$Num_COSMIC
summary_vars_bytiss<-summary_vars_bytiss[,c("ID","Gene","Var","driver","Perc_tiss", "Tiss", "tot_patients", "Num_Intogen", "Num_COSMIC",
                                            "Perc_Intogen", "Perc_COSMIC",
                                            "Num_tiss_Intogen", "Num_tiss_COSMIC",
                                            "Type_Intogen", "Type_COSMIC",
                                            "Match_tiss_Intogen", "Match_tiss_COSMIC",
                                            "Match_tiss_Intogen_name", "Match_tiss_COSMIC_name", "Matching_inboth", "Matching_union",
                                            "SIFT", "PolyPhen", "tissue")]
summary_vars_bytiss<-summary_vars_bytiss[order(summary_vars_bytiss$tot_patients, decreasing = T),]

##############
###DRUGS######
#############
setwd(home)

drugTargetInfo <- read.table(paste(home, path_data,'/raw/drug-target_data_hgvs_clean_Goncalves_et_all.txt', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)

drug_targets<-drugTargetInfo[, c("Drug.ID", "Gene.Target")]
colnames(drug_targets)<-c("drug_gdsc_id", "gene_id")
drug_targets<-drug_targets[drug_targets[,2]!="nan",]
write.xlsx(drug_targets, "drug_targets_gene.xlsx", rowNames=F)

gdsc1<-read.csv(paste(home,path_data,'/raw/GDSC1_fitted_dose_response_25Feb20.csv', sep=""),header = TRUE,stringsAsFactors = FALSE)
gdsc2<-read.csv(paste(home,path_data,'/raw/GDSC2_fitted_dose_response_25Feb20.csv', sep=""),header = TRUE,stringsAsFactors = FALSE)
colnames(gdsc1)[1]<-"DATASET"
colnames(gdsc2)[1]<-"DATASET"

gdscAll<-rbind(gdsc1,gdsc2)

drugs<-data.frame(drugTargetInfo$Name, drugTargetInfo$Drug.ID)
drugs<-unique(drugs)
colnames(drugs)<-c("drug_name", "drug_gdsc_id")
drugs<-drugs[drugs[,2] %in% drug_targets[,1],]
write.xlsx(drugs, "drug.xlsx", rowNames=F)

sensitivity<-data.frame(drug_gdsc_id=gdscAll$DRUG_ID, cell_line_id=gdscAll$SANGER_MODEL_ID,
                        gdsc=gdscAll$DATASET, value=gdscAll$LN_IC50)
sensitivity<-sensitivity[sensitivity[,1] %in% c(drugs[,2]),]
write.xlsx(sensitivity, "sensitivity.xlsx", rowNames=F)


path_data<-"/data"
path_results<-"/results/20220208/_DR_plots"
home<-"/Volumes/home/VUS/VUS/"

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

