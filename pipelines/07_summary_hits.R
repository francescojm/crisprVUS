library(reshape2)
library(ggplot2)

home<-"/Volumes/home/VUS/VUS/"
##driver genes
inTOgen_drivers<-read.table(paste(home, '/data/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
driver_genes<-unique(inTOgen_drivers$SYMBOL)


setwd("/Volumes/home/VUS/VUS/results/20220208")
load("hits_nodriver.RData")

#dev.off()

#####################################################################
################################### BY VAR
#########################################################################

#######################
##matching cancer types
##################
path_data<-"/data"
path_results<-"/results/20220208"
home<-"/Volumes/home/VUS/VUS/"

setwd(paste(home, "/", path_results, sep=""))

load("summary_byvar.RData")

##load data
tissues<-gsub("_results_ext.RData", "", list.files(pattern="results_ext.RData"))

results<-list()
ind_iter<-0
for(ctiss in tissues){
  ind_iter<-ind_iter+1
  load(paste(ctiss, "_results_ext.RData", sep=""))
  results[[ind_iter]] <- RESTOT
}
names(results)<-tissues

load("results_with_expr.RData")
##driver genes
inTOgen_drivers<-read.table(paste(home, '/data/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
driver_genes<-unique(inTOgen_drivers$SYMBOL)

library(tidyverse)
gene_annot <- read_csv(paste(home, "/",path_data, "/gene_identifiers_20191101.csv", sep=""))

CMP_annot <- read_csv(paste(home, "/",path_data,"/model_list_20210611.csv", sep="")) # from https://cog.sanger.ac.uk/cmp/download/model_list_20210611.csv

#################### match tissue names for COSMIC
CMP_tissue<-CMP_annot$tissue[match(tissues, CMP_annot$cancer_type)]
cosmic_match<-cbind(CMP_tissue,rep(NA, length(CMP_tissue)))

all_COSMIC_tissues<-unique(unlist(strsplit(summary_vars$Type_COSMIC, " \\| ")))
ind_mat<-which(toupper(cosmic_match[,1]) %in% toupper(gsub("_", " ", all_COSMIC_tissues)))
cosmic_match[ind_mat,2]<-gsub(" ", "_", tolower(cosmic_match[ind_mat,1]))

cosmic_match[cosmic_match[,1]=="Haematopoietic and Lymphoid",2]<-"haematopoietic_and_lymphoid_tissue"
cosmic_match[cosmic_match[,1]=="Esophagus",2]<-"oesophagus"
cosmic_match[5,2]<-"urinary_tract"
cosmic_match[cosmic_match[,1]=="Head and Neck",2]<-"upper_aerodigestive_tract | salivary_gland"
cosmic_match[20,2]<-"skin | eye" #melanoma pu� essere anche uveal melanoma
cosmic_match<-cbind(tissues,cosmic_match)

cosmic_ambigu<-cosmic_match[which(cosmic_match[,3]=="upper_aerodigestive_tract | salivary_gland"),]
cosmic_ambigu[,3]<-"salivary_gland"
cosmic_match<-rbind(cosmic_match, cosmic_ambigu)
cosmic_match[which(cosmic_match[,3]=="upper_aerodigestive_tract | salivary_gland"),3]<-"upper_aerodigestive_tract"

cosmic_ambigu<-cosmic_match[which(cosmic_match[,3]=="skin | eye"),]
cosmic_ambigu[3]<-"eye"
cosmic_match<-rbind(cosmic_match, cosmic_ambigu)
cosmic_match[which(cosmic_match[,3]=="skin | eye"),3]<-"skin"

cosmic_match<-rbind(cosmic_match, c("Esophageal Carcinoma", "Esophagus","upper_aerodigestive_tract"))
cosmic_match<-rbind(cosmic_match, c("Esophageal Squamous Cell Carcinoma", "Esophagus","upper_aerodigestive_tract"))
cosmic_match<-rbind(cosmic_match, c("Head and Neck Carcinoma", "Head and Neck","upper_aerodigestive_tract"))
cosmic_match<-rbind(cosmic_match, c("Oral Cavity Carcinoma", "Head and Neck","upper_aerodigestive_tract"))

write.xlsx(data.frame(cosmic_match), file="cosmic_match.xlsx")

#match tissue names for Intogen
load("/Volumes/home/VUS/VUS/results/20220208/cancer_match.RData")
cancer_match_ext<-matrix(ncol=2, byrow=F, unlist(cancer_match[-c(16,19,23),]))
cancer_match_ambigu<-matrix(ncol=2, c("Glioma", "HGG",
                                      "Glioma", "LGG",
                                      "Kidney Carcinoma", "RCCC",
                                      "Kidney Carcinoma", "RCH",
                                      "Kidney Carcinoma", "RPC",
                                      "Non-Small Cell Lung Carcinoma", "NSCLC",
                                      "Non-Small Cell Lung Carcinoma", "LUAD",
                                      "Esophageal Carcinoma", "HNSC",
                                      "Esophageal Squamous Cell Carcinoma", "HNSC"), byrow = T)
cancer_match_ext<-rbind(cancer_match_ext, cancer_match_ambigu)
cancer_match_add<-matrix(ncol=2, c("Gastric carcinoma", "STAD",
                                   "Gastric carcinoma", "GACA",
                                   "B-Lymphoblastic Leukemia", "CLLSLL",
                                   "Melanoma", "MEL",
                                   "B-Cell Non-Hodgkin's Lymphoma", "DLBCLNOS",
                                   "Hepatocellular Carcinoma","HCC",
                                   "Kidney Carcinoma", "CCRCC",
                                   "Kidney Carcinoma", "RCC",
                                   "Kidney Carcinoma", "CHRCC",
                                   "Kidney Carcinoma", "PRCC",
                                   "Colorectal Carcinoma", "COAD",
                                   "Colorectal Carcinoma", "COADREAD",
                                   "Biliary Tract Carcinoma", "CHOL",
                                   "Pancreatic Carcinoma", "PDAC",
                                   "Pancreatic Carcinoma", "PANCREAS",
                                   "Head and Neck Carcinoma", "NPC",
                                   "Neuroblastoma", "NBL",
                                   "Glioma", "HGGNOS",
                                   "Rhabdomyosarcoma", "RMS",
                                   "Rhabdomyosarcoma", "SOFT_TISSUE"

                                   ), byrow = T)
cancer_match_ext<-rbind(cancer_match_ext, cancer_match_add)
write.xlsx(data.frame(cancer_match_ext), file="cancer_match_ext.xlsx")

###################################



###########################
########### for each variant, indicate the number of tissues in intogen, cosmic
############################
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
}

summary_vars$tot_match<-(rowSums(summary_vars[,c("Match_tiss_COSMIC", "Match_tiss_Intogen")]))
summary_vars$tot_patients<-rowSums(summary_vars[,c("Num_COSMIC", "Num_Intogen")])

#write.csv(summary_vars, "summary_vars.csv")

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

summary_vars_bytiss$driver<-ifelse(summary_vars_bytiss$Gene %in% driver_genes, "Driver", "Non-driver")
summary_vars_bytiss$SIFT<-vars_vep$transcript_consequences_sift_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]
summary_vars_bytiss$PolyPhen<-vars_vep$transcript_consequences_polyphen_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]

summary_vars_bytiss$SIFT<-gsub("\\(.*\\)", "", summary_vars_bytiss$SIFT)
summary_vars_bytiss$PolyPhen<-gsub("\\(.*\\)", "", summary_vars_bytiss$PolyPhen)

summary_sel<-summary_vars_bytiss[which(summary_vars_bytiss$Perc_tiss>(0)),]
summary_sel<-summary_sel[which(summary_sel$tissue=="ALL_MATCHING"),]
summary_sel$Num_COSMIC<-rep(summary_sel$num_patients[summary_sel$database=="COSMIC"], each=2)
summary_sel$Num_Intogen<-rep(summary_sel$num_patients[summary_sel$database=="Intogen"], each=2)
summary_sel$tot_patients<-rowSums(summary_sel[,c("Num_COSMIC", "Num_Intogen")])
summary_sel<-summary_sel[which(summary_sel$Matching_union>0),]

summary_sel<-summary_sel[summary_sel$tissue=="ALL_MATCHING"&summary_sel$database=="COSMIC",]
summary_sel<-summary_sel[,c("ID", "Gene","Var","driver","Perc_tiss", "Tiss", "tot_patients", "Num_Intogen", "Num_COSMIC",
                            "Perc_Intogen", "Perc_COSMIC",
                            "Num_tiss_Intogen", "Num_tiss_COSMIC",
                            "Type_Intogen", "Type_COSMIC",
                            "Match_tiss_Intogen", "Match_tiss_COSMIC",
                            "Match_tiss_Intogen_name", "Match_tiss_COSMIC_name", "Matching_inboth", "Matching_union",
                            "SIFT", "PolyPhen")]
summary_sel<-summary_sel[order(summary_sel$tot_patients, decreasing = T),]

library(openxlsx)
write.xlsx(summary_sel, file="DAMs_patients_match.xlsx")


print(paste("Number of DAMs in at least one patient (any type)", length(which(rowSums(summary_vars[,c("Num_COSMIC","Num_Intogen")], na.rm=T)>0))))

print(paste("Number of non-driver genes with at least one matching patient", length(unique(summary_sel$Gene[summary_sel$driver=="Non-driver"]))))

###Print variants with expr rank ratio ==1

tissues<-names(results)

rr_exp_sel<-c()
rr_prot_sel<-c()
for(i in 1:nrow(summary_sel)){
  gene<-summary_sel$Gene[i]
  var<-summary_sel$Var[i]
  Tiss<-summary_sel$Tiss[i]
  for(ctiss in tissues){
    rr_exp_sel<-c(rr_exp_sel, results_expr[[ctiss]]$rr_expr[intersect(which(results_expr[[ctiss]]$GENE == gene) , grep(var, results_expr[[ctiss]]$var))])
    rr_prot_sel<-c(rr_prot_sel, results_expr[[ctiss]]$rr_prot[intersect(which(results_expr[[ctiss]]$GENE == gene) , grep(var, results_expr[[ctiss]]$var))])
    rr<-results_expr[[ctiss]]$rr_expr[intersect(which(results_expr[[ctiss]]$GENE == gene) , grep(var, results_expr[[ctiss]]$var))]
    tiss<-ctiss
    if(length(na.omit(rr))>0 ){
      if(length(grep(tiss, Tiss))>0){
        if(rr==1){
          print(gene)
          print(var)
          print(ctiss)
        }
      }
    }
  }
}

#############################
##### Circular plot
#############################

setwd("/Volumes/home/VUS/VUS/results/20220208")
summary_vars_bytiss$driver<-ifelse(summary_vars_bytiss$Gene %in% driver_genes, "Driver", "Non-driver")
summary_vars_bytiss$SIFT<-vars_vep$transcript_consequences_sift_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]
summary_vars_bytiss$PolyPhen<-vars_vep$transcript_consequences_polyphen_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]

summary_vars_bytiss$SIFT<-gsub("\\(.*\\)", "", summary_vars_bytiss$SIFT)
summary_vars_bytiss$PolyPhen<-gsub("\\(.*\\)", "", summary_vars_bytiss$PolyPhen)

write.xlsx(summary_vars_bytiss, file="Full_var_list.xlsx")

print(paste("Number of DAMs in at least one patient (same type) AND predicted deleterious", nrow(summary_vars_bytiss[summary_vars_bytiss$database=="COSMIC" & summary_vars_bytiss$tissue=="ALL" & ((summary_vars_bytiss$SIFT %in% c("deleterious", "deleterious_low_confidence"))|(summary_vars_bytiss$PolyPhen %in% c("possibly_damaging", "probably_damaging"))),])))

summary_vars_sel<-summary_vars_bytiss[which(summary_vars_bytiss$Perc_tiss>(0)),]
summary_vars_sel<-summary_vars_sel[which(summary_vars_sel$tissue=="ALL_MATCHING"),]
summary_vars_sel$Num_COSMIC<-rep(summary_vars_sel$num_patients[summary_vars_sel$database=="COSMIC"], each=2)
summary_vars_sel$Num_Intogen<-rep(summary_vars_sel$num_patients[summary_vars_sel$database=="Intogen"], each=2)
summary_vars_sel$tot_patients<-rowSums(summary_vars_sel[,c("Num_COSMIC", "Num_Intogen")])

summary_vars_sel$tot_match<-(-(summary_vars_sel[,c("Matching_union")]))
summary_vars_sel$driver<-ifelse(summary_vars_sel$Gene %in% driver_genes, "Driver", "Non-driver")

summary_vars_sel<-summary_vars_sel[which(summary_vars_sel$tot_match<0),]

print(paste("Number of DAMs in at least one patient (same type)", nrow(summary_vars_sel[summary_vars_sel$database=="COSMIC",])))

print(paste("Number of DAMs in at least one patient (same type) AND predicted deleterious", nrow(summary_vars_sel[summary_vars_sel$database=="COSMIC" & ((summary_vars_sel$SIFT %in% c("deleterious", "deleterious_low_confidence"))|(summary_vars_sel$PolyPhen %in% c("possibly_damaging", "probably_damaging"))),])))

ord<-order(summary_vars_sel[,"tot_patients"], decreasing=T)
summary_vars_sel$ID<-factor(summary_vars_sel$ID, levels=unique(summary_vars_sel$ID[ord]))

summary_vars_sel$SIFT<-vars_vep$transcript_consequences_sift_prediction[match(summary_vars_sel$ID, vars_vep$variant_id)]
summary_vars_sel$PolyPhen<-vars_vep$transcript_consequences_polyphen_prediction[match(summary_vars_sel$ID, vars_vep$variant_id)]

summary_vars_sel$SIFT<-gsub("\\(.*\\)", "", summary_vars_sel$SIFT)
summary_vars_sel$PolyPhen<-gsub("\\(.*\\)", "", summary_vars_sel$PolyPhen)

# ----- This section prepare a dataframe for labels ---- #
# Get the name and the y position of each label
summary_vars_sel2<-summary_vars_sel[which(summary_vars_sel$database=="COSMIC"),]

label_data <- summary_vars_sel[ord,]
label_data <-label_data[which(label_data$database=="COSMIC"),]
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  - 360 * (1:nrow(label_data)) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
angle <- angle - max(angle)
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90 & angle > -270 , 0, 1)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90 & angle > -270, angle+180, angle)
# ----- ------------------------------------------- ---- #
label_data$effect_sift<-ifelse(label_data$SIFT %in% c("deleterious", "deleterious_low_confidence"),1,NA)
label_data$effect_polyphen<-ifelse(label_data$PolyPhen %in% c("possibly_damaging", "probably_damaging"),1,NA)




#####sacked geom_col
#num patients transformed so that the sum of logarithms is the logarithm of tot patients 

summary_vars_sel$Num_Intogen_transf<-exp(log(summary_vars_sel$Num_Intogen+1)*log(summary_vars_sel$tot_patients+1)/(log(summary_vars_sel$Num_Intogen+1)+log(summary_vars_sel$Num_COSMIC+1)))-1
summary_vars_sel$Num_COSMIC_transf<-exp(log(summary_vars_sel$Num_COSMIC+1)*log(summary_vars_sel$tot_patients+1)/(log(summary_vars_sel$Num_Intogen+1)+log(summary_vars_sel$Num_COSMIC+1)))-1

melt_var<-summary_vars_sel
melt_var$type_patients<-melt_var$database
melt_var$num_patients[summary_vars_sel$database=="Intogen"]<-summary_vars_sel$Num_Intogen_transf[summary_vars_sel$database=="Intogen"]
melt_var$num_patients[summary_vars_sel$database=="COSMIC"]<-summary_vars_sel$Num_COSMIC_transf[summary_vars_sel$database=="COSMIC"]


melt_vep<-rbind.data.frame(summary_vars_sel2, summary_vars_sel2)
effect_sift<-rep(NA, length(label_data$SIFT))
effect_sift[which(summary_vars_sel2$SIFT %in% c("deleterious", "deleterious_low_confidence"))]<-1
effect_polyphen<-rep(NA, length(label_data$PolyPhen))
effect_polyphen[which(summary_vars_sel2$PolyPhen %in% c("possibly_damaging", "probably_damaging"))]<-1
melt_vep$effect<-c(effect_sift, effect_polyphen)
melt_vep$type_vep<-c(rep("SIFT", nrow(summary_vars_sel2)), rep("PolyPhen", nrow(summary_vars_sel2)))
melt_vep$hjust<-ifelse(melt_vep$type_vep=="SIFT",2, 3)


p <- ggplot(melt_var, aes(x=ID)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  
  geom_col(aes(y = log(num_patients+1), fill = "#F58C07", alpha=type_patients))+
  geom_col(data=label_data, aes(y = tot_match, fill = "#9232DB")) +
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  # This makes the coordinate polar instead of cartesian.
  scale_alpha_discrete(range = c(0.4, 0.8))+
  coord_polar(start = -pi/2, direction=1) +
  #ylim(-100,120) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) +
  
  
  
  scale_fill_identity(name = "",
                      breaks = c("#F58C07", "#9232DB"),
                      labels = c("# of Patients (log)", "# of Cancer Types as Hit"),
                      guide = "legend")+
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=ID, y=log(tot_patients+1)+0.5, label=tot_patients, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_point(data=melt_vep, aes(x=ID, y=log(tot_patients+1)+log10(tot_patients+1)/5 +hjust, size=effect, shape=type_vep), colour="#DB0952")+
  #geom_point(data=label_data, aes(x=ID, y=log(Num_Intogen+1)+log(Num_COSMIC+1)+log10(Num_Intogen+1)/5 +log10(Num_COSMIC+1)/5+3, size=effect_polyphen), fill="#DB0952", shape=24)+
  geom_text(data=label_data, aes(x=ID, y=tot_match-0.5, label= -round(tot_match), hjust=ifelse(hjust==1, 0, 1)), color="black", fontface="bold",alpha=0.6, size=1.5, angle= label_data$angle, inherit.aes = FALSE )+
  geom_text(data=label_data, aes(x=ID, y=log(tot_patients+1)+log10(tot_patients+1)/5 +4, label=ID, hjust=hjust, color=driver), fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  scale_y_continuous(limits = c(-15, 20))+ guides(alpha=guide_legend(title=""))+
  scale_size_area(max_size=1, guide="none", name="", labels="deleterious variants")+
  scale_shape_manual(name="", values=c("SIFT"=16, "PolyPhen"=17), labels=c("Deleterious variant (SIFT)", "Deleterious variant (PolyPhen)"))+
  scale_color_manual(name = "", values=c("Driver"="#DB0952", "Non-driver"="black"))


png("Circ_plot_byvar_stacked_log_vep_matched.png", res=600, 6000, 6000)
p+theme(legend.key.size = unit(3, 'mm'),
        legend.spacing.y = unit(-1, "mm"),
        legend.position = c(.52, .54),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.text = element_text(size=5),
        legend.margin = margin(0, 0, 0, 0)
)
dev.off()

pdf("Circ_plot_byvar_stacked_log_vep_matched.pdf", 10,10)
p+theme(legend.key.size = unit(3, 'mm'),
        legend.spacing.y = unit(-1, "mm"),
        legend.position = c(.52, .54),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.text = element_text(size=5),
        legend.margin = margin(0, 0, 0, 0)
)
dev.off()

##################show matches with drugs 
melt_var$drugs<-ifelse(melt_var$ID %in% summary_drugs$var_id, 1, 0)


p <- ggplot(melt_var, aes(x=ID)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_col(data=melt_var_match, aes(y = -num_match, fill = "#9232DB", alpha=type_match)) +
  geom_col(aes(y = log(num_patients+1), fill = "#F58C07", alpha=type_patients))+
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  # This makes the coordinate polar instead of cartesian.
  scale_alpha_discrete(range = c(0.4, 0.8))+
  coord_polar(start = -pi/2, direction=1) +
  #ylim(-100,120) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) +
  
  
  
  scale_fill_identity(name = "",
                      breaks = c("#F58C07", "#9232DB"),
                      labels = c("# of Patients (log)", "# of Cancer Types as Hit"),
                      guide = "legend")+
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=ID, y=log(tot_patients+1)+0.5, label=tot_patients, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_point(data=melt_var, aes(x=ID, y=log(tot_patients+1)+log10(tot_patients+1)/5 +3, size=drugs), colour="#DB0952")+
  #geom_point(data=label_data, aes(x=ID, y=log(Num_Intogen+1)+log(Num_COSMIC+1)+log10(Num_Intogen+1)/5 +log10(Num_COSMIC+1)/5+3, size=effect_polyphen), fill="#DB0952", shape=24)+
  geom_text(data=label_data, aes(x=ID, y=tot_match-0.5, label=ifelse(tot_match<(-1), -round(tot_match), ""), hjust=ifelse(hjust==1, 0, 1)), color="black", fontface="bold",alpha=0.6, size=2, angle= label_data$angle, inherit.aes = FALSE )+
  geom_text(data=label_data, aes(x=ID, y=log(tot_patients+1)+log10(tot_patients+1)/5 +4, label=ID, hjust=hjust, color=driver), fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  scale_y_continuous(limits = c(-15, 20))+ guides(alpha=guide_legend(title=""))+
  scale_size_area(max_size=1, guide="none", name="", labels="Impact on drugs")+
  scale_color_manual(name = "", values=c("Driver"="#DB0952", "Non-driver"="black"))


png("Circ_plot_byvar_stacked_log_drugs_matched.png", res=600, 5000, 5000)
p+theme(legend.key.size = unit(3, 'mm'),
        legend.spacing.y = unit(-1, "mm"),
        legend.position = c(.52, .54),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.text = element_text(size=5),
        legend.margin = margin(0, 0, 0, 0)
)
dev.off()

pdf("Circ_plot_byvar_stacked_log_drugs_matched.pdf", 10, 10)
p+theme(legend.key.size = unit(3, 'mm'),
        legend.spacing.y = unit(-1, "mm"),
        legend.position = c(.52, .54),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.text = element_text(size=5),
        legend.margin = margin(0, 0, 0, 0)
)
dev.off()

#####################
### Heatmap of variants by tissue
###################
summary_vars_sel<-summary_vars_bytiss[which(summary_vars_bytiss$Perc_tiss>(0)),]
summary_vars_sel2<-summary_vars_sel[which(summary_vars_sel$tissue %in% c("ALL_MATCHING")),]
summary_vars_sel2$Num_COSMIC<-rep(summary_vars_sel2$num_patients[summary_vars_sel2$database=="COSMIC"], each=2)
summary_vars_sel2$Num_Intogen<-rep(summary_vars_sel2$num_patients[summary_vars_sel2$database=="Intogen"], each=2)
summary_vars_sel2$tot_patients<-rowSums(summary_vars_sel2[,c("Num_COSMIC", "Num_Intogen")])
summary_vars_sel2$tot_match<-(-rowSums(summary_vars_sel2[,c("Match_tiss_COSMIC", "Match_tiss_Intogen")]))
summary_vars_sel2<-summary_vars_sel2[which(summary_vars_sel2$tot_match<0),]

summary_vars_sel<-summary_vars_sel[-which(summary_vars_sel$tissue %in% c("ALL","ALL_MATCHING")),]
summary_vars_sel$tot_match<-(-rowSums(summary_vars_sel[,c("Match_tiss_COSMIC", "Match_tiss_Intogen")]))
summary_vars_sel<-summary_vars_sel[which(summary_vars_sel$tot_match<0),]

ID_ord<-unique(summary_vars_sel2$ID[order(summary_vars_sel2$tot_patients,decreasing = T)])
forheat<-matrix(nrow=length(unique(summary_vars_sel$ID)), ncol=33)
colnames(forheat)<-tissues
rownames(forheat)<-unique(summary_vars_sel$ID)
for(ID in unique(summary_vars_sel$ID)){
  for(tiss in tissues){
    forheat[ID, tiss]<-sum(summary_vars_sel$num_patients[summary_vars_sel$ID==ID & summary_vars_sel$tissue==tiss], na.rm=T)
  }
}

library(pheatmap)
forheat<-forheat[,colSums(forheat)>0]
png("heatmap_vars_bytiss.png", res=300, 3000, 8000)
pheatmap(log(forheat[ID_ord,]+1), cluster_rows = F, cellwidth = 10, cellheight = 10)
dev.off()


