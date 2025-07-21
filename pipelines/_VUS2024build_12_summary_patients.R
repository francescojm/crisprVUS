library(reshape2)
library(ggplot2)
library(openxlsx)

#######################
##load data
##################
home<-"E:/crisprVUS"
path_data<-"/data"
path_results<-"/results/20250221"

##driver genes
inTOgen_drivers<-read.table(paste(home, '/data/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
### inTOgene drivers downloaded from https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip on 20241002
driver_genes<-unique(inTOgen_drivers$SYMBOL)

setwd(paste(home, "/", path_results, sep=""))
load("summary_byvar.RData")
load("all_tiss_comsic_byvar.RData")
load("all_tiss_intogen_byvar.RData")

tissues<-gsub("_results_ext.RData", "", list.files(pattern="results_ext.RData"))

results<-list()
ind_iter<-0
for(ctiss in tissues){
  ind_iter<-ind_iter+1
  load(paste(ctiss, "_results_ext.RData", sep=""))
  results[[ind_iter]] <- RESTOT
}
names(results)<-tissues

mapping<-read.csv(paste(home, path_data,'/raw/intOGen ctype mapping_AS.csv',sep=''),header = TRUE,row.names = 1, sep=";")

#reshape cancer-type mapping
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

###########################
########### for each variant, indicate the number of tissues in intogen, cosmic, each cancer type separately
############################

summary_vars$ID<-paste(summary_vars$Gene, summary_vars$Var, sep="-")
summary_vars$Num_tiss<-summary_vars$Perc_tiss*36
summary_vars$Num_Intogen[is.na(summary_vars$Num_Intogen)]<-0

summary_vars$Num_tiss_Intogen<-NA
summary_vars$Num_tiss_COSMIC<-NA
summary_vars$Match_tiss_Intogen<-NA
summary_vars$Match_tiss_COSMIC<-NA
summary_vars$Matching_inboth<-NA
summary_vars$Matching_union<-NA
summary_vars$Match_tiss_Intogen_name<-NA
summary_vars$Match_tiss_COSMIC_name<-NA

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

####include the number of patients separately for intOGen and COSMIC

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

#####save list of DAMs  with patient anno
write.xlsx(summary_vars_bytiss, file=paste(home, path_results, "/summary_DAMs_bytiss.xlsx", sep=""))


##select patients with DAMs with any cancer type (also non-matching)
summary_vars_all<-summary_vars_bytiss[summary_vars_bytiss$tissue=="ALL",c("Gene", "Var", "ID", "Num_Intogen","Num_COSMIC", "num_patients", "database")]

print(paste("Number of patients with a DAM, any cancer type:", sum(summary_vars_all$num_patients)))

ord_IDs<-order(summary_vars_all$Num_Intogen+summary_vars_all$Num_COSMIC, decreasing=T)
ordered_IDs<-unique(summary_vars_all$ID[ord_IDs])
summary_vars_all$ID<-factor(summary_vars_all$ID, levels=c(unique(summary_vars_all$ID[ord_IDs])))
summary_vars_all<-summary_vars_all[ord_IDs,]

ordered_IDs_known<-unique(summary_vars_all$ID[which(summary_vars_all$Gene %in% driver_genes)])
ordered_IDs_unreported<-unique(summary_vars_all$ID[-which(summary_vars_all$Gene %in% driver_genes)])

summary_vars_all_known<-summary_vars_all[which(summary_vars_all$Gene %in% driver_genes),]
summary_vars_all_unreported<-summary_vars_all[-which(summary_vars_all$Gene %in% driver_genes),]

pdf(paste(home, path_results, "/top_DAMs_all_known.pdf", sep=""), 7,6)
ggplot(subset(summary_vars_all_known, ID %in% ordered_IDs_known[1:20]), aes(x=ID, y=num_patients))+geom_bar(stat="identity")+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste(home, path_results, "/top_DAMs_all_unreported.pdf", sep=""), 7,6)
ggplot(subset(summary_vars_all_unreported, ID %in% ordered_IDs_unreported[1:20]), aes(x=ID, y=num_patients))+geom_bar(stat="identity")+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##select patients with DAMs with any matching cancer type 
summary_vars_matching<-summary_vars_bytiss[summary_vars_bytiss$tissue=="ALL_MATCHING",c("Gene", "Var", "ID", "Num_Intogen","Num_COSMIC", "num_patients", "database")]

print(paste("Number of patients with a DAM, matching cancer type:", sum(summary_vars_matching$num_patients)))

ordered_IDs<-names(sort(unlist(by(summary_vars_matching, summary_vars_matching$ID,function(x){sum(x$num_patients)})), decreasing=T))
summary_vars_matching$ID<-factor(summary_vars_matching$ID, levels=c(ordered_IDs))

ordered_IDs_known<-unique(summary_vars_matching$ID[match(ordered_IDs, summary_vars_matching$ID)][which(summary_vars_matching[match(ordered_IDs, summary_vars_matching$ID),"Gene"] %in% driver_genes)])
ordered_IDs_unreported<-unique(summary_vars_matching$ID[match(ordered_IDs, summary_vars_matching$ID)][-which(summary_vars_matching[match(ordered_IDs, summary_vars_matching$ID),"Gene"] %in% driver_genes)])

summary_vars_tiss<-summary_vars_bytiss[-which(summary_vars_bytiss$tissue %in% c("ALL_MATCHING","ALL")),c("Gene", "Var", "ID", "Num_Intogen","Num_COSMIC","tissue", "num_patients", "database")]
summary_vars_tiss$ID<-factor(summary_vars_tiss$ID, levels=c(ordered_IDs))

summary_vars_tiss_known<-summary_vars_tiss[which(summary_vars_tiss$Gene %in% driver_genes),]
summary_vars_tiss_unreported<-summary_vars_tiss[-which(summary_vars_tiss$Gene %in% driver_genes),]

CL_colors<-read.xlsx(paste(home, path_data, "/raw/CL_tissue_ctype_colors.xlsx", sep=""), sheet = 2)
CL_colors_v<-CL_colors$Color_hex_code
names(CL_colors_v)<-CL_colors$Cancer.Type

pdf(paste(home, path_results, "/top_DAMs_matching_known.pdf", sep=""), 10,6)
ggplot(subset(summary_vars_tiss_known, ID %in% ordered_IDs_known[1:20]), aes(x=ID, y=num_patients, fill=tissue))+geom_bar(stat="identity")+scale_fill_manual(values=CL_colors_v)+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste(home, path_results, "/top_DAMs_matching_unreported.pdf", sep=""), 8,6)
ggplot(subset(summary_vars_tiss_unreported, ID %in% ordered_IDs_unreported[1:24]), aes(x=ID, y=num_patients, fill=tissue))+geom_bar(stat="identity")+scale_fill_manual(values=CL_colors_v)+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


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


write.xlsx(summary_sel, file=paste(home, path_results, "/DAMs_patients_match.xlsx", sep=""))


print(paste("Number of DAMs in at least one patient (any type)", length(which(rowSums(summary_vars[,c("Num_COSMIC","Num_Intogen")], na.rm=T)>0))))
print(paste("Number of DAMs in at least one patient (same type)", length(which(summary_vars[,c("Matching_union")]>0))))

print(paste("Number of DAMbgs with a DAM in at least one patient (any type)", length(unique(summary_vars$Gene[which(rowSums(summary_vars[,c("Num_COSMIC","Num_Intogen")], na.rm=T)>0)]))))
print(paste("Number of DAMbgs with a DAM in at least one patient (same type)", length(unique(summary_vars$Gene[which(summary_vars[,c("Matching_union")]>0)]))))

print(paste("Number of non-driver genes with at least one matching patient", length(unique(summary_sel$Gene[summary_sel$driver=="Non-driver"]))))


#############################
##### Circular plot
#############################

summary_vars_bytiss$driver<-ifelse(summary_vars_bytiss$Gene %in% driver_genes, "Driver", "Non-driver")
summary_vars_bytiss$SIFT<-vars_vep_summ$sift_prediction[match(summary_vars_bytiss$ID, vars_vep_summ$var)]
summary_vars_bytiss$PolyPhen<-vars_vep_summ$polyphen_prediction[match(summary_vars_bytiss$ID, vars_vep_summ$var)]

summary_vars_bytiss$SIFT<-gsub("\\(.*\\)", "", summary_vars_bytiss$SIFT)
summary_vars_bytiss$PolyPhen<-gsub("\\(.*\\)", "", summary_vars_bytiss$PolyPhen)

#write.xlsx(summary_vars_bytiss, file="Full_var_list.xlsx")

print(paste("Number of DAMs in at least one patient (same type) AND predicted deleterious", length(which(summary_vars[,c("Matching_union")]>0 & ((summary_vars$SIFT %in% c("deleterious", "deleterious_low_confidence"))|(summary_vars$PolyPhen %in% c("possibly_damaging", "probably_damaging")))))))

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

summary_vars_sel$SIFT<-vars_vep_summ$sift_prediction[match(summary_vars_sel$ID, vars_vep_summ$var)]
summary_vars_sel$PolyPhen<-vars_vep_summ$polyphen_prediction[match(summary_vars_sel$ID, vars_vep_summ$var)]

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
melt_vep<-melt_vep[-which(is.na(melt_vep$effect)),]

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
  scale_shape_manual(name="", values=c("SIFT"=16, "PolyPhen"=17),breaks=c("SIFT", "PolyPhen"), labels=c("Deleterious variant (SIFT)", "Damaging variant (PolyPhen)"))+
  scale_color_manual(name = "", values=c("Driver"="#DB0952", "Non-driver"="black"))


png("Circ_plot_byvar_stacked_log_vep_matched.png", res=600, 10000, 10000)
p+theme(legend.key.size = unit(3, 'mm'),
        legend.spacing.y = unit(-1, "mm"),
        legend.position = c(.52, .54),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.text = element_text(size=5),
        legend.margin = margin(0, 0, 0, 0)
)
dev.off()

pdf("Circ_plot_byvar_stacked_log_vep_matched.pdf", 18,18)
p+theme(legend.key.size = unit(3, 'mm'),
        legend.spacing.y = unit(-1, "mm"),
        legend.position = c(.52, .54),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.text = element_text(size=5),
        legend.margin = margin(0, 0, 0, 0)
)
dev.off()

######cutting at 2 patients
summary_vars_sel<-summary_vars_sel[summary_vars_sel[,"tot_patients"]>=2,]

ord<-order(summary_vars_sel[,"tot_patients"], decreasing=T)

summary_vars_sel$SIFT<-vars_vep_summ$sift_prediction[match(summary_vars_sel$ID, vars_vep_summ$var)]
summary_vars_sel$PolyPhen<-vars_vep_summ$polyphen_prediction[match(summary_vars_sel$ID, vars_vep_summ$var)]

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
  geom_text(data=label_data, aes(x=ID, y=log(tot_patients+1)+0.5, label=tot_patients, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_point(data=melt_vep, aes(x=ID, y=log(tot_patients+1)+log10(tot_patients+1)/5 +hjust, size=effect, shape=type_vep), colour="#DB0952")+
  #geom_point(data=label_data, aes(x=ID, y=log(Num_Intogen+1)+log(Num_COSMIC+1)+log10(Num_Intogen+1)/5 +log10(Num_COSMIC+1)/5+3, size=effect_polyphen), fill="#DB0952", shape=24)+
  geom_text(data=label_data, aes(x=ID, y=tot_match-0.5, label= -round(tot_match), hjust=ifelse(hjust==1, 0, 1)), color="black", fontface="bold",alpha=0.6, size=2, angle= label_data$angle, inherit.aes = FALSE )+
  geom_text(data=label_data, aes(x=ID, y=log(tot_patients+1)+log10(tot_patients+1)/5 +4, label=ID, hjust=hjust, color=driver), fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +
  scale_y_continuous(limits = c(-15, 20))+ guides(alpha=guide_legend(title=""))+
  scale_size_area(max_size=1, guide="none", name="", labels="deleterious variants")+
  scale_shape_manual(name="", values=c("SIFT"=16, "PolyPhen"=17), breaks=c("SIFT", "PolyPhen"), labels=c("Deleterious variant (SIFT)", "Damaging variant (PolyPhen)"))+
  scale_color_manual(name = "", values=c("Driver"="#DB0952", "Non-driver"="black"))

png("Circ_plot_byvar_stacked_log_vep_matched_thr2.png", res=600, 5000, 5000)
p+theme(legend.key.size = unit(5, 'mm'),
        legend.spacing.y = unit(-1, "mm"),
        legend.position = c(.48, .56),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.text = element_text(size=8),
        legend.margin = margin(0, 0, 0, 0)
)
dev.off()

pdf("Circ_plot_byvar_stacked_log_vep_matched_thr2.pdf", 14,14)
p+theme(legend.key.size = unit(5, 'mm'),
        legend.spacing.y = unit(-1, "mm"),
        legend.position = c(.48, .56),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.text = element_text(size=8),
        legend.margin = margin(0, 0, 0, 0)
)
dev.off()


##################show matches with drugs 
load(file=paste(home, path_results, "/_DR_plots/summary_drugs.RData", sep=""))
write.csv(summary_drugs, file="summary_drugs.csv", row.names = F, quote = F)

print(paste("Number of SAMs that are also reliable DAMs", nrow(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id & summary_vars_bytiss$tissue=="ALL" & summary_vars_bytiss$database=="COSMIC"),])))

print(paste("Number of SAMs that are also reliable DAMs in at least one patient", nrow(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id & summary_vars_bytiss$tissue=="ALL" & summary_vars_bytiss$database=="COSMIC" & summary_vars_bytiss$num_patients>0 ),])))

write.xlsx(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id),], paste(home, path_results, "/vars_depANDdrug.xlsx", sep=""))

print(paste("Number of SAMs that are also reliable DAMs in at least one patient (same type)", nrow(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id & summary_vars_bytiss$tissue=="ALL_MATCHING" & summary_vars_bytiss$database=="COSMIC" & summary_vars_bytiss$num_patients>0 ),])))

print(paste("Number of SAMs that are also reliable DAMs in at least one patient (same type) and predicted deleterious", 
            nrow(summary_vars_bytiss[which(summary_vars_bytiss$ID %in% summary_drugs$var_id & summary_vars_bytiss$tissue=="ALL_MATCHING" & summary_vars_bytiss$database=="COSMIC" & summary_vars_bytiss$num_patients>0 &
                                             ((summary_vars_bytiss$SIFT %in% c("deleterious", "deleterious_low_confidence"))|(summary_vars_bytiss$PolyPhen %in% c("possibly_damaging", "probably_damaging") ))),])))

##SAMs in non-driver genes
summary_drugs[-which(summary_drugs$genes %in% driver_genes),]
##SAMs in driver genes
unique(summary_drugs$genes[which(summary_drugs$genes %in% driver_genes)])


melt_var$drugs<-ifelse(melt_var$ID %in% summary_drugs$var_id, 1, 0)



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
forheat<-matrix(nrow=length(unique(summary_vars_sel$ID)), ncol=length(tissues))
colnames(forheat)<-tissues
rownames(forheat)<-unique(summary_vars_sel$ID)
for(ID in unique(summary_vars_sel$ID)){
  for(tiss in tissues){
    forheat[ID, tiss]<-sum(summary_vars_sel$num_patients[summary_vars_sel$ID==ID & summary_vars_sel$tissue==tiss], na.rm=T)
  }
}

library(pheatmap)
forheat<-forheat[,colSums(forheat)>0]
png("heatmap_vars_bytiss.png", res=300, 3000, 20000)
pheatmap(log(forheat[ID_ord,]+1), cluster_rows = F, cellwidth = 10, cellheight = 10)
dev.off()


