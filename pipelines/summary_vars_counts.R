#Heatmap with top variants, number of patients in corresponding tumor types

library(openxlsx)
setwd("/Volumes/home/VUS/VUS/results/20220208")
summary_vars_bytiss$driver<-ifelse(summary_vars_bytiss$Gene %in% driver_genes, "Driver", "Non-driver")
summary_vars_bytiss$SIFT<-vars_vep$transcript_consequences_sift_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]
summary_vars_bytiss$PolyPhen<-vars_vep$transcript_consequences_polyphen_prediction[match(summary_vars_bytiss$ID, vars_vep$variant_id)]

summary_vars_bytiss$SIFT<-gsub("\\(.*\\)", "", summary_vars_bytiss$SIFT)
summary_vars_bytiss$PolyPhen<-gsub("\\(.*\\)", "", summary_vars_bytiss$PolyPhen)

summary_vars_sel<-summary_vars_bytiss[which(summary_vars_bytiss$Perc_tiss>(0)),]
summary_vars_sel<-summary_vars_sel[,-c(22, 23)]
summary_vars_sel<-unique(summary_vars_sel)
summary_vars_sel$tot_patients<-summary_vars_sel$Num_Intogen+summary_vars_sel$Num_COSMIC

summary_vars_sel<-summary_vars_sel[,c("ID", "Gene","Var","driver","Perc_tiss", "Tiss", "tot_patients", "Num_Intogen", "Num_COSMIC",
                            "Perc_Intogen", "Perc_COSMIC",
                            "Num_tiss_Intogen", "Num_tiss_COSMIC", 
                            "Type_Intogen", "Type_COSMIC",
                            "Match_tiss_Intogen", "Match_tiss_COSMIC",
                            "Match_tiss_Intogen_name", "Match_tiss_COSMIC_name", "Matching_inboth", "Matching_union",
                            "SIFT", "PolyPhen", "tissue")]
summary_vars_sel<-summary_vars_sel[order(summary_vars_sel$tot_patients, decreasing = T),]

#write.xlsx(summary_vars_sel, "Full_var_list.xlsx")

#####First test of web interface
summary_vars_sel<-summary_vars_bytiss[summary_vars_bytiss$Gene=="KRAS",]

count_mat_KRAS<-matrix(0, ncol=35, nrow=length(unique(summary_vars_sel$ID)))
rownames(count_mat_KRAS)<-unique(summary_vars_sel$ID)
colnames(count_mat_KRAS)<-c("ALL", "ALL_MATCHING", tissues)
for(var in unique(summary_vars_sel$ID)){
  counts_COSMIC<-summary_vars_sel$num_patients[summary_vars_sel$ID==var & summary_vars_sel$database=="COSMIC"]
  counts_Intogen<-summary_vars_sel$num_patients[summary_vars_sel$ID==var & summary_vars_sel$database=="Intogen"]
  tiss_COSMIC<-summary_vars_sel$tissue[summary_vars_sel$ID==var & summary_vars_sel$database=="COSMIC"]
  tiss_Intogen<-summary_vars_sel$tissue[summary_vars_sel$ID==var & summary_vars_sel$database=="Intogen"]
  
  count_mat_KRAS[var,tiss_COSMIC]<-count_mat_KRAS[var,tiss_COSMIC]+counts_COSMIC
  count_mat_KRAS[var,tiss_Intogen]<-count_mat_KRAS[var,tiss_Intogen]+counts_Intogen
 
}

write.csv(count_mat_KRAS, file="count_mat_KRAS.csv")


count_mat_sel<-count_mat[count_mat[,"ALL"]>0,]

count_mat_sel<-count_mat_sel[order(count_mat_sel[,"ALL"], decreasing = T),]



count_mat<-matrix(0, ncol=35, nrow=length(unique(summary_vars_sel$ID)))
rownames(count_mat)<-unique(summary_vars_sel$ID)
colnames(count_mat)<-c("ALL", "ALL_MATCHING", tissues)
for(var in unique(summary_vars_sel$ID)){
  counts<-summary_vars_sel$tot_patients[summary_vars_sel$ID==var]
  tiss<-summary_vars_sel$tissue[summary_vars_sel$ID==var]
  count_mat[var,tiss]<-counts
}

count_mat_sel<-count_mat[count_mat[,"ALL"]>0,]

count_mat_sel<-count_mat_sel[order(count_mat_sel[,"ALL"], decreasing = T),]

library(pheatmap)
library(colorBlindness)

pheatmap(log10(count_mat_sel+1)[1:50,], cluster_cols = F, cluster_rows = F)

toplot<-log10(count_mat_sel+1)[1:100,]
anno<-data.frame(Driver=summary_vars_sel$driver[match(rownames(toplot), summary_vars_sel$ID)], 
                 PolyPhen=summary_vars_sel$PolyPhen[match(rownames(toplot), summary_vars_sel$ID)],
                SIFT=summary_vars_sel$SIFT[match(rownames(toplot), summary_vars_sel$ID)])
rownames(anno)<-rownames(toplot)
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)
myBreaks <- c(seq(0,max(unlist(toplot), na.rm=T), length.out=floor(paletteLength)))
length(myBreaks) == length(paletteLength) + 1

anno_colors<-list(SIFT=c(deleterious=Blue2DarkRed18Steps[16], deleterious_low_confidence=Blue2DarkRed18Steps[13], 
                         tolerated=Blue2DarkRed18Steps[3], tolerated_low_confidence=Blue2DarkRed18Steps[6]), 
                  PolyPhen=c(benign=Blue2DarkRed18Steps[3], possibly_damaging=Blue2DarkRed18Steps[13], probably_damaging=Blue2DarkRed18Steps[16]),
                  Driver=c(Driver=Blue2DarkRed18Steps[16], "Non-driver"="darkgrey"))
pdf("figures/Pheatmap_counts_tiss.pdf",15,20)
pheatmap(toplot,   cellwidth=10, cellheight=10, breaks=myBreaks, color = myColor, 
         keep.dendro=T, annotation_row = anno, annotation_colors = anno_colors)
dev.off()
