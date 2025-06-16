library(pheatmap)
library(ggplot2)
library(openxlsx)
library(ReactomePA)
library(clusterProfiler)
library(tidyverse)

#set path
path_data<-"/data"
path_results<-"/results/20220208"
home<-"/Volumes/home/VUS/VUS"

setwd(paste(home, "/", path_results, sep=""))

##load DAM data
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




#################
####random p, rank ratio and median fitness effect distribution
###############
p<-c()
for(ctiss in tissues){
p<-c(p,results[[ctiss]]$pval_rand[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5])
}

df<-data.frame(pvalue=p)
pdf("figures/p_distribution.pdf",6,4)
ggplot(df, aes(x=pvalue))+geom_histogram(bins = 100)+theme_classic()+ geom_vline(xintercept = 0.2, linetype="dashed", color = "red", size=0.5)
dev.off()


rr<-c()
for(ctiss in tissues){
  rr<-c(rr,results[[ctiss]]$rank_ratio)
}

df<-data.frame(RankRatio=rr)
pdf("figures/rr_distribution.pdf",6,4)
ggplot(df, aes(x=RankRatio))+geom_histogram(bins = 100)+theme_classic()+ geom_vline(xintercept = 1.6, linetype="dashed", color = "red", size=0.5)
dev.off()

mfe<-c()
for(ctiss in tissues){
  mfe<-c(mfe,results[[ctiss]]$medFitEff)
}

df<-data.frame(medFitEff=mfe)
pdf("figures/medFitEff_distribution.pdf",6,4)
ggplot(df, aes(x=medFitEff))+geom_histogram(bins = 100)+theme_classic()+ geom_vline(xintercept = -0.5, linetype="dashed", color = "red", size=0.5)
dev.off()

####################
### Number of variants
#######################

num_vars<-c()
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)

  num_vars<-c(num_vars, length(unlist(strsplit(results[[ctiss]]$var[ind], " \\| "))))
}
sum(num_vars)

#############################
### How many hits per tissue?
##############################

num_hits<-c()
for(ctiss in tissues){
  num_hits<-c(num_hits, length(unique(results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2])))
}
names(num_hits)<-tissues

hist(num_hits, breaks=20)
barplot(sort(num_hits, decreasing = T)[1:10], las=2)

paste("the tissue with the highest number of hits is ", names(num_hits)[which.max(num_hits)],
      "with", num_hits[which.max(num_hits)], "hits", sep=" ")
paste("the tissue with the lowest number of hits is ", names(num_hits)[which.min(num_hits)],
      "with", num_hits[which.min(num_hits)], "hits", sep=" ")

##############################
##taking tumor burden into account
##############################
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

mut_burden_anno<-c()
for(ctiss in tissues){
  mut_burden_anno<-c(mut_burden_anno, mean(CMP_annot$mutational_burden[CMP_annot$cancer_type==ctiss], na.rm=T))
}
names(mut_burden_anno)<-tissues

mut_burden_new<-c()
for(ctiss in tissues){
  cl_variants_tmp<-cl_variants[cl_variants$cancer_type==ctiss,]
  num_mut<-c()
  for(model in unique(cl_variants_tmp$model_name)){
    num_mut<-c(num_mut,length(unique(cl_variants_tmp$gene_symbol[cl_variants_tmp$model_name==model])))
  }
  mut_burden_new<-c(mut_burden_new, mean(num_mut, na.rm=T))
}
names(mut_burden_new)<-tissues

plot(mut_burden_anno, mut_burden_new)#i due valori sono correlati

num_hits_norm<-num_hits/mut_burden_new
hist(num_hits_norm, breaks=20)


##########################
#############
### How many drivers per cancer type?
##########
cancer_match_long_CMP<-c()
cancer_match_long_into<-c()
for(i in 1:nrow(cancer_match)){
  cancer_match_long_into<-c(cancer_match_long_into, unlist(strsplit(cancer_match[i,2], " \\| ")))
  cancer_match_long_CMP<-c(cancer_match_long_CMP, rep(cancer_match[i,1],length(unlist(strsplit(cancer_match[i,2], " \\| ")))))
}
cancer_match_long<-cbind(cancer_match_long_CMP, cancer_match_long_into)
num_drivers<-table(unlist(strsplit(  unique(paste(inTOgen_drivers$SYMBOL, inTOgen_drivers$CANCER_TYPE)), " "))[seq(2, (2*nrow(inTOgen_drivers)),2)])
plot(num_hits[cancer_match_long[match(names(num_drivers), cancer_match_long[,2]),1]], num_drivers)
cor.test(num_hits[cancer_match_long[match(names(num_drivers), cancer_match_long[,2]),1]], num_drivers)

plot(num_hits_norm[cancer_match_long[match(names(num_drivers), cancer_match_long[,2]),1]], num_drivers)
cor.test(num_hits_norm[cancer_match_long[match(names(num_drivers), cancer_match_long[,2]),1]], num_drivers)

######################

paste("the tissue with the percentage of hits is ", names(num_hits_norm)[which.max(num_hits_norm)],
      "with", num_hits_norm[which.max(num_hits_norm)], "hits", sep=" ")
paste("the tissue with the lowest percentage hits is ", names(num_hits_norm)[which.min(num_hits_norm)],
      "with", num_hits_norm[which.min(num_hits_norm)], "hits", sep=" ")

#############################
### Hits frequency across tissues 
##############################

hits<-c()
for(ctiss in tissues){
  hits<-c(hits, unique(results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2]))
}

paste("Number of unique gene hits:", length(unique(hits)))

hits_2ormore<-sort(table(hits), decreasing=T)[1:sum(table(hits)>1)]
length(hits_2ormore) #203sono in pi� di un tessuto

paste("Number of gene hits in at least two cancer types:", length(hits_2ormore))

setdiff(names(hits_2ormore), driver_genes)

paste("of which not known to be drivers:", length(setdiff(names(hits_2ormore), driver_genes)))

paste("Number of gene hits known as driver:", length(intersect(hits, driver_genes)))

paste("Number of gene hits not known as driver:", length(setdiff(hits, driver_genes)))

df<-data.frame(counts=table(table(hits)))
df$counts.Freq<-log10(df$counts.Freq+1)
df$counts.Var1<-as.numeric(as.character(df$counts.Var1))

pdf("figures/num_tissues_perhit_log.pdf", 5, 5)
ggplot(df, aes(x=counts.Var1, y=counts.Freq))+geom_bar(stat="identity")+theme_classic()+
  ylab("log10(counts+1)")+xlab("Number of cancer types")+scale_x_continuous(breaks=c(1:15))
dev.off()


df<-data.frame(counts=table(table(hits)[setdiff(hits,driver_genes)]))
df$counts.Freq<-log10(df$counts.Freq+1)
df$counts.Var1<-as.numeric(as.character(df$counts.Var1))

pdf("figures/num_tissues_perhit_nondrivers_log.pdf", 2,5)
ggplot(df, aes(x=counts.Var1, y=counts.Freq))+geom_bar(stat="identity")+theme_classic()+
  ylab("log10(counts+1)")+xlab("Number of cancer types")+scale_x_continuous(breaks=c(1:10))
dev.off()


sort(table(hits)[setdiff(hits,driver_genes)])

###enrichment for drivers
a<-length(intersect(hits, driver_genes))
b<-length(setdiff(hits, driver_genes))
c<-length(setdiff(driver_genes, hits))
d<-length(unique(cl_variants$gene_symbol))

fisher.test(matrix(c(a,b,c,d), ncol=2), alternative="greater")

############
## Cell lines with the DAMs
############
lines_screened<-CMP_annot$model_name[which(CMP_annot$cancer_type %in% tissues)]

lines<-c()
for(ctiss in tissues){
ind<-which(results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)

for(i in 1:length(ind)){
gene_ctiss<-results[[ctiss]]$GENE[ind[i]]
vars_ctiss<-unlist(strsplit(results[[ctiss]]$var[ind[i]], " \\| "))
lines_var<-cl_variants$model_name[cl_variants$gene_symbol==gene_ctiss & cl_variants$protein_mutation %in% vars_ctiss & cl_variants$cancer_type==ctiss]
lines<-c(lines, intersect(lines_var, lines_screened))
}
}

hist(sort(table(lines)/(table(cl_variants$model_name)[table(lines)])), breaks=100)

df<-data.frame(line=lines, cancer_type=CMP_annot$cancer_type[match(lines, CMP_annot$model_name)])
df$cancer_type<-factor(df$cancer_type, levels=c(names(sort(table(df$cancer_type), decreasing=T))))
df$line_by_type<-NA
for(ctiss in tissues){
  line_by_type<-df$line[df$cancer_type==ctiss]
  df$line_by_type[df$cancer_type==ctiss]<-match(line_by_type, unique(line_by_type))
}

num_lines<-c()
unique_lines<-c()
for(ctiss in levels(df$cancer_type)){
  unique_lines<-c(unique_lines, max(df$line_by_type[df$cancer_type==ctiss]))
  num_lines<-c(num_lines, length(which(df$cancer_type==ctiss)))
}

df_text<-data.frame(cancer_type=levels(df$cancer_type), num_lines=num_lines, unique_lines=unique_lines)

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

set.seed(245374)
pdf("figures/cell_lines_DAMs.pdf",6,5)
ggplot(df, aes(x=cancer_type, fill=as.character(line_by_type)))+geom_bar()+theme_classic()+
  scale_fill_manual(values=sample(col_vector, 44))+geom_text(data=df_text, aes(x=cancer_type, y=num_lines+5, label=unique_lines), color="black", fontface="bold",alpha=0.6, size=2, inherit.aes = FALSE)+
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("cell lines")
dev.off()

print(paste("Number of cell lines with DAMs", sum(df_text$unique_lines) ))
sum(df_text$unique_lines)/length(lines_screened)

summary_lines<-data.frame(gene=NA, var=NA, line=NA, cancer_type=NA)
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)

  for(i in 1:length(ind)){
    gene_ctiss<-results[[ctiss]]$GENE[ind[i]]
    vars_ctiss<-unlist(strsplit(results[[ctiss]]$var[ind[i]], " \\| "))

    for(var in vars_ctiss){
    lines_tmp<-cl_variants$model_name[cl_variants$gene_symbol==gene_ctiss & cl_variants$protein_mutation %in% var & cl_variants$cancer_type==ctiss]
    summary_tmp<-data.frame(gene=rep(gene_ctiss, length(lines_tmp)), var=rep(var, length(lines_tmp)), line=lines_tmp, cancer_type=rep(ctiss, length(lines_tmp)))
    summary_lines<-rbind.data.frame(summary_lines, summary_tmp)
    }
  }
}

#####Are there drugs targeting a DAM-bearing-gene, not yet tested on the mutated cell line? 

drugTargetInfo <- read.table('../../Data/raw/drug-target_data_hgvs_clean_Goncalves_et_all.txt',sep='\t',stringsAsFactors = FALSE,header=TRUE)
gdsc1<-read.csv('../../Data/raw/GDSC1_fitted_dose_response_25Feb20.csv',header = TRUE,stringsAsFactors = FALSE)
gdsc2<-read.csv('../../Data/raw/GDSC2_fitted_dose_response_25Feb20.csv',header = TRUE,stringsAsFactors = FALSE)
colnames(gdsc1)[1]<-"DATASET"
colnames(gdsc2)[1]<-"DATASET"

gdscAll<-rbind(gdsc1,gdsc2)

drugs<-drugTargetInfo$Gene.Target %in% hits

summary_lines$tested_drug<-NA
summary_lines$tested_drug_line<-NA
for(i in 1:nrow(summary_lines)){
  if(summary_lines$gene[i] %in% drugTargetInfo$Gene.Target){
    drug<-drugTargetInfo$Drug.ID[which(drugTargetInfo$Gene.Target==summary_lines$gene[i])]
    lines_tested<-gdscAll$CELL_LINE_NAME[which(gdscAll$DRUG_ID %in% drug)]
    if(length(lines_tested)>0){
      summary_lines$tested_drug[i]<-"Yes"
      if(summary_lines$line[i] %in% lines_tested){
        summary_lines$tested_drug_line[i]<-"Yes"
      } else {
        summary_lines$tested_drug_line[i]<-"No"
      }
    }
  }
}

setdiff(summary_lines$gene[summary_lines$tested_drug_line=="No"], driver_genes)

########################################
#### Overlap between hits and driver genes
###################################
library(ggvenn)
pdf("figures/Venn_driver_hits.pdf", 5, 5)
ggvenn(data=list('Driver genes'=driver_genes, 'Hits'=hits),text_size=5, fill_color = c("#F0E442", "#0072B2"), show_percentage = F)
dev.off()
perc_int<-length(intersect(hits, driver_genes))/length(hits)
print(paste("Intersection", perc_int))

###enrichment for Act drivers
act_drivers<-unique(inTOgen_drivers$SYMBOL[inTOgen_drivers$ROLE=="Act"])
a<-length(intersect(hits, act_drivers))
b<-length(setdiff(hits, act_drivers))
c<-length(setdiff(act_drivers, hits))
d<-length(unique(cl_variants$gene_symbol))

fisher.test(matrix(c(a,b,c,d), ncol=2), alternative="greater")

###enrichment for LoF drivers
LoF_drivers<-unique(inTOgen_drivers$SYMBOL[inTOgen_drivers$ROLE=="LoF"])
a<-length(intersect(hits, LoF_drivers))
b<-length(setdiff(hits, LoF_drivers))
c<-length(setdiff(LoF_drivers, hits))
d<-length(unique(cl_variants$gene_symbol))

fisher.test(matrix(c(a,b,c,d), ncol=2), alternative="greater")


pdf("figures/Venn_driver_LoFAct_hits.pdf", 5, 5)
ggvenn(data=list( 'Act'=act_drivers, 'LoF'=LoF_drivers, 'Hits'=hits),text_size=5,  show_percentage = F)
dev.off()

######################
#### Other lists of drivers
######################
benchmark<-read.xlsx("../../data/benchmark-datasets.xlsx",1)

for(col in 1:8){
  dg1<-na.omit(benchmark[-1,col])

  a<-length(intersect(hits, dg1))
  b<-length(setdiff(hits, dg1))
  c<-length(setdiff(dg1, hits))
  d<-length(unique(cl_variants$gene_symbol))

  fisher.test(matrix(c(a,b,c,d), ncol=2), alternative="greater")

  pdf(paste("figures/Venn_driver_", col, ".pdf", sep=""), 5, 5)
  print(ggvenn(data=list('Driver genes'=dg1, 'Hits'=hits),text_size=5, fill_color = c("#F0E442", "#0072B2"), show_percentage = F))
  dev.off()

}


###enrichment for drivers

driver_union<-unique(c(driver_genes,c(na.omit(unlist(benchmark[-1,])))))
write.xlsx(data.frame(driver_union), "driver_union.xlsx")

a<-length(intersect(hits, driver_union))
b<-length(setdiff(hits, driver_union))
c<-length(setdiff(driver_union, hits))
d<-length(unique(cl_variants$gene_symbol))

fisher.test(matrix(c(a,b,c,d), ncol=2), alternative="greater")

pdf(paste("figures/Venn_driver_union.pdf", sep=""), 5, 5)
ggvenn(data=list('Driver genes'=driver_union, 'Hits'=hits),text_size=5, fill_color = c("#F0E442", "#0072B2"), show_percentage = F)
dev.off()

#######################
####How many hits in common between different tissues?  heatmap with the Jaccard index
#####################

JI <- matrix(NA, ncol=length(tissues), nrow=length(tissues))
rownames(JI) <- tissues
colnames(JI) <- tissues

num_int <- matrix(NA, ncol=length(tissues), nrow=length(tissues))
rownames(num_int) <- tissues
colnames(num_int) <- tissues


for(first_tiss in tissues){

  RESTOT<-results[[first_tiss]]
  ft_genes <- RESTOT$GENE[RESTOT$rank_ratio<1.6 & RESTOT$medFitEff< -.5]

  for(second_tiss in tissues){

    RESTOT<-results[[second_tiss]]
    st_genes <- RESTOT$GENE[RESTOT$rank_ratio<1.6 & RESTOT$medFitEff< -.5]

    JI[first_tiss, second_tiss]<-length(intersect(ft_genes, st_genes))/length(union(ft_genes, st_genes))
    num_int[first_tiss, second_tiss]<-length(intersect(ft_genes, st_genes))

  }
}

diag(JI)<-NA
diag(num_int)<-NA
pheatmap(JI)
pheatmap(num_int)

###################################
####How many drivers do we detect correctly? 
###################################

cancer_symbols<-c("AML", "DLBCL", "ALL",
                  "CH", "BLCA", "BRCA",
                  "CESC","MDPS", "COREAD",
                  "UCEC", "ESCA", "ESCA",
                  "EWS", "ST", "GBM",
                  "HGG | LGG", "HNSC", "HC",
                  "RCCC | RCH | RPC", "CM", "MESO",
                  "NB", "NSCLC | LUAD", "HNSC",
                  "OS", "OV", "PAAD",
                  "MM", "RHBDS", "SCLC",
                  "LUSC", "NHLY", "THCA")

cancer_match<-data.frame(tissues=(tissues), symbols=(cancer_symbols))
save(cancer_match, file="cancer_match.RData")

cancer_match<-read.xlsx("cancer_match_ext.xlsx")

known_drivers<-matrix(nrow=length(tissues), ncol=5)
rownames(known_drivers)<-tissues
colnames(known_drivers)<-c("Number known", "Of these found", "Total found", "Novel", "Not found")

for(ctiss in tissues){
  known<-unique(inTOgen_drivers$SYMBOL[inTOgen_drivers$CANCER_TYPE %in% (cancer_match[which(cancer_match[,1]==ctiss),2])])
  total_found<-intersect(results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2], driver_genes)

  known_drivers[ctiss,"Number known"]<-length(known)
  known_drivers[ctiss,"Total found"]<-length(total_found)
  known_drivers[ctiss,"Of these found"]<-length(intersect(known, total_found))
  known_drivers[ctiss,"Novel"]<-length(setdiff(total_found, known))
  known_drivers[ctiss,"Not found"]<-length(setdiff(known, total_found))
}


df<-data.frame(counts=c(known_drivers[,"Of these found"], known_drivers[,"Not found"]),
               type=rep(c("Found", "Not found"), each=length(tissues)), cancer_type=c(tissues, tissues))
ggplot(df, aes(x=cancer_type, y=counts, fill=type))+geom_col(position = position_dodge(0.8))+theme(axis.text.x = element_text(angle=90))

##FOR MOST CANCER DRIVERS, MUTATIONS DO NOT ALTER GENE'S ESSENTIALITY

df<-data.frame(counts=c(known_drivers[,"Novel"], known_drivers[,"Total found"]),
               type=rep(c("Novel", "Total found"), each=length(tissues)), cancer_type=c(tissues, tissues))
ggplot(df, aes(x=cancer_type, y=counts, fill=type))+geom_col(position = position_dodge(0.8))+theme(axis.text.x = element_text(angle=90))


##################################
##Are cancer drivers known in one tissue also hits in other tissues? 
#Heatmap with tissues in columns, drivers in rows, red if hit 
#Or histogram of the number of tissues where a gene is hit, blue if we already knew it was a driver, red if we didn't know 
########################################

summary_drivers<-matrix(nrow=length(driver_genes), ncol=length(tissues))
rownames(summary_drivers)<-driver_genes
colnames(summary_drivers)<-tissues

for(dg in driver_genes){
  for(ctiss in tissues){
    if(dg %in% results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2]){
      summary_drivers[dg,ctiss]<-"Novel"
      known<-cancer_match[which(cancer_match[,1]==ctiss),2]

    if(length(intersect(known, inTOgen_drivers$CANCER_TYPE[inTOgen_drivers$SYMBOL==dg]))>0){
      summary_drivers[dg,ctiss]<-"Known Found"
    }
    } else if(length(intersect(known, inTOgen_drivers$CANCER_TYPE[inTOgen_drivers$SYMBOL==dg]))>0){
      summary_drivers[dg,ctiss]<-"Known Not Found"
    }
  }
}

summary_drivers_bin<-summary_drivers
summary_drivers_bin[summary_drivers_bin=="Known Found"]<-1
summary_drivers_bin[summary_drivers_bin=="Novel"]<-0
class(summary_drivers_bin)<-"numeric"
pheatmap(summary_drivers_bin, cluster_rows = F, cluster_cols = F)

hist(rowSums(summary_drivers_bin==1, na.rm=T))#driver known in how many tissues?
hist(rowSums(summary_drivers_bin==0, na.rm=T))#driver novel in how many tissues?
hist(colSums(summary_drivers_bin==0, na.rm=T), breaks=20)#tissues with how many novel? #top breast carcinoma with 7

#for each cancer type, plot how many driver genes were known (and are hits) and how many novel ones we indicate
df<-data.frame(counts=c(colSums(summary_drivers=="Known Found", na.rm=T), colSums(summary_drivers=="Known Not Found", na.rm=T), colSums(summary_drivers=="Novel", na.rm=T)),
               type=rep(c("Known Found", "Known Not Found","Novel"), each=length(tissues)), cancer_type=c(tissues, tissues, tissues))
ggplot(df, aes(x=cancer_type, y=counts, fill=type))+geom_col(position = position_dodge(0.8))+theme(axis.text.x = element_text(angle=90))



#for each diver genes, plotted how many tumor types were known and how many novel ones we find
selhits<-names(which(rowSums(summary_drivers==c("Known Found"), na.rm=T)>0|rowSums(summary_drivers==c("Novel"), na.rm=T)>0))
df<-data.frame(counts=c(rowSums(summary_drivers=="Known Found", na.rm=T), rowSums(summary_drivers=="Known Not Found", na.rm=T), rowSums(summary_drivers=="Novel", na.rm=T)),
               type=rep(c("Found in a Known Cancer Type", "Not Found in a Known Cancer Type","Found in a Novel Cancer Type"), each=nrow(summary_drivers)), gene=c(rownames(summary_drivers), rownames(summary_drivers), rownames(summary_drivers)))

df$gene<-factor(df$gene, levels=c(df$gene[df$type=="Found in a Known Cancer Type"][order(df$counts[df$type=="Found in a Known Cancer Type"]+df$counts[df$type=="Found in a Novel Cancer Type"])]))

df$counts[df$type=="Not Found in a Known Cancer Type"]<- -df$counts[df$type=="Not Found in a Known Cancer Type"]
df$type<-factor(df$type, levels=c("Found in a Known Cancer Type", "Found in a Novel Cancer Type", "Not Found in a Known Cancer Type"))

png("figures/Drivers_barplot_mirror.png", res=300, 5000, 4000)
ggplot(subset(df, (gene %in% selhits)),aes(x=gene, y=counts, fill=type))+geom_bar(stat="identity")+theme_classic()+ylab("Number of cancer types")+xlab("")+theme(axis.text.x = element_text(angle=90, size=10))+
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#D55E00"))+ guides(fill=guide_legend(title="Cancer Type"))
dev.off()

pdf("figures/Drivers_barplot.pdf",20, 5)
ggplot(subset(df, (gene %in% selhits) & (type !="Not Found in a Known Cancer Type")),aes(x=gene, y=counts, fill=type))+geom_bar(stat="identity")+theme_classic()+ylab("Number of cancer types")+xlab("")+theme(axis.text.x = element_text(angle=90, size=10))+
  scale_fill_manual(values = c("#0072B2", "#56B4E9"), labels= c("Known", "Novel"))+ guides(fill=guide_legend(title="Cancer Type"))
dev.off()

df<-data.frame(counts=c(rowSums(summary_drivers_bin==1, na.rm=T), rowSums(summary_drivers_bin==0, na.rm=T)),
               type=rep(c("Known", "Novel"), each=nrow(summary_drivers_bin)), cancer_type=c(rownames(summary_drivers_bin), rownames(summary_drivers_bin)))

ggplot(df, aes(x=cancer_type, y=counts, fill=type))+geom_col(position = position_dodge(0.8))+theme(axis.text.x = element_text(angle=90))

sort(rowSums(summary_drivers==c("Novel"), na.rm=T), decreasing=T)[1:10]
#NRAS is the top one with 6 tissues where it wasn't annotated
table(summary_drivers["NRAS",])
which(summary_drivers["NRAS",]=="Novel")
which(summary_drivers["PIK3CA",]=="Novel")
#FLCN, HRAS, HSP90AB1, SMARCA4, SOS1, TSC2 with 3


#for each diver genes, plotted how many tumor types were known and how many novel ones we find
selhits<-names(which(rowSums(summary_drivers==c("Known Found")|summary_drivers==c( "Novel"), na.rm=T)>1))
df<-data.frame(counts=c(rowSums(summary_drivers=="Known Found", na.rm=T), rowSums(summary_drivers=="Known Not Found", na.rm=T), rowSums(summary_drivers=="Novel", na.rm=T)),
               type=rep(c("Found in a Known Cancer Type", "Not Found in a Known Cancer Type","Found in a Novel Cancer Type"), each=nrow(summary_drivers)), gene=c(rownames(summary_drivers), rownames(summary_drivers), rownames(summary_drivers)))

df$gene<-factor(df$gene, levels=c(df$gene[df$type=="Found in a Known Cancer Type"][order(df$counts[df$type=="Found in a Known Cancer Type"]+df$counts[df$type=="Found in a Novel Cancer Type"])]))

df$counts[df$type=="Not Found in a Known Cancer Type"]<- -df$counts[df$type=="Not Found in a Known Cancer Type"]
df$type<-factor(df$type, levels=c("Found in a Known Cancer Type", "Found in a Novel Cancer Type", "Not Found in a Known Cancer Type"))


pdf("figures/Drivers_barplot_thr2.pdf",10, 5)
ggplot(subset(df, (gene %in% selhits) & (type !="Not Found in a Known Cancer Type")),aes(x=gene, y=counts, fill=type))+geom_bar(stat="identity")+theme_classic()+ylab("Number of cancer types")+xlab("")+theme(axis.text.x = element_text(angle=90, size=10))+
  scale_fill_manual(values = c("#0072B2", "#56B4E9"), labels= c("Known", "Novel"))+ guides(fill=guide_legend(title="Cancer Type"))
dev.off()

###############################
### Genes not known to be drivers: in how many cancer types 
################################
hits_nodriver<-setdiff(unique(hits), driver_genes)
save(hits_nodriver, file="hits_nodriver.RData")
hits<-unique(hits)
save(hits, file="hits.RData")

summary_nodrivers<-matrix(0, nrow=length(hits_nodriver), ncol=length(tissues))
rownames(summary_nodrivers)<-hits_nodriver
colnames(summary_nodrivers)<-tissues

for(dg in hits_nodriver){
  for(ctiss in tissues){
    if(dg %in% results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2]){
      summary_nodrivers[dg,ctiss]<-1
    }

  }
}

pheatmap(summary_nodrivers, cluster_rows = F, cluster_cols = F)

hist(rowSums(summary_nodrivers==1, na.rm=T))#any novel in how many tissues? 
hist(colSums(summary_nodrivers==1, na.rm=T), breaks=20)#in one tissue how many novel drivers?

df<-data.frame(counts=c(colSums(summary_drivers=="Known Found", na.rm=T), colSums(summary_drivers=="Novel", na.rm=T), colSums(summary_nodrivers==1)),
               type=rep(c("Known Driver, right tissue", "Known Driver, wrong tissue","Not Known as Driver"),each=33), tissue=tissues)
df$counts<-df$counts/rep(mut_burden_new, 3)
df$tissue<-factor(df$tissue, levels=c(df$tissue[order(num_hits_norm, decreasing=T)]))
df_sel<-df[df$tissue %in% df$tissue[order(num_hits_norm, decreasing=T)][1:10],]

pdf("figures/Hits_tissue_stack.pdf", 5, 5)
ggplot(df_sel, aes(x=tissue, y=counts, fill=type))+geom_bar(stat="identity")+theme_classic()+ylab("% of hits")+xlab("")+theme(axis.text.x = element_text(angle=90, size=10))+
scale_fill_manual(values = c("#0072B2", "#56B4E9", "#E69F00"))+ guides(fill=guide_legend(title="Hit"))
dev.off()


df<-data.frame(counts=c(colSums(summary_drivers=="Known Found", na.rm=T), colSums(summary_drivers=="Novel", na.rm=T), colSums(summary_nodrivers==1)),
               type=rep(c("Known Driver", "Known Driver","Not Known as Driver"),each=33), tissue=tissues)
df$counts<-df$counts/rep(mut_burden_new, 3)
df$tissue<-factor(df$tissue, levels=c(df$tissue[order(num_hits_norm, decreasing=T)]))
df_sel<-df[df$tissue %in% df$tissue[order(num_hits_norm, decreasing=T)][1:10],]

pdf("figures/Hits_tissue_stack_simpl.pdf", 5, 5)
ggplot(df_sel, aes(x=tissue, y=counts, fill=type))+geom_bar(stat="identity")+theme_classic()+ylab("% of hits")+xlab("")+theme(axis.text.x = element_text(angle=90, size=10))+
  scale_fill_manual(values = c("#0072B2",  "#E69F00"))+ guides(fill=guide_legend(title="Hit"))
dev.off()


selhits<-names(which(rowSums(summary_nodrivers, na.rm=T)>1))
df<-data.frame(counts=c(rowSums(summary_nodrivers, na.rm=T)), gene=c(rownames(summary_nodrivers)))

df$gene<-factor(df$gene, levels=c(df$gene[order(df$counts)]))

pdf("figures/NODrivers_barplot.pdf",20, 5)
ggplot(subset(df, (gene %in% selhits)),aes(x=gene, y=counts))+geom_bar(stat="identity")+theme_classic()+ylab("Number of cancer types")+xlab("")+theme(axis.text.x = element_text(angle=90, size=10))
dev.off()

#genes with mutations changig function in many tumor types
barplot(sort(rowSums(summary_nodrivers==1, na.rm=T), decreasing=T)[1:10], las=2)

for(i in 1:10){
gene<-names(sort(rowSums(summary_nodrivers==1, na.rm=T), decreasing=T))[i]
variants<-c()
for(ctiss in tissues){
  if(gene %in% unique(results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2])){

    variants<-c(variants, strsplit(results[[ctiss]][results[[ctiss]]$GENE==gene, "var"], " \\| "))
  }
}
#all different variants
variants<-unlist(variants)
print(gene)
print(table(variants))
}

##mutations' position 
pos<-as.numeric(sub(".*?([0-9]+).*", "\\1",cl_variants$protein_mutation[cl_variants$gene_symbol=="DIDO1"]))

###variants treated separately, which are hits in two tissues? 
vars_tot<-c()
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$GENE %in% hits_nodriver & results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)
vars<-strsplit(results[[ctiss]]$var[ind], " \\| ")
genes<-rep(results[[ctiss]]$GENE[ind], lapply(vars, length))
vars_tot<-c(vars_tot, paste(genes, unlist(vars), sep="-"))
}

sort(table(vars_tot), decreasing=T)[1:10]

###how many variants on the same ammminoacid are hits in two tissues? 
vars_tot<-c()
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$GENE %in% hits_nodriver & results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)
  vars<-strsplit(results[[ctiss]]$var[ind], " \\| ")

  genes<-rep(results[[ctiss]]$GENE[ind], lapply(vars, length))

  vars<-gsub("\\*.*","", unlist(vars))
  vars<-as.numeric(gsub("\\D", "", unlist(vars)))
  vars_tot<-c(vars_tot, paste(genes, unlist(vars), sep="-"))
}

sort(table(vars_tot), decreasing=T)[1:10]

vars_tot<-c()
vars_tot_df<-data.frame(gene_id=NA, gene_symbol_2019=NA, gene_symbol=NA, var=NA)
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$GENE %in% hits & results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)
  vars<-strsplit(results[[ctiss]]$var[ind], " \\| ")
  genes<-rep(results[[ctiss]]$GENE[ind], lapply(vars, length))
  vars_tot<-c(vars_tot, paste(genes, gsub("p.", "", unlist(vars)), sep="-"))
  gene_id<-cl_variants$gene_id[match(genes, cl_variants$gene_symbol_2019)]
  gene_new<-cl_variants$gene_symbol[match(genes, cl_variants$gene_symbol_2019)]
  vars_tot_df<-rbind.data.frame(vars_tot_df, data.frame(gene_id=gene_id, gene_symbol_2019=genes, gene_symbol=gene_new, var=gsub("p.", "", unlist(vars))))
}
vars_tot_df<-vars_tot_df[-1,]

write.csv(vars_tot, "allvars.csv")
write.csv(vars_tot_df, "allvars_id.csv")

print(paste("Total number of variants:", length(unique(vars_tot))))

##in which tissues are the genes mutated? how?
for(ctiss in tissues){
  ind<-which(results[[ctiss]]$GENE %in% "ANKRD49" & results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand<0.2)
  vars<-strsplit(results[[ctiss]]$var[ind], " \\| ")

  if(length(vars)>0){
  print(ctiss)
  print(vars)
  }
}


############
###EGO nrichment test
############
#testare l'arricchimento degli hits in un tessuto per GO
#come background, prendere tutti i geni nella lista di HUGO symbols, meno quelli esclusi perch� troppo mutati o perch� in core fitness
#quindi il background � tissue-specific
#meglio scegliere come bckground solo le varianti testate?
library(ReactomePA)
library(clusterProfiler)

toexclude<-unique(c(ADaM, Perc_AUC))
background_all<- setdiff(setdiff(gene_annot$hgnc_symbol, toexclude),driver_genes)

ego<-enrichPathway(gene=bitr(hits_nodriver, "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db")[,2],
                   pvalueCutoff = 0.05, readable=TRUE, universe=bitr(background_all, "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db")[,2])

write.xlsx(summary(ego), paste("Reactome_nodriver_all.xlsx"))

library(enrichplot)
require(DOSE)
pdf("figures/Reactome_nodriver.pdf", 7 ,6)
dotplot(ego, showCategory=10)
dev.off()

background_all<- setdiff(gene_annot$hgnc_symbol, toexclude)

ego<-enrichPathway(gene=bitr(hits, "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db")[,2],
                   pvalueCutoff = 0.05, readable=TRUE, universe=bitr(background_all, "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db")[,2])

write.xlsx(summary(ego), paste("Reactome_all.xlsx"))

pdf("figures/Reactome_all.pdf", 7 ,6)
dotplot(ego, showCategory=10)
dev.off()



REACT_allcat<-c()
for (ctiss in tissues){
  ts_depFC<-scaled_depFC[,CMP_annot$model_name[CMP_annot$cancer_type==ctiss]]
  ts_bdep<-bdep[,CMP_annot$model_name[CMP_annot$cancer_type==ctiss]]

  ts_cl_variants<-cl_variants[which(is.element(cl_variants$model_name,colnames(ts_depFC))),]
  ts_cl_variants<-ts_cl_variants[which(is.element(ts_cl_variants$gene_symbol,rownames(ts_depFC))),]
  genesToTest<-unique(ts_cl_variants$gene_symbol)
  vs_spec_cardinality<-unlist(lapply(lapply(1:length(genesToTest),
                                            function(x){
                                              dd<-variantSpectrum(cl_var = ts_cl_variants,gene = genesToTest[x])
                                              dd<-dd$protein_mutation}),'length'))

  names(vs_spec_cardinality)<-genesToTest

  toexclude<-names(which(vs_spec_cardinality>=10))
  ##### questo esclude gene pan-cancer cf and common-essential genes
  toexclude<-unique(c(toexclude, ADaM, Perc_AUC))

  background_all<- setdiff(gene_annot$hgnc_symbol, toexclude)

  RESTOT<-results[[ctiss]]


  ego<-enrichPathway(gene=bitr(RESTOT$GENE[RESTOT$rank_ratio<1.6 & RESTOT$medFitEff< -.5 & RESTOT$pval_rand<0.2], "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db")[,2],
                     pvalueCutoff = 0.05, readable=TRUE, universe=bitr(background_all, "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db")[,2])

  write.csv(summary(ego), paste("Reactome_all",ctiss,".csv"))

  REACT_allcat<-c(REACT_allcat, data.frame(ego)[,2])
}

sort(table(REACT_allcat), decreasing = T)[1:50]

paths_totest<-names(table(REACT_allcat))[which(table(REACT_allcat)>6)]
#tolto "MAPK1/MAPK3 signaling" perch� d� errore
paths_totest<-paths_totest[-which(paths_totest=="MAPK1/MAPK3 signaling")]
paths_totest<-paths_totest[-which(paths_totest=="Signaling by the B Cell Receptor (BCR)")]

REACT_allcat_nd<-c()
for (ctiss in tissues){
  ts_depFC<-scaled_depFC[,CMP_annot$model_name[CMP_annot$cancer_type==ctiss]]
  ts_bdep<-bdep[,CMP_annot$model_name[CMP_annot$cancer_type==ctiss]]

  ts_cl_variants<-cl_variants[which(is.element(cl_variants$model_name,colnames(ts_depFC))),]
  ts_cl_variants<-ts_cl_variants[which(is.element(ts_cl_variants$gene_symbol,rownames(ts_depFC))),]
  genesToTest<-unique(ts_cl_variants$gene_symbol)
  vs_spec_cardinality<-unlist(lapply(lapply(1:length(genesToTest),
                                            function(x){
                                              dd<-variantSpectrum(cl_var = ts_cl_variants,gene = genesToTest[x])
                                              dd<-dd$protein_mutation}),'length'))

  names(vs_spec_cardinality)<-genesToTest

  toexclude<-names(which(vs_spec_cardinality>=10))
  ##### questo esclude gene pan-cancer cf and common-essential genes
  toexclude<-unique(c(toexclude, ADaM, Perc_AUC, driver_genes))


  background_all<- setdiff(gene_annot$hgnc_symbol, toexclude)

  RESTOT<-results[[ctiss]]


  ego<-enrichPathway(gene=bitr(setdiff(RESTOT$GENE[RESTOT$rank_ratio<1.6 & RESTOT$medFitEff< -.5 & RESTOT$pval_rand<0.2], driver_genes), "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db")[,2],
                     pvalueCutoff = 0.05, readable=TRUE, universe=bitr(background_all, "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db")[,2])

  write.csv(summary(ego), paste("Reactome_all_nd",ctiss,".csv"))

  REACT_allcat_nd<-c(REACT_allcat_nd, data.frame(ego)[,2])
}

sort(table(REACT_allcat_nd), decreasing = T)[1:10]

###Reactome pathway pi� frequenti, summary
#Arriccimenti mettendo tutti i geni dei vari tumori insieme
#Arricchimento ClinVar (tipo di effect annotato?)

path_newdrivers<-matrix(nrow=length(paths_totest), ncol=1)
rownames(path_newdrivers)<-paths_totest


path_newdrivers_uniq<-matrix(nrow=length(paths_totest), ncol=1)
rownames(path_newdrivers_uniq)<-paths_totest

for(path_name in paths_totest){
  print(path_name)
  path<-viewPathway(path_name)
  genes_path<-matrix(0, ncol=length(tissues), nrow=length(path$data$name), dimnames = list(c(path$data$name), c(tissues)))

  for(ctiss in tissues){
    RESTOT<-results[[ctiss]]
    genes<-RESTOT$GENE[RESTOT$rank_ratio<1.6 & RESTOT$medFitEff< -.5]
    genes_path[intersect(genes, path$data$name), ctiss]<-1
  }


  path_nodriver<-setdiff(rownames(genes_path), driver_genes)
  path_newdrivers[path_name,1]<-paste(path_nodriver[rowSums(genes_path[path_nodriver,])>0], collapse="|")

  genes_path<-t(genes_path)
  selgenes<-path_nodriver[colSums(genes_path[,path_nodriver])>0]
  mut_binary[selgenes,]
  drivers_path<-intersect(colnames(genes_path), driver_genes)
  selmatrix<-genes_path[rowSums(genes_path[,drivers_path])==0,]
  path_newdrivers_uniq[path_name,1]<-paste(colnames(selmatrix)[which(selmatrix==1, arr.ind=T)[,2]], collapse="|")

}


path_name<-"Oncogenic MAPK signaling"
pheatmap(t(genes_path))

vars<-c()
for(ctiss in tissues){
  if("ITPR3" %in% unique(results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5])){

    vars<-c(vars,(results[[ctiss]][results[[ctiss]]$GENE=="ITPR3", "var"]))
  }
}

#SRC p.R163W 1 melanoma
#ITPR3 p.Q2135H 1 cervix TCGA-EK-A2RJ-01

ind<-grep("^ITPR3", COSMIC[,"Gene name"])
COSMIC_sel<-COSMIC[ind,]
ind_AA<-grep(gsub(" ", "", paste(vars, collapse="|")), COSMIC_sel[, "Mutation AA"])
COSMIC_sel<-COSMIC_sel[ind_AA, ]
COSMIC_sel<-COSMIC_sel[COSMIC_sel[,"Sample Type"]!="cell-line",]
