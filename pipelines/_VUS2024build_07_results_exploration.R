library(pheatmap)
library(ggplot2)
library(openxlsx)
library(ReactomePA)
library(clusterProfiler)
library(tidyverse)

#set path
path_data<-"data"
pathdata<-"data"
path_results<-"/results/20250221"
home<-"E:/VUS_2024build"

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

####################################################################################################
#### load results with empirical pvalues
####################################################################

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
inTOgen_drivers<-read.table(paste(home, '/data/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
### inTOgene drivers downloaded from https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip on 20241002
driver_genes<-unique(inTOgen_drivers$SYMBOL)

setwd(home)

#################
####random p, rank ratio and median fitness effect distribution
###############
p<-c()
for(ctiss in tissues){
p<-c(p,results[[ctiss]]$pval_rand[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5])
}

df<-data.frame(pvalue=p)
pdf(paste(home, path_results, "/exploration/figures/p_distribution.pdf", sep=""),6,4)
ggplot(df, aes(x=pvalue))+geom_histogram(bins = 100)+theme_classic()+ geom_vline(xintercept = 0.2, linetype="dashed", color = "red", size=0.5)
dev.off()


rr<-c()
for(ctiss in tissues){
  rr<-c(rr,results[[ctiss]]$rank_ratio)
}

df<-data.frame(RankRatio=rr)
pdf(paste(home, path_results, "/exploration/figures/rr_distribution.pdf", sep=""),6,4)
ggplot(df, aes(x=RankRatio))+geom_histogram(bins = 100)+theme_classic()+ geom_vline(xintercept = 1.6, linetype="dashed", color = "red", size=0.5)
dev.off()

mfe<-c()
for(ctiss in tissues){
  mfe<-c(mfe,results[[ctiss]]$medFitEff)
}

df<-data.frame(medFitEff=mfe)
pdf(paste(home, path_results, "/exploration/figures/medFitEff_distribution.pdf", sep=""),6,4)
ggplot(df, aes(x=medFitEff))+geom_histogram(bins = 100)+theme_classic()+ geom_vline(xintercept = -0.5, linetype="dashed", color = "red", size=0.5)
dev.off()


##############################
##taking tumor burden into account
##############################

#compute the tumor burden from the number of mutations in cell lines, compare it with tumor burden annotated in CMP
mut_burden_anno<-c()
for(ctiss in tissues){
  mut_burden_anno<-c(mut_burden_anno, mean(CMP_annot$mutational_burden[CMP_annot$cancer_type==ctiss], na.rm=T))
}
names(mut_burden_anno)<-tissues

mut_burden_new<-c()
for(ctiss in tissues){
  cl_variants_tmp<-cl_variants[cl_variants$model_id %in% CMP_annot$model_id[CMP_annot$cancer_type==ctiss],]
  num_mut<-c()
  for(model in unique(cl_variants_tmp$model_id)){
    num_mut<-c(num_mut,length(unique(cl_variants_tmp$gene_symbol_2023[cl_variants_tmp$model_id==model])))
  }
  mut_burden_new<-c(mut_burden_new, mean(num_mut, na.rm=T))
}
names(mut_burden_new)<-tissues

#plot(mut_burden_anno, mut_burden_new)#the two values are correlated

num_hits<-c()
for(ctiss in tissues){
  num_hits<-c(num_hits, length(unique(results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand < 0.2])))
}
names(num_hits)<-tissues

num_hits_norm<-num_hits/mut_burden_new


##########################
#############
### How many drivers per cancer type?
##########
cancer_match<-read.csv(paste(pathdata,'/raw/intOGen ctype mapping_AS.csv',sep=''),header = TRUE,row.names = 1, sep=";")

cancer_match_long_CMP<-c()
cancer_match_long_into<-c()
for(i in 1:nrow(cancer_match)){
  cancer_match_long_into<-c(cancer_match_long_into, unlist(strsplit(cancer_match[i,1], " \\| ")))
  cancer_match_long_CMP<-c(cancer_match_long_CMP, rep(rownames(cancer_match)[i],length(unlist(strsplit(cancer_match[i,1], " \\| ")))))
}
cancer_match_long<-cbind(cancer_match_long_CMP, cancer_match_long_into)
num_drivers<-table(unlist(strsplit(  unique(paste(inTOgen_drivers$SYMBOL, inTOgen_drivers$CANCER_TYPE)), " "))[seq(2, (2*nrow(inTOgen_drivers)),2)])

######################

paste("the tissue with the highest percentage of DAM-bearing genes (considering mutational burden) is ", names(num_hits_norm)[which.max(num_hits_norm)],
      "with", num_hits_norm[which.max(num_hits_norm)], "DAM-bearing genes", sep=" ")
paste("the tissue with the lowest percentage DAM-bearing genes (considering mutational burden) is ", names(num_hits_norm)[which.min(num_hits_norm)],
      "with", num_hits_norm[which.min(num_hits_norm)], "DAM-bearing genes", sep=" ")

#############################
### DAM-bearing genes frequency across tissues
##############################

hits<-c()
for(ctiss in tissues){
  hits<-c(hits, unique(results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand < 0.2]))
}

paste("Number of unique DAM-bearing genes:", length(unique(hits)))

hits_2ormore<-sort(table(hits), decreasing=T)[1:sum(table(hits)>1)]

paste("Number of DAM-bearing genes in at least two cancer types:", length(hits_2ormore))

paste("of which not known to be drivers:", length(setdiff(names(hits_2ormore), driver_genes)))

paste("Number of DAM-bearing genes known as driver:", length(intersect(hits, driver_genes)))

paste("Number of DAM-bearing genes not known as driver:", length(setdiff(hits, driver_genes)))


df<-data.frame(counts=table(table(hits)))
df$counts.Freq<-log10(df$counts.Freq+1)
df$counts.Var1<-as.numeric(as.character(df$counts.Var1))

pdf(paste(home, path_results, "/exploration/figures/num_tissues_perhit_log.pdf", sep=""), 5, 5)
ggplot(df, aes(x=counts.Var1, y=counts.Freq))+geom_bar(stat="identity")+theme_classic()+
  ylab("log10(counts+1)")+xlab("Number of cancer types")+scale_x_continuous(breaks=c(1:15))
dev.off()


df<-data.frame(counts=table(table(hits)[setdiff(hits,driver_genes)]))
df$counts.Freq<-log10(df$counts.Freq+1)
df$counts.Var1<-as.numeric(as.character(df$counts.Var1))

pdf(paste(home, path_results, "/exploration/figures/num_tissues_perhit_nondrivers_log.pdf", sep=""), 2,5)
ggplot(df, aes(x=counts.Var1, y=counts.Freq))+geom_bar(stat="identity")+theme_classic()+
  ylab("log10(counts+1)")+xlab("Number of cancer types")+scale_x_continuous(breaks=c(1:10))
dev.off()


############
## Cell lines with the DAMs
############
lines_screened<-CMP_annot$model_id[which(CMP_annot$cancer_type %in% tissues)]

lines<-c()
for(ctiss in tissues){
ind<-which(results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand < 0.2)

for(i in 1:length(ind)){
print(i)
  gene_ctiss<-results[[ctiss]]$GENE[ind[i]]
vars_ctiss<-unlist(strsplit(results[[ctiss]]$var[ind[i]], " \\| "))
lines_var<-cl_variants$model_id[cl_variants$gene_symbol_2023==gene_ctiss & cl_variants$protein_mutation %in% vars_ctiss & CMP_annot$cancer_type[match(cl_variants$model_id, CMP_annot$model_id)]==ctiss]
lines<-c(lines, intersect(lines_var, lines_screened))
}
}

df<-data.frame(line=lines, cancer_type=CMP_annot$cancer_type[match(lines, CMP_annot$model_id)])
df$cancer_type<-factor(df$cancer_type, levels=c(names(sort(table(df$cancer_type), decreasing=T))))
df$line_by_type<-NA
for(ctiss in tissues){
  print(ctiss)
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
pdf(paste(home, path_results, "/exploration/figures/cell_lines_DAMs.pdf", sep=""),6,5)
ggplot(df, aes(x=cancer_type, fill=as.character(line_by_type)))+geom_bar()+theme_classic()+
  scale_fill_manual(values=sample(col_vector, 78, replace=T))+geom_text(data=df_text, aes(x=cancer_type, y=num_lines+5, label=unique_lines), color="black", fontface="bold",alpha=0.6, size=2, inherit.aes = FALSE)+
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("cell lines")
dev.off()

print(paste("Number of cell lines with DAMs", sum(df_text$unique_lines) ))
print(paste("Fraction of cell lines with DAMs", sum(df_text$unique_lines)/length(lines_screened) ))


########################################
#### Overlap between DAM-bearing genes and driver genes
###################################
library(ggvenn)
pdf(paste(home, path_results, "/exploration/figures/Venn_driver_hits.pdf", sep=""), 5, 5)
ggvenn(data=list('Driver genes'=driver_genes, 'Hits'=hits),text_size=5, fill_color = c("#F0E442", "#0072B2"), show_percentage = F)
dev.off()
perc_int<-length(intersect(hits, driver_genes))/length(hits)
print(paste("Intersection", perc_int))


######################
#### Other lists of drivers
######################
benchmark<-read.xlsx(paste(path_data,"/raw/benchmark-datasets.xlsx", sep=""),1)

for(col in 1:8){
  dg1<-na.omit(benchmark[-1,col])

  a<-length(intersect(hits, dg1))
  b<-length(setdiff(hits, dg1))
  c<-length(setdiff(dg1, hits))
  d<-length(unique(cl_variants$gene_symbol))

  fisher.test(matrix(c(a,b,c,d), ncol=2), alternative="greater")

  pdf(paste(home, path_results, "/exploration/figures/Venn_driver_", col, ".pdf", sep=""), 5, 5)
  print(ggvenn(data=list('Driver genes'=dg1, 'Hits'=hits),text_size=5, fill_color = c("#F0E442", "#0072B2"), show_percentage = F))
  dev.off()

}


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
    ct_into<-cancer_match_long_into[which(cancer_match_long_CMP==ctiss)]

    if(dg %in% results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand < 0.2]){
      summary_drivers[dg,ctiss]<-"Novel"

    if(length(intersect(ct_into, inTOgen_drivers$CANCER_TYPE[inTOgen_drivers$SYMBOL==dg]))>0){
      summary_drivers[dg,ctiss]<-"Known Found"
    }
    } else if(length(intersect(ct_into, inTOgen_drivers$CANCER_TYPE[inTOgen_drivers$SYMBOL==dg]))>0){
      summary_drivers[dg,ctiss]<-"Known Not Found"
    }
  }
}

summary_drivers_bin<-summary_drivers
summary_drivers_bin[summary_drivers_bin=="Known Found"]<-1
summary_drivers_bin[summary_drivers_bin=="Novel"]<-0
class(summary_drivers_bin)<-"numeric"
save(summary_drivers_bin, file=paste(home, path_results, "summary_drivers_bin.RData", sep=""))


#for each diver genes, plotted how many tumor types were known and how many novel ones we find
selhits2<-names(which(rowSums(summary_drivers==("Known Found")|summary_drivers==("Novel"), na.rm=T)>1))

toplot<-summary_drivers_bin[selhits2,]
toplot<-toplot[order(rowSums(!is.na(toplot)), decreasing = T),order(colSums(!is.na(toplot)), decreasing = T)]

pdf(paste(home, path_results, "/exploration/figures/Drivers_pheat_novelvsknown.pdf", sep=""), 15, 15)
pheatmap(toplot, cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15)
dev.off()
graphics.off()


#for each diver genes, plotted how many tumor types were known and how many novel ones we find
selhits<-names(which(rowSums(summary_drivers==c("Known Found")|summary_drivers==c( "Novel"), na.rm=T)>1))
df<-data.frame(counts=c(rowSums(summary_drivers=="Known Found", na.rm=T), rowSums(summary_drivers=="Known Not Found", na.rm=T), rowSums(summary_drivers=="Novel", na.rm=T)),
               type=rep(c("Found in a Known Cancer Type", "Not Found in a Known Cancer Type","Found in a Novel Cancer Type"), each=nrow(summary_drivers)), gene=c(rownames(summary_drivers), rownames(summary_drivers), rownames(summary_drivers)))

df$gene<-factor(df$gene, levels=c(df$gene[df$type=="Found in a Known Cancer Type"][order(df$counts[df$type=="Found in a Known Cancer Type"]+df$counts[df$type=="Found in a Novel Cancer Type"])]))

df$counts[df$type=="Not Found in a Known Cancer Type"]<- -df$counts[df$type=="Not Found in a Known Cancer Type"]
df$type<-factor(df$type, levels=c("Found in a Known Cancer Type", "Found in a Novel Cancer Type", "Not Found in a Known Cancer Type"))


pdf(paste(home, path_results, "/exploration/figures/Drivers_barplot_thr2.pdf", sep=""),10, 5)
ggplot(subset(df, (gene %in% selhits) & (type !="Not Found in a Known Cancer Type")),aes(x=gene, y=counts, fill=type))+geom_bar(stat="identity")+theme_classic()+ylab("Number of cancer types")+xlab("")+theme(axis.text.x = element_text(angle=90, size=10))+
  scale_fill_manual(values = c("#0072B2", "#56B4E9"), labels= c("Known", "Novel"))+ guides(fill=guide_legend(title="Cancer Type"))
dev.off()

###############################
### Genes not known to be drivers: in how many cancer types
################################
hits_nodriver<-setdiff(unique(hits), driver_genes)
save(hits_nodriver, file=paste(home, path_results, "/exploration/hits_nodriver.RData", sep=""))
hits<-unique(hits)
save(hits, file=paste(home, path_results, "/exploration/hits.RData", sep=""))

summary_nodrivers<-matrix(0, nrow=length(hits_nodriver), ncol=length(tissues))
rownames(summary_nodrivers)<-hits_nodriver
colnames(summary_nodrivers)<-tissues

for(dg in hits_nodriver){
  for(ctiss in tissues){
    if(dg %in% results[[ctiss]]$GENE[results[[ctiss]]$rank_ratio<1.6 & results[[ctiss]]$medFitEff< -.5 & results[[ctiss]]$pval_rand < 0.2]){
      summary_nodrivers[dg,ctiss]<-1
    }

  }
}



df<-data.frame(counts=c(colSums(summary_drivers=="Known Found", na.rm=T), colSums(summary_drivers=="Novel", na.rm=T), colSums(summary_nodrivers==1)),
               type=rep(c("Known Driver", "Known Driver","Not Known as Driver"),each=length(tissues)), tissue=tissues)
df$counts<-df$counts/rep(mut_burden_new, 3)
df$tissue<-factor(df$tissue, levels=c(df$tissue[order(num_hits_norm, decreasing=T)]))
df_sel<-df[df$tissue %in% df$tissue[order(num_hits_norm, decreasing=T)][1:10],]

pdf(paste(home, path_results, "/exploration/figures/Hits_tissue_stack_simpl.pdf", sep=""), 5, 5)
ggplot(df_sel, aes(x=tissue, y=counts, fill=type))+geom_bar(stat="identity")+theme_classic()+ylab("% of hits")+xlab("")+theme(axis.text.x = element_text(angle=90, size=10))+
  scale_fill_manual(values = c("#0072B2",  "#E69F00"))+ guides(fill=guide_legend(title="Hit"))
dev.off()