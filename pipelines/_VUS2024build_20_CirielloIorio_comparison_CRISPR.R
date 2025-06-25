library(ggplot2)
library(ggvenn)
resultPath<-'results/20250221/'
home<-"E:/VUS_2024build"
load(paste(resultPath,'_totalTestedVariants.RData',sep=''))
load(paste(resultPath,'_allHits.RData',sep=''))
load(paste(resultPath,'_allDAMs.RData',sep=''))
load(paste(resultPath,'_allDAM_bearing_genes.RData',sep=''))

load("data/Ciriello/mut_summary_merged_4comparisonDAMs.RData")
mutC<-mut_summary_merged
#####intersection of tested variants
Iorio_tested<-unique(paste(totalTestedVariants$gene_symbol, gsub("p.", "", totalTestedVariants$protein_mutation)))
Ciriello_tested<-unique(paste(mutC$gene, mutC$mutation)[which(!is.na(mutC$level_CRISPR))])
#fs variants are labelled differently, change Ter in * as in Iorios
Ciriello_tested<-gsub("Ter","*", Ciriello_tested)
##remove a letter in Ciriello before fs*
Ciriello_tested[grep("fs\\*",Ciriello_tested)]<-gsub("[[:upper:]]fs","fs", Ciriello_tested[grep("fs\\*",Ciriello_tested)])
##inframe deletions: in Ciriello not indicated the aminoacid deleted after "del", only before
#remove uppercase letters after "del" in Iorios
Iorio_tested[grep("del[[:upper:]]+$",Iorio_tested)]<-gsub("del[[:upper:]]+$","del",Iorio_tested[grep("del[[:upper:]]+$",Iorio_tested)])

##frameshift with fs*1 in Iorio terminate with * in Ciriello as the meaning of the notation is:
#p.R83SfsX1	arginine (R) is the first amino acid changed, it is in position 83, it makes serine (S) instead, the length of the shift frame is 1, including the stop codon (X)
#equivalent to R83*, as the only included "AA" is the stop codon
Iorio_tested[grep("fs\\*1$", Iorio_tested)]<-gsub("fs\\*1$", "\\*", Iorio_tested[grep("fs\\*1$", Iorio_tested)])



tested_ic<-intersect(Iorio_tested,Ciriello_tested)
Ciriello_unique<-setdiff(Ciriello_tested, Iorio_tested)
Iorio_unique<-setdiff(Iorio_tested, Ciriello_tested)

df<-data.frame(class=c("Iorio", "Ciriello", "Intersection"),
               num=c(length(Iorio_tested), length(Ciriello_tested), length(tested_ic)))

pdf(paste(resultPath, "Background_variants_comparison.pdf", sep="/"), 8.5,1.5)
ggplot(df, aes(x=class, y=num))+geom_bar(stat="identity")+theme_classic()+theme(axis.text=element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))+ coord_flip()
dev.off()

#######################
### Enrichment test at the level of variants DEP_I
##########################
Iorio_DAMs<-unique(paste(allDAMs$GENE, gsub("p.", "", allDAMs$var)))
Ciriello_DAMs<-unique(mutC[which(mutC$ann_CRISPR_refined=="Dep_I"),])
Ciriello_DAMs<-paste(Ciriello_DAMs$gene, Ciriello_DAMs$mutation)
#fs variants are labelled differently, change Ter in * as in Iorios
Ciriello_DAMs<-gsub("Ter","*", Ciriello_DAMs)
##remove a letter in Ciriello before fs*
Ciriello_DAMs[grep("fs\\*",Ciriello_DAMs)]<-gsub("[[:upper:]]fs","fs", Ciriello_DAMs[grep("fs\\*",Ciriello_DAMs)])
##inframe deletions: in Ciriello not indicated the aminoacid deleted after "del", only before
#remove uppercase letters after "del" in Iorios
Iorio_DAMs[grep("del[[:upper:]]+$",Iorio_DAMs)]<-gsub("del[[:upper:]]+$","del",Iorio_DAMs[grep("del[[:upper:]]+$",Iorio_DAMs)])

##frameshift with fs*1 in Iorio terminate with * in Ciriello as the meaning of the notation is:
#p.R83SfsX1	arginine (R) is the first amino acid changed, it is in position 83, it makes serine (S) instead, the length of the shift frame is 1, including the stop codon (X)
#equivalent to R83*, as the only included "AA" is the stop codon
Iorio_DAMs[grep("fs\\*1$", Iorio_DAMs)]<-gsub("fs\\*1$", "\\*", Iorio_DAMs[grep("fs\\*1$", Iorio_DAMs)])

Iorio_DAMs_NI<-setdiff(Iorio_DAMs, tested_ic)
Ciriello_DAMs_NI<-setdiff(Ciriello_DAMs, tested_ic)

##DAMs must be in the common set of tested variants
Iorio_DAMs<-intersect(Iorio_DAMs, tested_ic)
Ciriello_DAMs<-intersect(Ciriello_DAMs, tested_ic)

bothIC<-intersect(Iorio_DAMs, Ciriello_DAMs)

library(openxlsx)
write.xlsx(data.frame(var=bothIC), "CirielloIorioshared.xlsx")
df<-data.frame(class=c("Iorio", "Ciriello"),
               shared_background=c("Shared", "Shared", "Non-shared", "Non-shared"),
               num=c(length(Iorio_DAMs), length(Ciriello_DAMs), length(Iorio_DAMs_NI), length(Ciriello_DAMs_NI)))

pdf(paste(resultPath, "Num_variants_comparison.pdf", sep="/"), 10,1.2)
ggplot(df, aes(x=class, y=num, fill=shared_background, label=num))+geom_bar(stat="identity")+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+theme_classic()+theme(axis.text=element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))+ coord_flip()
dev.off()

library(ggvenn)
pdf(paste(resultPath, "Venn_comparison.pdf", sep=""), 5, 5)
ggvenn(data=list('Iorio'=Iorio_DAMs, 'Ciriello'=Ciriello_DAMs),text_size=5, fill_color = c("#F0E442", "#0072B2"), show_percentage = F)
dev.off()

###enrichment for drivers
a<-length(intersect(Iorio_DAMs, Ciriello_DAMs))
b<-length(setdiff(Iorio_DAMs, Ciriello_DAMs))
c<-length(setdiff(Ciriello_DAMs, Iorio_DAMs))
d<-length(setdiff(tested_ic,intersect(Iorio_DAMs, Ciriello_DAMs)))

fisher.test(matrix(c(a,b,c,d), ncol=2), alternative="greater")

###do we both find mutations not known to be drivers?
##driver genes
inTOgen_drivers<-read.table(paste(home, '/data/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
### inTOgene drivers downloaded from https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip on 20241002
driver_genes<-unique(inTOgen_drivers$SYMBOL)

genes_bothIC<-unique(unlist(strsplit(bothIC, " "))[seq(1,length(bothIC)*2,2)])
paste("Number of gene hits known as driver:", length(intersect(genes_bothIC, driver_genes)))
intersect(genes_bothIC, driver_genes)
paste("Number of gene hits not known as driver:", length(setdiff(genes_bothIC, driver_genes)))
setdiff(genes_bothIC, driver_genes)

df<-data.frame(class=c("Driver", "Non-driver"),
               num=c(length(intersect(genes_bothIC, driver_genes)), length(setdiff(genes_bothIC, driver_genes))))

pdf(paste(resultPath, "Num_variants_comparison_driver.pdf", sep="/"), 7, 1.5)
ggplot(df, aes(x=class, y=num))+geom_bar(stat="identity")+theme_classic()+theme(axis.text=element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))+ coord_flip()
dev.off()


all_genes_bothIC<-(unlist(strsplit(bothIC, " "))[seq(1,length(bothIC)*2,2)])
driver_status<-ifelse(all_genes_bothIC %in% driver_genes, "Driver", "Non-driver")
df<-data.frame(gene=all_genes_bothIC,
               driver_status=driver_status)
df$gene<-factor(df$gene, levels=names(sort(table(df$gene), decreasing=T)))

pdf(paste(resultPath, "Num_occurrences.pdf", sep=""), 7, 5)
ggplot(df, aes(x=gene, fill=driver_status))+geom_bar()+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))+scale_fill_manual(values = c("#0072B2",  "#E69F00"))
dev.off()


#######################################################
#### GENES
##########################################################

#####intersection of tested genes
Iorio_tested<-unique(totalTestedVariants$gene_symbol)
Ciriello_tested<-unique(mutC$gene)
tested_ic<-intersect(Iorio_tested,Ciriello_tested)
Ciriello_unique<-setdiff(Ciriello_tested, Iorio_tested)
Iorio_unique<-setdiff(Iorio_tested, Ciriello_tested)

df<-data.frame(class=c("Iorio", "Ciriello", "Intersection"),
               num=c(length(Iorio_tested), length(Ciriello_tested), length(tested_ic)))

pdf(paste(resultPath, "Background_genes_comparison.pdf", sep="/"), 7,2)
ggplot(df, aes(x=class, y=num))+geom_bar(stat="identity")+theme_classic()+theme(axis.text=element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))+ coord_flip()
dev.off()

#######################
### Enrichment test at the level of variants
##########################
Iorio_DAMs<-unique(allDAMs$GENE)
Ciriello_DAMs<-mutC[which(mutC$ann_CRISPR_refined=="Dep_I"),]
Ciriello_DAMs<-unique(Ciriello_DAMs$gene)

##DAMs must be in the common set of tested variants
Iorio_DAMs<-intersect(Iorio_DAMs, tested_ic)
Ciriello_DAMs<-intersect(Ciriello_DAMs, tested_ic)

bothIC<-intersect(Iorio_DAMs, Ciriello_DAMs)

df<-data.frame(class=c("Iorio", "Ciriello"),
               shared_background=c("Shared", "Shared", "Non-shared", "Non-shared"),
               num=c(length(Iorio_DAMs), length(Ciriello_DAMs), length(Iorio_DAMs_NI), length(Ciriello_DAMs_NI)))

pdf(paste(resultPath, "Num_variants_comparison_gene.pdf", sep="/"), 7,1.5)
ggplot(df, aes(x=class, y=num, fill=shared_background))+geom_bar(stat="identity")+theme_classic()+theme(axis.text=element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))+ coord_flip()
dev.off()

library(ggvenn)
pdf(paste(resultPath, "Venn_comparison_gene.pdf", sep=""), 5, 5)
ggvenn(data=list('Iorio'=Iorio_DAMs, 'Ciriello'=Ciriello_DAMs),text_size=5, fill_color = c("#F0E442", "#0072B2"), show_percentage = F)
dev.off()

###enrichment for drivers
a<-length(intersect(Iorio_DAMs, Ciriello_DAMs))
b<-length(setdiff(Iorio_DAMs, Ciriello_DAMs))
c<-length(setdiff(Ciriello_DAMs, Iorio_DAMs))
d<-length(setdiff(tested_ic,intersect(Iorio_DAMs, Ciriello_DAMs)))

fisher.test(matrix(c(a,b,c,d), ncol=2), alternative="greater")

###do we both find mutations not known to be drivers?
##driver genes
inTOgen_drivers<-read.table(paste(home, '/data/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
### inTOgene drivers downloaded from https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip on 20241002
driver_genes<-unique(inTOgen_drivers$SYMBOL)

paste("Number of gene hits known as driver:", length(intersect(bothIC, driver_genes)))
intersect(bothIC, driver_genes)
paste("Number of gene hits not known as driver:", length(setdiff(bothIC, driver_genes)))
setdiff(bothIC, driver_genes)

df<-data.frame(class=c("Driver", "Non-driver"),
               num=c(length(intersect(bothIC, driver_genes)), length(setdiff(bothIC, driver_genes))))

pdf(paste(resultPath, "Num_variants_comparison_driver_gene.pdf", sep="/"), 7,1.5)
ggplot(df, aes(x=class, y=num))+geom_bar(stat="identity")+theme_classic()+theme(axis.text=element_text(size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))+ coord_flip()
dev.off()


driver_status<-ifelse(bothIC %in% driver_genes, "Driver", "Non-driver")
df<-data.frame(gene=bothIC,
               driver_status=driver_status)
df$gene<-factor(df$gene, levels=names(sort(table(df$gene), decreasing=T)))

pdf(paste(resultPath, "Num_occurrences_gene.pdf", sep=""), 7, 5)
ggplot(df, aes(x=driver_status))+geom_bar(stat="count")+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12), legend.title = element_text(size=12), legend.text = element_text(size=12))+
  scale_fill_manual(values = c("#0072B2",  "#E69F00"))+ coord_flip()
dev.off()
