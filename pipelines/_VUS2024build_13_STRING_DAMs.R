library(ggplot2)
library(openxlsx)
library(tidyverse)
library(STRINGdb)

path_data<-"/data"
path_results<-"/results/20250221"
home<-"E:/VUS_2024build"

setwd(paste(home, "/", path_results, sep=""))
load("hits.RData")

##driver genes
inTOgen_drivers<-read.table(paste(home, '/data/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
driver_genes<-unique(inTOgen_drivers$SYMBOL)

##get STRING interactions
string_db <- STRINGdb$new( version="12", species=9606,
                           score_threshold=200, input_directory="")

#compute an interaction score between DAM-bearing genes and driver genes by summing all STRING scores
df<-data.frame(gene=unique(c(hits, driver_genes)))
Modmapped <- string_db$map( df, "gene", removeUnmappedRows = TRUE )
int<-string_db$get_interactions(Modmapped$STRING_id)
scoremod<-sum(int[int[,1] %in% Modmapped[Modmapped[,1] %in% hits,2] & int[,2] %in% Modmapped[Modmapped[,1] %in% driver_genes,2],3])+
  sum(int[int[,2] %in% Modmapped[Modmapped[,1] %in% hits,2] & int[,1] %in% Modmapped[Modmapped[,1] %in% driver_genes,2],3])

## empirical pvalue: sample random genes in the background of all tested genes and compute the STRING interaction score with known drivers
library(tidyverse)
load("E:/VUS_2024build/results/20250221/_totalTestedVariants.RData")
all_genes<-unique(totalTestedVariants$gene_symbol)

set.seed(84905750)
scorerand<-c()

  for(j in 1:1000){
    print(j)
    rand_genes<-sample(all_genes, length(hits), replace = F)
    dfrand<-data.frame(gene=unique(c(rand_genes, driver_genes)))
    Randmapped <- string_db$map( dfrand, "gene", removeUnmappedRows = TRUE )
    intrand<-string_db$get_interactions(Randmapped$STRING_id)
    scorerand<-c(scorerand, sum(intrand[intrand[,1] %in% Randmapped[Randmapped[,1] %in% rand_genes,2] & intrand[,2] %in% Randmapped[Randmapped[,1] %in% driver_genes,2],3])+
      sum(intrand[intrand[,2] %in% Randmapped[Randmapped[,1] %in% rand_genes,2] & intrand[,1] %in% Randmapped[Randmapped[,1] %in% driver_genes,2],3]))

  }

save(scorerand, file="scorerand.RData")

df<-data.frame(First_neighbours=scorerand)

pdf("PPI_neigghbours.pdf",5,4)
ggplot(df, aes(x=First_neighbours))+geom_density()+ geom_vline(xintercept = scoremod, linetype="dashed", 
                                                               color = "red", size=1)+theme_classic()

dev.off()

