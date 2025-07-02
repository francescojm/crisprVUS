library(ggplot2)
library(openxlsx)
library(tidyverse)
library(STRINGdb)


path_data<-"/data"
path_results<-"/results/20250221"
#home<-"E:/VUS_2024build"
home<-"/Users/francesco.iorio/Dropbox/CODING/vus/VUS 2025/"

setwd(paste(home, "/", path_results, sep=""))

#load("hits.RData")
####### I couldn't fine Hits.Rdata #########
load("_allHits.RData")
hits<-unique(allHits$GENE)

##driver genes
inTOgen_drivers<-read.table(paste(home, '/data/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv', sep=""),sep='\t',stringsAsFactors = FALSE,header=TRUE)
driver_genes<-unique(inTOgen_drivers$SYMBOL)
hits_nodriver<-setdiff(hits, driver_genes)


##get STRING interactions
string_db <- STRINGdb$new( version="12", species=9606,
                           score_threshold=200, input_directory="")

#compute an interaction score between DAM-bearing genes and driver genes by summing all STRING scores
df<-data.frame(gene=unique(c(hits_nodriver, driver_genes)))
Genesmapped <- string_db$map( df, "gene", removeUnmappedRows = TRUE )
int<-string_db$get_interactions(Genesmapped$STRING_id)
score<-sum(int[int[,1] %in% Genesmapped[Genesmapped[,1] %in% hits_nodriver,2] & int[,2] %in% Genesmapped[Genesmapped[,1] %in% driver_genes,2],3])+
  sum(int[int[,2] %in% Genesmapped[Genesmapped[,1] %in% hits_nodriver,2] & int[,1] %in% Genesmapped[Genesmapped[,1] %in% driver_genes,2],3])

library(tidyverse)

## empirical pvalue: sample random genes in the background of all tested genes not in drivers and compute the STRING interaction score with known drivers
load("_totalTestedVariants.RData")
all_genes_nodriver<-setdiff(unique(totalTestedVariants$gene_symbol), driver_genes)

set.seed(84905750)
scorerand<-c()

  for(j in 1:1000){
    print(j)
    rand_genes<-sample(all_genes_nodriver, length(hits_nodriver), replace = F)
    dfrand<-data.frame(gene=unique(c(rand_genes, driver_genes)))
    Randmapped <- string_db$map( dfrand, "gene", removeUnmappedRows = TRUE )
    intrand<-string_db$get_interactions(Randmapped$STRING_id)
    scorerand<-c(scorerand, sum(intrand[intrand[,1] %in% Randmapped[Randmapped[,1] %in% rand_genes,2] & intrand[,2] %in% Randmapped[Randmapped[,1] %in% driver_genes,2],3])+
      sum(intrand[intrand[,2] %in% Randmapped[Randmapped[,1] %in% rand_genes,2] & intrand[,1] %in% Randmapped[Randmapped[,1] %in% driver_genes,2],3]))
    hist(scorerand,xlim=c(min(c(scorerand,score)),max(c(scorerand,score))))
    abline(v=score)
    print(length(which(scorerand>=score))/j)
  }

save(scorerand, file="scorerand_novel_FI.RData")

df<-data.frame(First_neighbours=scorerand)

pdf("PPI_neigghbours_novel_FI.pdf",5,4)
ggplot(df, aes(x=First_neighbours))+geom_density()+
  geom_vline(xintercept = score, linetype="dashed", 
        color = "red", size=1)+theme_classic()

dev.off()




