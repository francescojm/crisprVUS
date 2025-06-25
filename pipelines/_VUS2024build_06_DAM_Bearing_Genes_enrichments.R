library(ReactomePA)
library(biomaRt)
library(clusterProfiler)
library(reactome.db)
library(AnnotationDbi)


####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250221/'

my.hypTest<-function(x,k,n,N){
  
  PVALS<-phyper(x-1,n,N-n,k,lower.tail=FALSE)
  
  return(PVALS)
}

# loading all DAM bearing genes and inTOgen drivers
load(paste(resultPath,'/_allDAM_bearing_genes.RData',sep=''))
inTOgen_drivers<-read.table(paste(pathdata,"/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", sep=""), sep='\t',stringsAsFactors = FALSE,header=TRUE)
inTOgen_drivers<-unique(sort(inTOgen_drivers$SYMBOL))

#defining the background
load(paste(resultPath,'/_totalTestedVariants.RData',sep=''))
background<-unique(totalTestedVariants$gene_symbol)

#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",mirror = "asia")

# Converting all dam bearing genes into entrez id
gene_symbols <- unique(allDAM_bearing_genes$allDAM_bearing_genes)
allDAM_bearing_entrez<- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "hgnc_symbol",
  values = gene_symbols,
  mart = ensembl
)

# Converting all background genes into entrez id
background_entrez<- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "hgnc_symbol",
  values = background,
  mart = ensembl
)

allDAM_bearing_entrez<-as.character(allDAM_bearing_entrez$entrezgene_id)
background_entrez<-as.character(background_entrez$entrezgene_id)

# Performing Enrichment analysis of reactome pathways in the DAM bearing genes (all of them)
all_DAM_enrichments <- enrichPathway(gene=allDAM_bearing_entrez, pvalueCutoff = 0.05, readable=TRUE, universe=background_entrez)

pdf(paste(resultPath, "exploration/figures/REACTOME_allDAMbearing_AS2.pdf",sep=""), 7,6)
dotplot(all_DAM_enrichments,showCategory=10)
dev.off()

write.table(all_DAM_enrichments,quote=FALSE,sep='\t',dec=",", row.names = F, file=paste(resultPath,'/_pathEnrichments_all_DAMbearing_genes_AS2.txt',sep=''))

# selecting only new DAM bearing genes
new_gene_symbols <- setdiff(unique(allDAM_bearing_genes$allDAM_bearing_genes),inTOgen_drivers)
new_DAM_bearing_entrez <- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "hgnc_symbol",
  values = new_gene_symbols,
  mart = ensembl
)
new_DAM_bearing_entrez<-as.character(new_DAM_bearing_entrez$entrezgene_id)

# Performing Enrichment analysis of reactome pathways in the DAM bearing genes (new only of them)
new_DAM_enrichments <- enrichPathway(gene=new_DAM_bearing_entrez, pvalueCutoff = 0.05, readable=TRUE, universe=background_entrez)

write.table(new_DAM_enrichments,quote=FALSE,sep='\t',dec=",",row.names = F, file=paste(resultPath,'/_pathEnrichments_new_DAMbearing_genes_AS2.txt',sep=''))


IntoGEN_Drivers_entrez<-result <- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "hgnc_symbol",
  values = inTOgen_drivers,
  mart = ensembl
)
IntoGEN_Drivers_entrez<-as.character(IntoGEN_Drivers_entrez$entrezgene_id)

# selecting only known DAM bearing genes
known_DAM_bearing_entrez<-intersect(allDAM_bearing_entrez,IntoGEN_Drivers_entrez)

# Performing Enrichment analysis of reactome pathways in the DAM bearing genes (new only of them)
known_DAM_enrichments <- enrichPathway(gene=known_DAM_bearing_entrez, pvalueCutoff = 0.05, readable=TRUE, universe=background_entrez)

write.table(known_DAM_enrichments,quote=FALSE,sep='\t',dec=",",row.names = F, file=paste(resultPath,'/_pathEnrichments_known_DAMbearing_genes_AS2.txt',sep=''))

known_DAM_enrichments<-as.data.frame(known_DAM_enrichments)
all_DAM_enrichments<-as.data.frame(all_DAM_enrichments)

allpaths<-intersect(known_DAM_enrichments$Description,all_DAM_enrichments$Description)

all_genes <- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  mart = ensembl
)

all_genes<-all_genes$entrezgene_id
all_genes<-all_genes[!is.na(all_genes)]

N<-length(all_genes)

######genes in the background not known to be drivers
new_background<-setdiff(background, inTOgen_drivers)
# Converting all background genes into entrez id
new_background_entrez<- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters = "hgnc_symbol",
  values = new_background,
  mart = ensembl
)

# uncomment to recompute cleng_
pathway_to_genes <- as.list(reactomePATHID2EXTID)

cleng_<-NULL

n_new_DAM_entrez_in_Reactome<-length(intersect(new_DAM_bearing_entrez,unlist(pathway_to_genes)))
n_new_DAM_entrez_out_Reactome<-length(setdiff(new_DAM_bearing_entrez,unlist(pathway_to_genes)))

allGenes_in_reactome<-intersect(all_genes,unlist(pathway_to_genes))
allGenes_out_reactome<-setdiff(all_genes,unlist(pathway_to_genes))

set.seed(1234)

for (i in 1:1000){
   print(i)

   randomGenesInPathways<-as.character(sample(intersect(new_background_entrez, allGenes_in_reactome),n_new_DAM_entrez_in_Reactome))
   randomGenesOutPathways<-as.character(sample(intersect(new_background_entrez, allGenes_out_reactome),n_new_DAM_entrez_out_Reactome))

   randomGenes<-c(randomGenesInPathways,randomGenesOutPathways)
   random_DAM_enrichments <- enrichPathway(gene=c(known_DAM_bearing_entrez,
                                                  randomGenes),
                                                  pvalueCutoff = 0.05, readable=TRUE, universe=background_entrez)
   commonRandpaths<-intersect(known_DAM_enrichments$Description,random_DAM_enrichments$Description)
   cleng_[i]<-length(commonRandpaths)
   hist(cleng_)
   print(length(which(cleng_>=length(allpaths)))/i)
 }
 save(cleng_,file=paste(resultPath,'/_all_vs_known_DAM_bearingG_enrichedPathways_overlap_random_AS2.RData',sep=''))
# uncomment to recompute cleng_
 
load(file=paste(resultPath,'/_all_vs_known_DAM_bearingG_enrichedPathways_overlap_random_AS2.RData',sep=''))


hist(cleng_)
abline(v=length(allpaths),col='red')

pval<-length(which(cleng_>=length(allpaths)))/length(cleng_)

print(paste('probability of having',length(allpaths),
            'enriched pathways conserved when adding to the DAM bearing genes known to be cancer deriver genes an additional',
      length(new_DAM_bearing_entrez),'genes selected by random chance =',pval))


pp1<-rep(1,length(allpaths))
pp2<-rep(1,length(allpaths))

pp1<-known_DAM_enrichments$p.adjust[match(allpaths,known_DAM_enrichments$Description)]
names(pp1)<-allpaths
pp2<-all_DAM_enrichments$p.adjust[match(allpaths,all_DAM_enrichments$Description)]
names(pp2)<-allpaths




has_a_cancer_driver<-unlist(lapply(pathway_to_genes,
                                   function(x){length(intersect(x,IntoGEN_Drivers_entrez))>0}
)
)
all_genes <- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  mart = ensembl
)

all_genes<-all_genes$entrezgene_id
all_genes<-all_genes[!is.na(all_genes)]

N<-length(all_genes)

genes_in_pathways_with_a_cancer_driver<-
  unique(unlist(pathway_to_genes[which(has_a_cancer_driver)]))

genes_in_pathways_with_a_cancer_driver<-intersect(genes_in_pathways_with_a_cancer_driver,all_genes)

new_DAM_bearing_in_a_pathway_with_a_cancer_driver<-
  intersect(new_DAM_bearing_entrez,genes_in_pathways_with_a_cancer_driver)

n<-length(genes_in_pathways_with_a_cancer_driver)
k<-length(new_DAM_bearing_entrez)
x<-length(new_DAM_bearing_in_a_pathway_with_a_cancer_driver)


res<-NULL
for (i in 1:1000){
  nn<-sample(all_genes,length(new_DAM_bearing_entrez))  
  res[i]<-length(intersect(nn,genes_in_pathways_with_a_cancer_driver))
}

my.hypTest(x,k,n,N)


Path_increasedCoverage<-
  cbind(known_DAM_enrichments[match(allpaths,known_DAM_enrichments$Description),'Count'],
        all_DAM_enrichments[match(allpaths,all_DAM_enrichments$Description),'Count'])

rownames(Path_increasedCoverage)<-allpaths

Path_increasedCoverage<-Path_increasedCoverage[Path_increasedCoverage[,2]-Path_increasedCoverage[,1]>0,]

oo<-order(Path_increasedCoverage[,2]-Path_increasedCoverage[,1])
pdf(paste(resultPath, "exploration/figures/path_increasedcoverage_AS2.pdf",sep=""), 10,15)
par(mar=c(4,26,0,0.5))
barplot(t(cbind(Path_increasedCoverage[oo,2]-Path_increasedCoverage[oo,1],Path_increasedCoverage[oo,2])),
        beside = FALSE,horiz = TRUE,las=2,xlab='n. genes',xlim=c(0,130),cex.names=0.5,border=FALSE,
        col=c('purple','orange'))
dev.off()


