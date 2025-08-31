set.seed(123)
library(CELLector)
library(tidyverse)
library(readr)
library(data.table)

####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250808_bugFixed_and_RR_th.1.71_wr/'


# loading all DAM bearing genes and inTOgen drivers
load(paste(resultPath,'/_allDAM_bearing_genes.RData',sep=''))
inTOgen_drivers<-read.table(paste(pathdata,"/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", sep=""), sep='\t',stringsAsFactors = FALSE,header=TRUE)
inTOgen_drivers<-unique(sort(inTOgen_drivers$SYMBOL))


f <- file.path(pathdata, "raw", "OpenTargets_tractability_20250829.tsv")
## downloaded from http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/input/target/tractability/tractability.tsv (on the 29 of August 2025)

tract <- data.table::fread(
  f, sep = "\t", header = TRUE, quote = "", fill = TRUE,
  showProgress = FALSE
)


load(paste(resultPath,'_allDAMs_with_SIFT_PolyPhen_scores.RData',sep=''))
load(paste(resultPath,'_totalTestedVariants.RData',sep=''))

allDAMs<-allDAMs_with_SIFT_Polyphen

DAM_unique_ensIds<-
  unique(totalTestedVariants$ensembl_gene_id[match(allDAMs$GENE,totalTestedVariants$gene_symbol)])

DAM_unique_genes<-
  totalTestedVariants$gene_symbol[match(DAM_unique_ensIds,totalTestedVariants$ensembl_gene_id)]

tractability<-tract[match(DAM_unique_ensIds,tract$ensembl_gene_id),
      c('Top_bucket_sm','Top_bucket_ab')]

knownDAMbgs<-is.element(DAM_unique_genes,inTOgen_drivers)

TractabilityTable<-cbind(DAM_unique_ensIds,DAM_unique_genes,knownDAMbgs,tractability)
minBucket<-apply(TractabilityTable[,c(4,5)],MARGIN = 1,FUN = 'min',na.rm=TRUE)

Tractability_tiers <- cut(
  minBucket,
  breaks = c(0, 3, 6, 10),   # boundaries
  labels = c(1, 2, 3),       # bin labels
  right = TRUE,              # intervals closed on the right
  include.lowest = TRUE
)

TractabilityTable<-cbind(TractabilityTable,minBucket,Tractability_tiers)

BAR1<-summary(as.factor(TractabilityTable$knownDAMbgs[which(TractabilityTable$Tractability_tiers==1)]))
BAR2<-summary(as.factor(TractabilityTable$knownDAMbgs[which(TractabilityTable$Tractability_tiers==2)]))
BAR3<-summary(as.factor(TractabilityTable$knownDAMbgs[which(TractabilityTable$Tractability_tiers==3)]))

toPlot<-rbind(BAR1,BAR2,BAR3)

colnames(toPlot)<-c('unreported','known')


DAMactionability<-TractabilityTable$Tractability_tiers[match(allDAMs$GENE,TractabilityTable$DAM_unique_genes)]

ii<-which(allDAMs$impact_call!='Low impact' &
            allDAMs$impact_call!='Uknwon impact')

allDAMs<-allDAMs[ii,]

DAMactionability<-DAMactionability[ii]

DAMactionabilityKnown<-DAMactionability[which(allDAMs$impact_call!='Low impact' &
                                           allDAMs$impact_call!='Uknwon impact' &
                                             is.element(allDAMs$GENE,inTOgen_drivers))]


DAMactionabilityUnknown<-DAMactionability[which(allDAMs$impact_call!='Low impact' &
                                                allDAMs$impact_call!='Uknwon impact' &
                                                !is.element(allDAMs$GENE,inTOgen_drivers))]


