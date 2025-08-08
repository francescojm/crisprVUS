
inTOgen_drivers<-read.table("data/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", sep='\t',stringsAsFactors = FALSE,header=TRUE)
load('results/20250221/_allHits.RData')
load('results/20250221/_allDAM_bearing_genes.RData')
load('results/20250221/_allDAMs.RData')
driver_genes<-inTOgen_drivers$SYMBOL

unreportedDAMs<-allDAMs[allDAMs$GENE %in% setdiff(allDAMs$GENE, driver_genes),]
knownDAMs<-allDAMs[allDAMs$GENE %in% driver_genes,]

unreportedDAMs$type<-"unreported"
knownDAMs$type<-"known"
DAMs<-rbind(knownDAMs,unreportedDAMs)

write.csv(unreportedDAMs, file="unreportedDAMs.csv")
write.csv(knownDAMs, file="knownDAMs.csv")
write.csv(DAMs, file="DAMs.csv")

library(tidyverse)

####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250221/'

###loading input data
gene_annot <- read_csv(paste(pathdata, "/raw/gene_identifiers_20241212.csv", sep=""))
### gene_identifiers_20191101 downloaded from https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20241212.csv on 20250221

CMP_annot <- read_csv(paste(pathdata,"/raw/model_list_20250630.csv", sep="")) 
### model_list_20250630.csv downloaded from https://cog.sanger.ac.uk/cmp/download/model_list_20250630.csv on 20250724

cl_variants <- read_csv(paste(pathdata,'/raw/mutations_all_20250318.csv', sep=""))
### mutations_all_20250318 downloaded from  https://cog.sanger.ac.uk/cmp/download/mutations_all_20250318.zip on 20250724

cl_variants <- cbind(cl_variants,gene_annot$hgnc_symbol[match(cl_variants$gene_id,gene_annot$gene_id)])
colnames(cl_variants)[ncol(cl_variants)]<-'gene_symbol_2023'


########organoids
organoid_id<-CMP_annot$model_id[CMP_annot$model_type=="Organoid"]

org_variants<-cl_variants[cl_variants$model_id %in% organoid_id,]
org_variants$tissue<-CMP_annot$cancer_type[match(org_variants$model_id, CMP_annot$model_id)]

##look for organoids with a DAM in matching cancer type
string_tosearch_org<-paste(org_variants$gene_symbol_2023, org_variants$protein_mutation, org_variants$tissue, sep="\t")
string_tosearch_unrDAMs<-paste(unreportedDAMs$GENE, unreportedDAMs$var, unreportedDAMs$ctype, sep="\t")
string_tosearch_knDAMs<-paste(knownDAMs$GENE, knownDAMs$var, knownDAMs$ctype, sep="\t")

found_unr<-intersect(string_tosearch_org, string_tosearch_unrDAMs)
found_kn<-intersect(string_tosearch_org, string_tosearch_knDAMs)

organoid_ind<-match(c(found_unr, found_kn), string_tosearch_org)
unr_ind<-match(c(found_unr), string_tosearch_unrDAMs)
kn_ind<-match(c(found_kn), string_tosearch_knDAMs)
status<-c(rep("Unreported", length(found_unr)), rep("Known", length(found_kn)))
  
found<-data.frame(gene=c(unreportedDAMs$GENE[unr_ind], knownDAMs$GENE[kn_ind]),
                  var=c(unreportedDAMs$var[unr_ind], knownDAMs$var[kn_ind]),
                  cancer_type=c(unreportedDAMs$ctype[unr_ind], knownDAMs$ctype[kn_ind]),
                  status=c(unreportedDAMs$type[unr_ind], knownDAMs$type[kn_ind]),
                  organoid_id=org_variants$model_id[organoid_ind])

library(openxlsx)
write.xlsx(found, "organoids_with_DAMs.xlsx")
