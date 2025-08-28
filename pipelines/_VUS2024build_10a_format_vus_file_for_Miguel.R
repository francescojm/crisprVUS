
library(data.table)
library(stringr)
library(biomaRt)


####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250808_bugFixed_and_RR_th.1.71_wr/'

cl_variants <- read_csv(paste(pathdata,'/raw/mutations_all_20241212.csv', sep=""))
### mutations_all_20230202 downloaded from  https://cog.sanger.ac.uk/cmp/download/mutations_all_20241212.zip on 20250221

load(paste(resultPath,'_allDAMs.RData',sep=''))


variant_signature<-paste(cl_variants$gene_symbol,cl_variants$protein_mutation)
variant_signature1<-paste(allDAMs$GENE,allDAMs$var)

ii<-which(is.element(variant_signature,variant_signature1))
variant_signature<-variant_signature[ii]
cl_variants <- cl_variants[ii,]

uvariant_signature<-unique(variant_signature)

cl_variants<-cl_variants[match(uvariant_signature,variant_signature),]

cl_variants<-cl_variants[order(cl_variants$gene_symbol),]

cl_variants<-as.data.frame(cl_variants)
cl_variants<-cl_variants[,-c(3,4,5,11,12,13,14,16,17,18,19)]
cl_variants<-cl_variants[,c(1,3,2,6,7,5,4,8)]

colnames(cl_variants)<-c('gene_name','protein_mutation','gene_id','chromosome','position','cdna_mutation','rna_mutation','effect')

dt <- as.data.table(cl_variants)  # your table

# 1) seq_region_name
dt[, seq_region_name := str_replace(chromosome, "^chr", "")]

# 2) start/end (SNVs)
dt[, `:=`(start = position, end = position)]

# 3) allele_string from cDNA HGVS like "c.1799T>A"
#    If not present, try rna_mutation similarly; fall back to NA
parse_allele <- function(hgvs) {
  if (is.na(hgvs)) return(NA_character_)
  m <- str_match(hgvs, "c\\.[0-9+_-]+([ACGT])>([ACGT])")
  if (!is.na(m[1,1])) return(paste0(m[1,2], "/", m[1,3]))
  # simple del/ins/dup examples (you’ll want VEP for full reliability)
  m <- str_match(hgvs, "c\\.([0-9_]+)del([ACGT]+)?")
  if (!is.na(m[1,1])) return(paste0(ifelse(is.na(m[1,3]), "", m[1,3]), "/-"))
  m <- str_match(hgvs, "c\\.([0-9_]+)ins([ACGT]+)")
  if (!is.na(m[1,1])) return(paste0("-/", m[1,3]))
  NA_character_
}
dt[, allele_string := ifelse(!is.na(cdna_mutation),
                             vapply(cdna_mutation, parse_allele, character(1)),
                             vapply(rna_mutation,  parse_allele, character(1)))]

# 4) strand from Ensembl gene_id (GRCh38)
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = NULL)
strand_map <- getBM(
  attributes = c("ensembl_gene_id", "strand"),
  filters    = "ensembl_gene_id",
  values     = unique(dt$gene_id),
  mart       = mart
)
dt <- merge(dt, as.data.table(strand_map), by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE)

allDAMs_Positions<-as.data.frame(dt)

allDAMs_Positions<-allDAMs_Positions[order(allDAMs_Positions$gene_name),]

allDAMs_Positions<-allDAMs_Positions[,c(2,3,1,4,5,6,7,8,9,10,11,13,12)]

rownames(allDAMs_Positions)<-NULL


save(allDAMs_Positions,file=paste(resultPath,'/_allDAMs_Positions.RData',sep=''))

write.table(allDAMs_Positions,file=paste(resultPath,'/_allDAMs_Positions.tsv',sep=''),sep='\t',quote=FALSE,row.names = FALSE)

