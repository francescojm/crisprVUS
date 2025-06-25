#find patients with a DAM but without mutations in any known cancer driver

library(tidyverse)

path_data<-"/data"
path_results<-"/results/20250221"
home<-"E:/VUS_2024build"

###load patient data
COSMIC<-read_tsv(paste(home, path_data,"/raw/Cosmic_GenomeScreensMutant_Tsv_v101_GRCh38/Cosmic_GenomeScreensMutant_v101_GRCh38.tsv.gz", sep=""))
COSMIC<-data.frame(COSMIC)
COSMICsamples<-read_tsv(paste(home, path_data,"/raw/Cosmic_Sample_Tsv_v101_GRCh38/Cosmic_Sample_v101_GRCh38.tsv.gz", sep=""))
COSMICsamples<-data.frame(COSMICsamples)
COSMICclassification<-read_tsv(paste(home, path_data,"/raw/Cosmic_Classification_Tsv_v101_GRCh38/Cosmic_Classification_v101_GRCh38.tsv.gz", sep=""))
COSMICclassification<-data.frame(COSMICclassification)
COSMICprimarysite<-(COSMICclassification$PRIMARY_SITE[match(COSMIC$COSMIC_PHENOTYPE_ID, COSMICclassification$COSMIC_PHENOTYPE_ID)])


##clinical variants
#downloaded from civicdb on the 16th of May 2025 https://civicdb.org/downloads/01-May-2025/01-May-2025-ClinicalEvidenceSummaries.tsv
civicdb_var<-read_tsv("data/raw/01-May-2025-ClinicalEvidenceSummaries.tsv")
civicdb_var_sel<-civicdb_var[-grep("fusion|mutation|expression|deletion|duplication|methylation|Loss-of-function|amplification", civicdb_var$molecular_profile, ignore.case = T),]
civicdb_var_sel$molecular_profile<-gsub(" ", "-", civicdb_var_sel$molecular_profile)

civicdb_mapping<-read.csv("data/raw/CivicDB_mapping.csv", sep=";")

##load DAMs data
inTOgen_drivers<-read.table("data/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", sep='\t',stringsAsFactors = FALSE,header=TRUE)
drivers<-unique(inTOgen_drivers$SYMBOL)
act_drivers<-unique(inTOgen_drivers$SYMBOL[inTOgen_drivers$ROLE=="Act"])

load(file=paste(home, "/", path_results, "/_allDAMs.RData", sep=""))
vars<-paste(allDAMs$GENE, gsub("p.", "", allDAMs$var), sep="-")

COSMIC$MUTATION_AA<-gsub("p.", "", COSMIC$MUTATION_AA)

##to uniform COSMIC nomenclature
##remove a letter in COSMIC before fs*
COSMIC$MUTATION_AA[grep("fs\\*",COSMIC$MUTATION_AA)]<-gsub("[[:upper:]]fs","fs", COSMIC$MUTATION_AA[grep("fs\\*",COSMIC$MUTATION_AA)])
#remove uppercase letters after "del" in DAMs, but also in intogen (no del mapped)
vars[grep("del[[:upper:]]+$",vars)]<-gsub("del[[:upper:]]+$","del",vars[grep("del[[:upper:]]+$",vars)])
allDAMs$var<-vars

mapping<-read.csv(paste(home, path_data,'/raw/intOGen ctype mapping_AS.csv',sep=''),header = TRUE,row.names = 1, sep=";")

patients_drivers<-unique(COSMIC$COSMIC_SAMPLE_ID[which(COSMIC$GENE_SYMBOL %in% act_drivers)])

patients_sel<-c()
ct_sel<-c()
var_sel<-c()
gene_sel<-c()
for(i in 1:length(allDAMs$var)){
  print(i)
  gene<-allDAMs$GENE[i]
  mut<-unlist(strsplit(allDAMs$var[i], "-"))[[2]]


  ct<-strsplit(mapping[allDAMs$ctype[i],2], " \\| ")
  #COSMIC num
  ind_sel<-which(COSMIC$GENE_SYMBOL==gene & COSMIC$MUTATION_AA==mut & COSMICprimarysite %in% ct)
  patients_sel<-c(patients_sel, unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))
  ct_sel<-c(ct_sel, rep(allDAMs$ctype[i], length(unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))) )
  var_sel<-c(var_sel, rep(allDAMs$var[i], length(unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))))
  gene_sel<-c(gene_sel, rep(allDAMs$GENE[i], length(unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))))
}

patients_civic<-c()
for(i in 1:length(civicdb_var_sel$molecular_profile)){
  print(i)
  gene<-unlist(strsplit(civicdb_var_sel$molecular_profile[i], "-"))[1]
  mut<-unlist(strsplit(civicdb_var_sel$molecular_profile[i], "-"))[2]

  if(!is.na(civicdb_var_sel$disease[i])){
if (civicdb_var_sel$disease[i]=="Cancer"){
  ind_sel<-which(COSMIC$GENE_SYMBOL==gene & COSMIC$MUTATION_AA==mut)
  patients_civic<-c(patients_civic, unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))
} else {
  ct<-civicdb_mapping[match(civicdb_var_sel$disease[i],civicdb_mapping[,1]) ,2]
  #COSMIC num
  ind_sel<-which(COSMIC$GENE_SYMBOL==gene & COSMIC$MUTATION_AA==mut & COSMICprimarysite %in% ct)
  patients_civic<-c(patients_civic, unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))
}
}
}

save(patients_civic, file="patients_civic.RData")

for(i in 1:length(civicdb_var_sel$molecular_profile)){
  print(i)
  gene<-unlist(strsplit(civicdb_var_sel$molecular_profile[i], "-"))[1]
  mut<-unlist(strsplit(civicdb_var_sel$molecular_profile[i], "-"))[2]

  if(!is.na(civicdb_var_sel$disease[i])){
if (civicdb_var_sel$disease[i]=="Cancer"){
  ind_sel<-which(COSMIC$GENE_SYMBOL==gene & COSMIC$MUTATION_AA==mut)
  patients_civic<-c(patients_civic, unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))
} else {
  ct<-civicdb_mapping[match(civicdb_var_sel$disease[i],civicdb_mapping[,1]) ,2]
  #COSMIC num
  ind_sel<-which(COSMIC$GENE_SYMBOL==gene & COSMIC$MUTATION_AA==mut & COSMICprimarysite %in% ct)
  patients_civic<-c(patients_civic, unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))
}
}
}

##include patients with any mutation in variants labelled as gene + "Mutation"
civicdb_var_mut<-civicdb_var[grep("mutation", civicdb_var$molecular_profile, ignore.case = T),]
civicdb_var_mut<-civicdb_var_mut[-grep("exon", civicdb_var_mut$molecular_profile, ignore.case = T),]
civicdb_var_mut<-civicdb_var_mut[-grep("intron", civicdb_var_mut$molecular_profile, ignore.case = T),]

for(i in 1:length(civicdb_var_mut$molecular_profile)){
  print(i)
  gene<-gsub(" Mutation", "", civicdb_var_mut$molecular_profile[i],ignore.case = T)
  
  if(!is.na(civicdb_var_mut$disease[i])){
    if (civicdb_var_mut$disease[i]=="Cancer"){
      ind_sel<-which(COSMIC$GENE_SYMBOL==gene & COSMIC$MUTATION_AA==mut)
      patients_civic<-c(patients_civic, unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))
    } else {
      ct<-civicdb_mapping[match(civicdb_var_mut$disease[i],civicdb_mapping[,1]) ,2]
      #COSMIC num
      ind_sel<-which(COSMIC$GENE_SYMBOL==gene & COSMIC$MUTATION_AA==mut & COSMICprimarysite %in% ct)
      patients_civic<-c(patients_civic, unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))
    }
  }
}

save(patients_civic, file="patients_civic_all.RData")


patients_civic<-unique(patients_civic)


patients_diff<-setdiff(patients_sel, patients_civic)

df<-data.frame(patient=patients_sel, tissue=ct_sel, variant=var_sel, cooccurrent=ifelse(patients_sel %in% patients_civic, "Co-occurrent", "Non-co-occurrent"))
df$tissue<-factor(df$tissue, levels=c(names(sort(table(df$tissue), decreasing=T))))

pdf("results/20250221/co_occurrence_patients_civic.pdf", 10, 7)
ggplot(df, aes(x=tissue, fill=cooccurrent))+geom_bar()+ theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

######################
## same with DAMs that are also SAMs
#####################

load(file=paste(home, path_results, "/_DR_plots/summary_drugs.RData", sep=""))

patients_sel_drugs<-c()
ct_sel_drugs<-c()
var_sel_drugs<-c()
gene_sel_drugs<-c()
for(i in 1:length(summary_drugs$vars)){
  print(i)
  gene<-summary_drugs$genes[i]
  mut<-summary_drugs$vars[i]
  
  
  ct<-strsplit(mapping[summary_drugs$tiss[i],2], " \\| ")
  #COSMIC num
  ind_sel<-which(COSMIC$GENE_SYMBOL==gene & COSMIC$MUTATION_AA==mut & COSMICprimarysite %in% ct)
  patients_sel_drugs<-c(patients_sel_drugs, unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))
  ct_sel_drugs<-c(ct_sel_drugs, rep(summary_drugs$tiss[i], length(unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))) )
  var_sel_drugs<-c(var_sel_drugs, rep(summary_drugs$vars[i], length(unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))))
  gene_sel_drugs<-c(gene_sel_drugs, rep(summary_drugs$genes[i], length(unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel]))))
}

patients_diff<-setdiff(patients_sel_drugs, patients_civic)

df<-data.frame(patient=patients_sel_drugs, tissue=ct_sel_drugs, variant=var_sel_drugs, cooccurrent=ifelse(patients_sel_drugs %in% patients_civic, "Co-occurrent", "Non-co-occurrent"))
df$tissue<-factor(df$tissue, levels=c(names(sort(table(df$tissue), decreasing=T))))

pdf("results/20250221/co_occurrence_patients_civic_drugs.pdf", 7, 5)
ggplot(df, aes(x=tissue, fill=cooccurrent))+geom_bar()+ theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


##################
## N patients per cancer type
##################

patients_ct<-c()
for(i in 1:nrow(mapping)){
  print(i)
    ct<-strsplit(mapping[i,2], " \\| ")
  #COSMIC num
  ind_sel<-which(COSMICprimarysite %in% ct)
  patients_ct<-c(patients_ct, length(unique(COSMIC$COSMIC_SAMPLE_ID[ind_sel])))
}

names(patients_ct)<-rownames(mapping)


df<-data.frame(patient=patients_sel_drugs, tissue=ct_sel_drugs, variant=var_sel_drugs, cooccurrent=ifelse(patients_sel_drugs %in% patients_civic, "Co-occurrent", "Non-co-occurrent"))
df$tissue<-factor(df$tissue, levels=c(names(sort(table(df$tissue), decreasing=T))))

num<-table(paste(df$tissue, df$cooccurrent))
cooccurr<-rep("Co-occurrent", length(num))
cooccurr[grep("Non-co-occurrent", names(num))]<-"Non-co-occurrent"
  
df_prop<-data.frame(prop=num/patients_ct[match(gsub(" Non-co-occurrent| Co-occurrent", "", names(num)), names(patients_ct))], 
                    tissue=gsub(" Non-co-occurrent| Co-occurrent", "", names(num)),
                    cooccurrent=cooccurr)


df_prop<-df_prop[,-1]
colnames(df_prop)[1]<-"prop"
prop_tot<-by(df_prop$prop, df_prop$tissue, sum)
df_prop$tissue<-factor(df_prop$tissue, levels=names(sort(prop_tot, decreasing=T)))

pdf("results/20250221/co_occurrence_patients_civic_drugs_prop.pdf", 7, 5)
ggplot(df_prop, aes(x=tissue, fill=cooccurrent, y=prop))+geom_bar(stat="identity")+ theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




df<-data.frame(patient=patients_sel, tissue=ct_sel, variant=var_sel, cooccurrent=ifelse(patients_sel %in% patients_civic, "Co-occurrent", "Non-co-occurrent"))
df$tissue<-factor(df$tissue, levels=c(names(sort(table(df$tissue), decreasing=T))))

num<-table(paste(df$tissue, df$cooccurrent))
cooccurr<-rep("Co-occurrent", length(num))
cooccurr[grep("Non-co-occurrent", names(num))]<-"Non-co-occurrent"

df_prop<-data.frame(prop=num/patients_ct[match(gsub(" Non-co-occurrent| Co-occurrent", "", names(num)), names(patients_ct))], 
                    tissue=gsub(" Non-co-occurrent| Co-occurrent", "", names(num)),
                    cooccurrent=cooccurr)


df_prop<-df_prop[,-1]
colnames(df_prop)[1]<-"prop"
prop_tot<-by(df_prop$prop, df_prop$tissue, sum)
df_prop$tissue<-factor(df_prop$tissue, levels=names(sort(prop_tot, decreasing=T)))

pdf("results/20250221/co_occurrence_patients_civic_prop.pdf", 7, 5)
ggplot(df_prop, aes(x=tissue, fill=cooccurrent, y=prop))+geom_bar(stat="identity")+ theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

###################################
#### excluding patients with DAMs in known drivers
###################################

gene_sel_nodriver<-gene_sel[-which(gene_sel %in% drivers)]
ct_sel_nodriver<-ct_sel[-which(gene_sel %in% drivers)]
var_sel_nodriver<-var_sel[-which(gene_sel %in% drivers)]
patients_sel_nodriver<-patients_sel[-which(gene_sel %in% drivers)]

df<-data.frame(patient=patients_sel_nodriver, tissue=ct_sel_nodriver, variant=var_sel_nodriver, cooccurrent=ifelse(patients_sel_nodriver %in% patients_civic, "Co-occurrent", "Non-co-occurrent"))
df$tissue<-factor(df$tissue, levels=c(names(sort(table(df$tissue), decreasing=T))))

pdf("results/20250221/co_occurrence_patients_civic_nodriver.pdf", 10, 7)
ggplot(df, aes(x=tissue, fill=cooccurrent))+geom_bar()+ theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

df<-data.frame(patient=patients_sel_nodriver, tissue=ct_sel_nodriver, variant=var_sel_nodriver, cooccurrent=ifelse(patients_sel_nodriver %in% patients_civic, "Co-occurrent", "Non-co-occurrent"))
df$tissue<-factor(df$tissue, levels=c(names(sort(table(df$tissue), decreasing=T))))

num<-table(paste(df$tissue, df$cooccurrent))
cooccurr<-rep("Co-occurrent", length(num))
cooccurr[grep("Non-co-occurrent", names(num))]<-"Non-co-occurrent"

df_prop<-data.frame(prop=num/patients_ct[match(gsub(" Non-co-occurrent| Co-occurrent", "", names(num)), names(patients_ct))], 
                    tissue=gsub(" Non-co-occurrent| Co-occurrent", "", names(num)),
                    cooccurrent=cooccurr)


df_prop<-df_prop[,-1]
colnames(df_prop)[1]<-"prop"
prop_tot<-by(df_prop$prop, df_prop$tissue, sum)
df_prop$tissue<-factor(df_prop$tissue, levels=names(sort(prop_tot, decreasing=T)))

pdf("results/20250221/co_occurrence_patients_civic_prop_nodriver.pdf", 7, 5)
ggplot(df_prop, aes(x=tissue, fill=cooccurrent, y=prop))+geom_bar(stat="identity")+ theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##################################################################
gene_sel_drugs_nodriver<-gene_sel_drugs[-which(gene_sel_drugs %in% drivers)]
ct_sel_drugs_nodriver<-ct_sel_drugs[-which(gene_sel_drugs %in% drivers)]
var_sel_drugs_nodriver<-var_sel_drugs[-which(gene_sel_drugs %in% drivers)]
patients_sel_drugs_nodriver<-patients_sel_drugs[-which(gene_sel_drugs %in% drivers)]

df<-data.frame(patient=patients_sel_drugs_nodriver, tissue=ct_sel_drugs_nodriver, variant=var_sel_drugs_nodriver, cooccurrent=ifelse(patients_sel_drugs_nodriver %in% patients_civic, "Co-occurrent", "Non-co-occurrent"))
df$tissue<-factor(df$tissue, levels=c(names(sort(table(df$tissue), decreasing=T))))

pdf("results/20250221/co_occurrence_patients_civic_drugs_nodriver.pdf", 10, 7)
ggplot(df, aes(x=tissue, fill=cooccurrent))+geom_bar()+ theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

df<-data.frame(patient=patients_sel_drugs_nodriver, tissue=ct_sel_drugs_nodriver, variant=var_sel_drugs_nodriver, cooccurrent=ifelse(patients_sel_drugs_nodriver %in% patients_civic, "Co-occurrent", "Non-co-occurrent"))
df$tissue<-factor(df$tissue, levels=c(names(sort(table(df$tissue), decreasing=T))))

num<-table(paste(df$tissue, df$cooccurrent))
cooccurr<-rep("Co-occurrent", length(num))
cooccurr[grep("Non-co-occurrent", names(num))]<-"Non-co-occurrent"

df_prop<-data.frame(prop=num/patients_ct[match(gsub(" Non-co-occurrent| Co-occurrent", "", names(num)), names(patients_ct))], 
                    tissue=gsub(" Non-co-occurrent| Co-occurrent", "", names(num)),
                    cooccurrent=cooccurr)


df_prop<-df_prop[,-1]
colnames(df_prop)[1]<-"prop"
prop_tot<-by(df_prop$prop, df_prop$tissue, sum)
df_prop$tissue<-factor(df_prop$tissue, levels=names(sort(prop_tot, decreasing=T)))

pdf("results/20250221/co_occurrence_patients_civic_drugs_prop_nodriver.pdf", 7, 5)
ggplot(df_prop, aes(x=tissue, fill=cooccurrent, y=prop))+geom_bar(stat="identity")+ theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


