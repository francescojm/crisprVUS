#### DAMs that could not be tested as SAMs
load("results/20250221/_allHits.RData")

targets<-c()
tissues<-c()
inGDSC<-c()

for(x in 1:nrow(allHits)){
  
  target<-allHits$GENE[x]
  variant<-allHits$var[x]
  tissue<-allHits$ctype[x]
  cellLines<-unlist(str_split(allHits$ps_cl[x],', '))
  
  #find DAMs with an existing drug targeting the corresponding DAMbgs
  ids<-which(drugTargetInfo$Gene.Target==target)
  drug_ids<-drugTargetInfo$Drug.ID[ids]
  drug_names<-drugTargetInfo$Name[ids]
  
  #if there is a drug, has this been used in the GDSC and on the cell line with the variant?
  if (length(ids)>0){
    
    print(x)
    data1<-gdsc1[which(is.element(gdsc1$DRUG_ID,drug_ids) & is.element(gdsc1$SANGER_MODEL_ID,cellLines)),]
    data2<-gdsc2[which(is.element(gdsc2$DRUG_ID,drug_ids) & is.element(gdsc2$SANGER_MODEL_ID,cellLines)),]
    data1<-rbind(data1,data2)
    
    targets<-c(targets, target)
    tissues<-c(tissues, tissue)
    inGDSC<-c(inGDSC, ifelse(nrow(data1)>0, "Yes", "No"))
  
  }
}

d<-data.frame(targets, tissues, inGDSC)


tract<-read.csv("data/raw/Tractability_pacini_et_al_2024.txt", sep="\t")

pdf("results/20250221/tractability_bucket.pdf",5,7)
barplot(table(tract$min_bucket[match(unique(allDAMs$GENE), tract$id)]))
dev.off()


OT<-read.csv("data/raw/OT_targets_combined_df 1.csv", row.names=NULL)

intersect(OT$approvedSymbol, allDAMs$GENE)
intersect(OT$targetName, allDAMs$GENE)




targets<-c()
tissues<-c()
inGDSC<-c()

for(x in 1:nrow(allHits)){
  
  target<-allHits$GENE[x]
  variant<-allHits$var[x]
  tissue<-allHits$ctype[x]
  cellLines<-unlist(str_split(allHits$ps_cl[x],', '))
  
  #find DAMs with an existing drug targeting the corresponding DAMbgs
  ids<-which(OT$approvedSymbol==target)
  drug_names<-OT$prefName[ids]
  
  #if there is a drug, has this been used in the GDSC and on the cell line with the variant?
  if (length(ids)>0){
    
    print(x)
    data1<-gdsc1[which(is.element(toupper(gdsc1$DRUG_NAME),toupper(drug_names)) & is.element(gdsc1$SANGER_MODEL_ID,cellLines)),]
    data2<-gdsc2[which(is.element(toupper(gdsc2$DRUG_NAME),toupper(drug_names)) & is.element(gdsc2$SANGER_MODEL_ID,cellLines)),]
    data1<-rbind(data1,data2)
    
    targets<-c(targets, target)
    tissues<-c(tissues, tissue)
    inGDSC<-c(inGDSC, ifelse(nrow(data1)>0, "Yes", "No"))
    
  }
}

dOT<-data.frame(targets, tissues, inGDSC)


#using the union of Goncalves and OT
drug_anno<-rbind.data.frame(data.frame(drug=drugTargetInfo$Name, target=drugTargetInfo$Gene.Target),
                            data.frame(drug=OT$prefName, target=OT$approvedSymbol))
  
targets<-c()
tissues<-c()
inGDSC<-c()

for(x in 1:nrow(allHits)){
  
  target<-allHits$GENE[x]
  variant<-allHits$var[x]
  tissue<-allHits$ctype[x]
  cellLines<-unlist(str_split(allHits$ps_cl[x],', '))
  
  #find DAMs with an existing drug targeting the corresponding DAMbgs
  ids<-which(drug_anno$target==target)
  drug_names<-drug_anno$drug[ids]
  
  #if there is a drug, has this been used in the GDSC and on the cell line with the variant?
  if (length(ids)>0){
    
    print(x)
    data1<-gdsc1[which(is.element(toupper(gdsc1$DRUG_NAME),toupper(drug_names)) & is.element(gdsc1$SANGER_MODEL_ID,cellLines)),]
    data2<-gdsc2[which(is.element(toupper(gdsc2$DRUG_NAME),toupper(drug_names)) & is.element(gdsc2$SANGER_MODEL_ID,cellLines)),]
    data1<-rbind(data1,data2)
    
    targets<-c(targets, target)
    tissues<-c(tissues, tissue)
    inGDSC<-c(inGDSC, ifelse(nrow(data1)>0, "Yes", "No"))
    
  }
}

dall<-data.frame(targets, tissues, inGDSC)


