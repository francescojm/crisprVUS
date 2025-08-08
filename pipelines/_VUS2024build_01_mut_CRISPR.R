library(sets)
library(binaryLogic)#installed from github
library(tidyverse)
library(CoRe)

RR_th<-1.71

##functions called in the following code
variantSpectrum<-function(cl_var,gene){
  
  tmp<-cl_var[cl_var$gene_symbol_2023==gene,c('protein_mutation','gene_symbol_2023')]
  
  aa<-sort(unique(tmp$protein_mutation))
  aa<-setdiff(aa, c("-", "p.?"))
  
  aa<-tmp[match(aa,tmp$protein_mutation),]
  
  return(aa)
}

positives<-function(cl_var,gene,variants){
  return(sort(unique(cl_var$model_id[cl_var$gene_symbol_2023==gene & is.element(cl_var$protein_mutation,variants)])))
}

clasTest<-function(ess_scores,bess_scores,cl_var,gene,vs,vs_cds,display=TRUE,save=NULL){
  
  ps_cl<-positives(cl_var,gene,variants=vs)
  ps_cl<-intersect(ps_cl,colnames(ess_scores))
  
  FCp<-unlist(ess_scores[gene,])
  hits<-match(ps_cl,names(FCp)[order(FCp)])
  
  n<-length(hits)
  
  if (sum(hits)>0){
    dependent_cls<-names(which(bess_scores[gene,]>0))
    
    JI_essMut<-length(intersect(dependent_cls,ps_cl))/length(union(dependent_cls,ps_cl))
    rankRatio<-sum(hits)/((n*(n+1))/2)
    medEff<-median(FCp[ps_cl])
    
    if(length(ps_cl)!=length(FCp)){
      mostDepNegative<-sort(FCp[setdiff(names(FCp),ps_cl)])[1]  
    }else{
      mostDepNegative<-NA
    }
    
    if (is.element(gene,rownames(basal_exp))){
      expPattern<-basal_exp[gene,intersect(colnames(basal_exp),names(FCp))]
      avgEXP_fpkm_in_positive_cl<-mean(expPattern[ps_cl])
      
      if(sum(!is.na(expPattern))>=10){
        cdf <- ecdf(expPattern)
        perc_basal_exp_of_Avg_positive_cl<-round(100*cdf(mean(expPattern[ps_cl]))) 
      }else{
        
        DD<-sort(expPattern)
        testval<-mean(expPattern[ps_cl])
        perc_basal_exp_of_Avg_positive_cl<-100*max(which(DD<=testval))/length(DD)
      }
    }else{
      avgEXP_fpkm_in_positive_cl<-NA
      perc_basal_exp_of_Avg_positive_cl<-NA
    }
    
            RES<- data.frame(GENE=gene,var=c(paste(vs,collapse=' | ')),var_cs=c(paste(vs_cds,collapse=' | ')),
                           npos=length(ps_cl),
                     medFitEff=medEff,
                     leap=medEff-mostDepNegative,
                     totcl=length(FCp),
                     rank_ratio = rankRatio,
                     matching=JI_essMut,
                     ps_cl=paste(ps_cl,collapse=', '),highest_rank=max(hits),
                     avgBasalEXP_fpkm_in_ps_cl=avgEXP_fpkm_in_positive_cl,
                     percBasalEXP_of_ps_cl=perc_basal_exp_of_Avg_positive_cl,
                     stringsAsFactors = FALSE)
    
    if (display){
      FCp<-sort(FCp)
      
      pcpatt<-rep(16,length(FCp))
      pcpatt[which(bess_scores[gene,names(FCp)]>0)]<-18
      
      bg<-rep(rgb(0,0,255,alpha = 110,maxColorValue = 255),length(FCp))
      names(bg)<-names(FCp)
      
      bg[ps_cl]<-'red'
      
      #if(rankRatio==1 & JI_essMut>0 & medEff< -1){
        if(length(save)>0){
          pdf(paste(save,gene,'_',paste(gsub("\\?|\\*|!|>", "", vs),collapse='AND'),'.pdf',sep=''),7.60,7.20)
        }
        
        par(mfrow=c(2,1))
        plot(FCp,col=bg,pch=pcpatt,xlab='cell lines',
             ylab=paste(gene,'scaled fitness effect'),
             main=paste('Rank ratio =',format(sum(hits)/((n*(n+1))/2),digits=3),', mut/ess match = ',
                                                                 format(100*JI_essMut,digits=3),'%'))
        abline(h= -1,col='gray')
        abline(h= -0.5,col='gray',lty=2)
        legend('bottomright',pch=c(15,18,16),col=c('red','gray','gray'),
               legend = c(paste(vs,collapse=' | '),'significant fitness effect','non significant fitness effect'),bg = 'white')
        
        hist(FCp,100,col=rgb(0,0,255,alpha = 110,maxColorValue = 255),border=NA,xlab=paste(gene,'fitness effect'),main='')
        
        points(FCp[ps_cl],rep(0,length(ps_cl)),col='red',pch=15,cex=1.5)
        abline(v= -1,col='gray')
        abline(v= -0.5,col='gray',lty=2)
        
        if(length(save)>0){
          dev.off()
        } 
      }    
    #}
    
  }
  
  
  
  return(RES)
  
}

Optimal_clasTest<-function(ess_scores,bess_scores,cl_var,gene,vs,vs_cds,display=TRUE,save=NULL){
  
  if(length(vs)>1){
    vs_subset_idx<-do.call(rbind,as.binary(1:(2^length(vs)-1),n = length(vs)))+0  
    
    res<-do.call(rbind,lapply(1:nrow(vs_subset_idx),function(x){
      
      curr_var<-vs[which(vs_subset_idx[x,]>0)]
      curr_var_cds<-vs_cds[which(vs_subset_idx[x,]>0)]
      clasTest(ess_scores = ess_scores,
               bess_scores = bess_scores,
               cl_var = cl_var,
               gene = gene,
               vs=curr_var,
               vs_cds = curr_var_cds,
               display = FALSE)
    }))
    
    #res$rank_ratio[res$rank_ratio<RR_th]<-1 #serviva perché altrimenti si selezionano sempre varianti presenti in una sola linea
    if(length(which(res$rank_ratio<RR_th & res$medFitEff))==0){
      res<-res[which(res$rank_ratio==min(res$rank_ratio)),]
    } else {
      res<-res[which(res$rank_ratio<RR_th & res$medFitEff< -.5),]
      #res<-res[which(res$rank_ratio<RR_th,]
    }
    
    if(nrow(res)>1){
      ind<-which.max(unlist(lapply(res$var, function(x){length(unlist(str_split(x," \\| ")))})))
      res<-res[ind,]
    }
  }else{
    vs_subset_idx<-1
    
    res<-
      clasTest(ess_scores = ess_scores,
               bess_scores = bess_scores,
               cl_var = cl_var,
               gene = gene,vs=vs,vs_cds = vs_cds,display = FALSE)
    }
  
  
  if(res$rank_ratio<RR_th & res$medFitEff< -.5){
    curr_var<-setdiff(unlist(str_split(res$var,'[ | ]')),'')
    curr_var_cds<-setdiff(unlist(str_split(res$var_cds,'[ | ]')),'')
    
    clasTest(ess_scores = ess_scores,
             bess_scores = bess_scores,
             cl_var = cl_var,
             gene = gene,
             vs=curr_var,vs_cds = curr_var_cds,display = TRUE,save=save)
  }
  
  
  return(res) 
}

getCOSMICfreqs<-function(COSMICf,GENE,Mutation_CDS,Mutation_AA){
  RES<-COSMICf[which(COSMICf$GENE_NAME==GENE & COSMICf$`Mutation CDS`==Mutation_CDS & COSMICf$`Mutation AA`==Mutation_AA),]

  tmp<-unlist(str_split(RES$DISEASE,';'))
  tmp<-str_split(tmp,'=')
  
  disease<-unlist(lapply(tmp,function(x){x[1]}))
  
  n_positives<-unlist(lapply(str_split(unlist(lapply(tmp,function(x){x[2]})),'/'),function(x){as.numeric(x[1])}))
  n_tested<-unlist(lapply(str_split(unlist(lapply(tmp,function(x){x[2]})),'/'),function(x){as.numeric(x[2])}))
  
  percs<-round(100*n_positives/n_tested,2)
  
  
  
  tmp_wgs<-unlist(str_split(RES$WGS_DISEASE,';'))
  tmp_wgs<-str_split(tmp_wgs,'=')
  
  disease_wgs<-unlist(lapply(tmp_wgs,function(x){x[1]}))
  
  n_positives_wgs<-unlist(lapply(str_split(unlist(lapply(tmp_wgs,function(x){x[2]})),'/'),function(x){as.numeric(x[1])}))
  n_tested_wgs<-unlist(lapply(str_split(unlist(lapply(tmp_wgs,function(x){x[2]})),'/'),function(x){as.numeric(x[2])}))
  
  percs_wgs<-round(100*n_positives_wgs/n_tested_wgs,2)
  
  if(!is.na(tmp_wgs[1])){
    disease<-c(disease,disease_wgs)
    n_positives<-c(n_positives,n_positives_wgs)
    n_tested<-c(n_tested,n_tested_wgs)
    percs<-c(percs,percs_wgs)
    
  }
  
  
  ncas<-length(percs)
  
  
  RES<-data.frame(
    GENE=rep(GENE,ncas),
    ROLE = rep(RES$ONC_TSG, ncas),
    CDS = rep(RES$`Mutation CDS`, ncas),
    AA = rep(RES$`Mutation AA`, ncas),
    desc_CDS = rep(RES$`Mutation Description CDS`, ncas),
    desc_AA = rep(RES$`Mutation Description AA`, ncas),
    N_POS_TOTAL = rep(RES$COSMIC_SAMPLE_MUTATED,ncas),
    N_TESTED_TOTAL = rep(RES$COSMIC_SAMPLE_TESTED,ncas),
    DISEASE = disease,
    nMutant = n_positives,
    nTested = n_tested,
    perc = percs, stringsAsFactors = FALSE)
  
  RES<-RES[order(RES$perc,decreasing=TRUE),]
  
  return(RES)
  
}

multipleGeneClasTests<-function(ess_scores,bess_scores,cl_var,genes,path){
  
  genes<-intersect(genes,cl_var$gene_symbol_2023)
  ngenes<-length(genes)
  
  resTOT<-NULL
  for (i in 1:ngenes){
    print(paste(genes[i],i,'of',ngenes))
    vs<-variantSpectrum(cl_var = cl_var,gene = genes[i])
    vs_cds<-vs$cDNA_mutation
    vs<-vs$protein_mutation
    curRES<-Optimal_clasTest(ess_scores = ess_scores,
                             bess_scores = bess_scores,
                             cl_var = cl_var,
                             gene = genes[i],
                             vs = vs,
                             vs_cds = vs_cds,
                             save = path)
    
    ## attach hypergeometric test p
    curRES<-cbind(curRES,hypTest_p=my.hypTest(x=curRES$npos,k = curRES$highest_rank,n = curRES$npos,N = curRES$totcl))
    resTOT<-rbind(resTOT,curRES)
    if(length(resTOT)>0){resTOT<-resTOT[order(resTOT$rank_ratio),]}
  }
  
  return(resTOT)
}
my.hypTest<-function(x,k,n,N){
  
  PVALS<-phyper(x-1,n,N-n,k,lower.tail=FALSE)
  
  return(PVALS)
}

###### ANALYSIS

# select the tissue and filter the cell lines accordingly
ts_depFC<-scaled_depFC[,CMP_annot$model_id[CMP_annot$cancer_type==ctiss]]
ts_bdep<-bdep[,CMP_annot$model_id[CMP_annot$cancer_type==ctiss]]

ts_cl_variants<-cl_variants[which(is.element(cl_variants$model_id,colnames(ts_depFC))),]

# excluding non-coding mutations
ts_cl_variants<-ts_cl_variants[which(ts_cl_variants$coding),]
# excluding start_lost mutations
ts_cl_variants<-ts_cl_variants[which(ts_cl_variants$effect!='start_lost'),]
# excluding silent mutations
ts_cl_variants<-ts_cl_variants[which(ts_cl_variants$effect!='silent'),]
# excluding nonsense mutations 
ts_cl_variants<-ts_cl_variants[which(ts_cl_variants$effect!='nonsense'),]
# excluding stop_lost mutations 
ts_cl_variants<-ts_cl_variants[which(ts_cl_variants$effect!='stop_lost'),]
# excluding ess_splice mutations 
ts_cl_variants<-ts_cl_variants[which(ts_cl_variants$effect!='ess_splice'),]

# select the genes to analyze (i.e. those with cancer dependency data available)
ts_cl_variants<-ts_cl_variants[which(is.element(ts_cl_variants$gene_symbol_2023,rownames(ts_depFC))),]

# select the genes to analyze
genesToTest<-unique(ts_cl_variants$gene_symbol_2023)

#counts the number of observed syntactically different mutations per gene
vs_spec_cardinality<-unlist(lapply(lapply(1:length(genesToTest),
                                          function(x){
                                            dd<-variantSpectrum(cl_var = ts_cl_variants,gene = genesToTest[x])
                                            dd<-dd$protein_mutation}),'length'))

names(vs_spec_cardinality)<-genesToTest

##### exclude highly mutated genes (e.g. TTN and Tumor suppressor genes) 
genesToTest<-sort(names(which(vs_spec_cardinality<=10 & vs_spec_cardinality>0)))

##### exclude pan-cancer core fitness genes and common-essential genes
genesToTest<-setdiff(genesToTest,ADaM)
genesToTest<-setdiff(genesToTest,Perc_AUC)

#create output folder
if(!dir.exists(paste(resultPath,'/_CRISPR_plots/',sep=''))){
  dir.create(paste(resultPath,'/_CRISPR_plots/',sep=''))  
}

dir.create(paste(resultPath,'/_CRISPR_plots/',ctiss,sep=''))

#run main function and generate a data.frame with DAMs
RESTOT<-multipleGeneClasTests(ess_scores = ts_depFC,
                              bess_scores = ts_bdep,
                              cl_var = ts_cl_variants,
                              genes = genesToTest,
                              path = paste(resultPath,'/_CRISPR_plots/',ctiss,'/',sep=''))

RESTOT<-data.frame(ctype=rep(ctiss,nrow(RESTOT)),RESTOT,stringsAsFactors = FALSE)

#annotate genes with info on core fitness genes
data(curated_BAGEL_essential)
cfg<-CoRe.ADaM(ts_bdep,TruePositives = curated_BAGEL_essential,display = FALSE,verbose = FALSE)

is_cts_cfg<-is.element(RESTOT$GENE,cfg)
RESTOT<-cbind(RESTOT,is_cts_cfg)

save(RESTOT,file=paste(resultPath,'/',ctiss,'_results.RData',sep=''))
save(ts_cl_variants,file=paste(resultPath,'/',ctiss,'_testedVariants.RData',sep=''))
write.table(RESTOT,quote=FALSE,sep='\t',row.names = FALSE,file=paste(resultPath,'/',ctiss,'_results.tsv',sep=''))
write.table(ts_cl_variants,quote=FALSE,sep='\t',row.names = FALSE,file=paste(resultPath,'/',ctiss,'_testedVariants.tsv',sep=''))

####Save all genes with DAMs in an RData object
DAM_bearing_genes<-unique(RESTOT$GENE[RESTOT$rank_ratio<RR_th & RESTOT$medFitEff< -.5 ])
DAM_bearing_genes<-unique(DAM_bearing_genes)

save(DAM_bearing_genes, file=paste(resultPath,'/',ctiss,'_DAM_bearing_genes.RData',sep=''))

