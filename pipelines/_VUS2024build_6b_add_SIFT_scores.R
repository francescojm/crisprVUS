library(readr)
library(data.table)
library(stringr)
library(biomaRt)
library(httr)
library(jsonlite)



####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250808_bugFixed_and_RR_th.1.71_wr/'

# loading all DAM bearing genes and inTOgen drivers
load(paste(resultPath,'/_allDAM_bearing_genes.RData',sep=''))
inTOgen_drivers<-read.table(paste(pathdata,"/raw/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv", sep=""), sep='\t',stringsAsFactors = FALSE,header=TRUE)
inTOgen_drivers<-unique(sort(inTOgen_drivers$SYMBOL))

# #Uncomment to recompute SIFT and PolyPhen scores
# load(paste(resultPath,'_allDAMs.RData',sep=''))
# 
# cl_variants <- read_csv(paste(pathdata,'/raw/mutations_all_20241212.csv', sep=""))
# ### mutations_all_20230202 downloaded from  https://cog.sanger.ac.uk/cmp/download/mutations_all_20241212.zip on 20250221
# 
# cl_variants<-cl_variants[which(!is.na(cl_variants$position)),]
# 
# variant_signature<-paste(cl_variants$gene_symbol,cl_variants$protein_mutation)
# variant_signature1<-paste(allDAMs$GENE,allDAMs$var)
# 
# ii<-which(is.element(variant_signature,variant_signature1))
# variant_signature<-variant_signature[ii]
# cl_variants <- cl_variants[ii,]
# 
# uvariant_signature<-unique(variant_signature)
# 
# cl_variants<-cl_variants[match(uvariant_signature,variant_signature),]
# 
# cl_variants<-cl_variants[order(cl_variants$gene_symbol),]
# 
# dt <- as.data.table(cl_variants)  # replace with your object
# 
# ## 1) seq_region_name
# dt[, seq_region_name := str_replace(chromosome, "^chr", "")]
# 
# ## 2) start / end (SNVs). For indels, adjust end if you wish.
# dt[, `:=`(start = as.integer(position),
#           end   = as.integer(position))]
# 
# ## 3) allele_string from REF/ALT
# dt[, allele_string := ifelse(!is.na(reference) & !is.na(alternative),
#                              paste0(reference, "/", alternative), NA_character_)]
# 
# ## 4) strand via Ensembl (GRCh38 by default)
# ens_ids <- unique(na.omit(ifelse(!is.na(dt$ensembl_gene_id), dt$ensembl_gene_id, dt$gene_id)))
# 
# mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# strand_map <- getBM(
#   attributes = c("ensembl_gene_id","strand"),
#   filters    = "ensembl_gene_id",
#   values     = ens_ids,
#   mart       = mart
# )
# strand_map <- as.data.table(strand_map)
# dt <- strand_map[dt, on = c("ensembl_gene_id" = "ensembl_gene_id")]
# setnames(dt, "strand", "strand_numeric")
# dt[, strand := as.integer(strand_numeric)][, strand_numeric := NULL]
# 
# ## 5) final table
# out <- dt[, .(gene_symbol, model_id, seq_region_name, start, end, strand, allele_string)]
# head(out)
# 
# dt<-dt[,-c(4,5,13,14)]
# 
# vep_sift_polyphen_safe <- function(df,
#                                    batch_size  = 50,
#                                    max_retries = 3,
#                                    sleep_sec   = 1,
#                                    verbose     = TRUE) {
#   # ---- Input checks & sanitisation -----------------------------------------
#   req_cols <- c("seq_region_name","start","end","reference","alternative")
#   stopifnot(all(req_cols %in% names(df)))
#   
#   x <- as.data.table(df)
#   x <- x[!is.na(seq_region_name) & !is.na(start) & !is.na(end) &
#            !is.na(reference) & !is.na(alternative)]
#   x[, `:=`(start = as.integer(start), end = as.integer(end))]
#   x <- x[start > 0 & end >= start]
#   
#   if (!nrow(x)) {
#     warning("No valid variants after basic filtering.")
#     return(list(results = data.table(), failed = data.table()))
#   }
#   
#   mk_variant <- function(chr, st, en, ref, alt) sprintf("%s %s %s %s/%s", chr, st, en, ref, alt)
#   x[, vep_input := mk_variant(seq_region_name, start, end, reference, alternative)]
#   
#   batches   <- split(x, ceiling(seq_len(nrow(x))/batch_size))
#   n_batches <- length(batches)
#   
#   if (verbose) {
#     message(sprintf("Submitting %d variants in %d batch(es) to Ensembl VEP REST…",
#                     nrow(x), n_batches))
#     pb <- utils::txtProgressBar(min = 0, max = n_batches, style = 3)
#     on.exit(try(close(pb), silent = TRUE), add = TRUE)
#   }
#   
#   # ---- Helpers --------------------------------------------------------------
#   parse_one <- function(rec) {
#     # Return NULL if no transcript consequences
#     if (is.null(rec$transcript_consequences)) return(NULL)
#     tc <- rec$transcript_consequences
#     dt <- rbindlist(lapply(tc, as.data.frame), fill = TRUE)
#     
#     # Carry over input + core locus fields
#     dt[, `:=`(input          = rec$input,
#               seq_region_name= rec$seq_region_name,
#               start          = rec$start,
#               end            = rec$end,
#               allele_string  = rec$allele_string)]
#     dt
#   }
#   
#   request_batch <- function(dt_batch, depth = 0, batch_id = NA_integer_) {
#     if (nrow(dt_batch) == 0) return(list(ok = NULL, bad = NULL))
#     
#     variants <- as.list(dt_batch$vep_input)
#     attempt  <- 0
#     
#     repeat {
#       attempt <- attempt + 1
#       if (verbose) {
#         msg <- sprintf("Batch %s | size=%d | attempt=%d%s",
#                        ifelse(is.na(batch_id), "?", batch_id),
#                        nrow(dt_batch), attempt,
#                        if (depth > 0) sprintf(" | depth=%d (bisect)", depth) else "")
#         message(msg)
#       }
#       
#       res <- try({
#         POST(
#           url = "https://rest.ensembl.org/vep/human/region",
#           add_headers(`Content-Type` = "application/json", `Accept` = "application/json"),
#           body = toJSON(list(variants = variants), auto_unbox = TRUE),
#           timeout(60)
#         )
#       }, silent = TRUE)
#       
#       # Connection error
#       if (inherits(res, "try-error")) {
#         if (verbose) message("  -> HTTP request failed (connection).")
#         if (attempt <= max_retries) { Sys.sleep(sleep_sec * attempt); next } else break
#       }
#       
#       code <- status_code(res)
#       
#       if (code == 200) {
#         ans <- content(res, as = "parsed", type = "application/json")
#         ok  <- rbindlist(lapply(ans, parse_one), fill = TRUE)
#         return(list(ok = ok, bad = NULL))
#       }
#       
#       if (code %in% c(429, 503)) {
#         if (verbose) message(sprintf("  -> Server says %d; backing off…", code))
#         if (attempt <= max_retries) { Sys.sleep(sleep_sec * attempt); next } else break
#       }
#       
#       if (code >= 500) {
#         if (verbose) message(sprintf("  -> HTTP %d; bisecting batch to isolate offending records…", code))
#         if (nrow(dt_batch) == 1 || depth > 10) {
#           if (verbose) message("  -> Reached single record or max depth; marking as failed.")
#           return(list(ok = NULL, bad = dt_batch))
#         }
#         mid   <- nrow(dt_batch) %/% 2
#         left  <- request_batch(dt_batch[1:mid], depth + 1, batch_id)
#         right <- request_batch(dt_batch[(mid+1):nrow(dt_batch)], depth + 1, batch_id)
#         return(list(ok = rbindlist(list(left$ok, right$ok), fill = TRUE),
#                     bad = rbindlist(list(left$bad, right$bad), fill = TRUE)))
#       }
#       
#       # Other client error (400 etc.) — mark whole batch as failed
#       if (verbose) message(sprintf("  -> HTTP %d (client error); marking batch as failed.", code))
#       return(list(ok = NULL, bad = dt_batch))
#     }
#     
#     # Exit without success
#     list(ok = NULL, bad = dt_batch)
#   }
#   
#   # ---- Process batches ------------------------------------------------------
#   ok_list <- list(); bad_list <- list()
#   for (i in seq_along(batches)) {
#     ans <- request_batch(batches[[i]], batch_id = i)
#     if (!is.null(ans$ok))  ok_list[[length(ok_list)+1]]  <- ans$ok
#     if (!is.null(ans$bad)) bad_list[[length(bad_list)+1]] <- ans$bad
#     if (verbose) utils::setTxtProgressBar(pb, i)
#   }
#   
#   if (verbose) message("\nParsing and summarising transcript-level predictions…")
#   resdt <- rbindlist(ok_list, fill = TRUE)
#   
#   # ---- Robust field harmonisation & summarisation ---------------------------
#   if (nrow(resdt)) {
#     # Harmonise gene symbol naming
#     if (!"gene_symbol" %in% names(resdt) && "hgnc_symbol" %in% names(resdt)) {
#       resdt[, gene_symbol := hgnc_symbol]
#     }
#     # Ensure consequence_terms exists
#     if (!"consequence_terms" %in% names(resdt)) {
#       if ("consequence_term" %in% names(resdt)) {
#         resdt[, consequence_terms := consequence_term]
#       } else if ("consequence" %in% names(resdt)) {
#         resdt[, consequence_terms := consequence]
#       } else {
#         resdt[, consequence_terms := NA_character_]
#       }
#     }
#     # Ensure prediction/score columns exist
#     for (cn in c("sift_score","sift_prediction","polyphen_score","polyphen_prediction",
#                  "gene_symbol","transcript_id")) {
#       if (!cn %in% names(resdt)) resdt[[cn]] <- NA
#     }
#     
#     # Build join key and merge back to inputs
#     resdt[, join_key := sprintf("%s %s %s %s", seq_region_name, start, end, allele_string)]
#     x[, join_key := sprintf("%s %s %s %s/%s", seq_region_name, start, end, reference, alternative)]
#     
#     keep <- resdt[, .(join_key, gene_symbol, transcript_id, consequence_terms,
#                       sift_score, sift_prediction, polyphen_score, polyphen_prediction)]
#     
#     suppressWarnings({
#       keep[, sift_score := as.numeric(sift_score)]
#       keep[, polyphen_score := as.numeric(polyphen_score)]
#     })
#     
#     summary <- keep[
#       , .(sift_min          = if (all(is.na(sift_score))) NA_real_ else min(sift_score, na.rm = TRUE),
#           sift_pred_any     = paste(na.omit(unique(sift_prediction)), collapse = "|"),
#           polyphen_max      = if (all(is.na(polyphen_score))) NA_real_ else max(polyphen_score, na.rm = TRUE),
#           polyphen_pred_any = paste(na.omit(unique(polyphen_prediction)), collapse = "|")),
#       by = join_key
#     ]
#     
#     annotated <- merge(x, summary, by = "join_key", all.x = TRUE)
#   } else {
#     annotated <- x
#     annotated[, c("sift_min","sift_pred_any","polyphen_max","polyphen_pred_any") := NA]
#   }
#   
#   failed <- rbindlist(bad_list, fill = TRUE)
#   
#   if (verbose) {
#     message(sprintf("Done. Annotated: %d | Failed: %d", nrow(annotated), nrow(failed)))
#     if (nrow(failed)) {
#       message("A few examples of failed inputs:")
#       print(head(failed[, .(seq_region_name, start, end, reference, alternative)], 5))
#       message("Common causes: wrong build, REF/ALT mismatch, large/complex indels.")
#     }
#   }
#   
#   list(
#     results = annotated[, !c("join_key","vep_input")],
#     failed  = failed
#   )
# }
# 
# 
# annotated <- vep_sift_polyphen_safe(dt)
# 
# allDAMs_sigs<-paste(allDAMs$GENE,allDAMs$var)
# 
# eff_sigs<-paste(annotated$results$gene_symbol,annotated$results$protein_mutation)
# 
# nDAMs<-length(allDAMs_sigs)
# 
# ii<-match(allDAMs_sigs,eff_sigs)
# 
# allDAMs_with_SIFT_Polyphen<-cbind(allDAMs,annotated$results[ii,c('sift_min',
#                                           'sift_pred_any',
#                                           'polyphen_max',
#                                           'polyphen_pred_any')])
# 
# 
# 
# `%||%` <- function(a,b) if (is.null(a) || is.na(a)) b else a
# has_label <- function(x, lab) grepl(sprintf("\\b%s\\b", lab), x %||% "", ignore.case = TRUE)
# 
# is_sift_del   <- function(s) !is.na(s) & s < 0.05
# is_sift_tol   <- function(s) !is.na(s) & s >= 0.05
# 
# is_pp_prob    <- function(p) !is.na(p) & p >= 0.909
# is_pp_poss    <- function(p) !is.na(p) & p >= 0.446 & p < 0.909
# is_pp_benign  <- function(p) !is.na(p) & p <= 0.445
# 
# classify_variant <- function(sift_min, sift_pred_any, polyphen_max, polyphen_pred_any) {
#   # numeric cues
#   s_del <- is_sift_del(sift_min)
#   s_tol <- is_sift_tol(sift_min)
#   
#   p_prob <- is_pp_prob(polyphen_max)
#   p_poss <- is_pp_poss(polyphen_max)
#   p_ben  <- is_pp_benign(polyphen_max)
#   
#   # text cues (any transcript label)
#   s_del_l <- has_label(sift_pred_any, "deleterious")
#   s_tol_l <- has_label(sift_pred_any, "tolerated")
#   
#   p_prob_l <- has_label(polyphen_pred_any, "probably")
#   p_poss_l <- has_label(polyphen_pred_any, "possibly")
#   p_ben_l  <- has_label(polyphen_pred_any, "benign")
#   
#   # collapse signals
#   any_del <- s_del || s_del_l || p_prob || p_prob_l || p_poss || p_poss_l
#   any_ben <- s_tol || s_tol_l || p_ben || p_ben_l
#   
#   # UNKNOWN: no numbers and no labels at all
#   if (all(is.na(sift_min), is.na(polyphen_max)) &&
#       !any(nzchar(sift_pred_any), nzchar(polyphen_pred_any))) {
#     return("Unknown impact")
#   }
#   
#   # HIGH / MODERATE / POSSIBLE impact strata
#   if ((s_del || s_del_l) && (p_prob || p_prob_l))            return("High impact")
#   if ((s_del || s_del_l) && (p_poss || p_poss_l))            return("Moderate impact")
#   if ((p_prob || p_prob_l) && !(s_del || s_del_l))           return("Moderate impact")
#   if ((s_del || s_del_l) && !(p_prob || p_prob_l))           return("Possible impact")
#   if ((p_poss || p_poss_l) && !(s_del || s_del_l))           return("Possible impact")
#   
#   # LOW: at least one benign/tolerated cue and no damaging cues
#   if (any_ben && !any_del)                                   return("Low impact")
#   
#   # Fallback when mixed weak cues exist but not matching above
#   "Unknown impact"
# }
# 
# # Apply
# allDAMs_with_SIFT_Polyphen$impact_call <- mapply(classify_variant,
#                                                  allDAMs_with_SIFT_Polyphen$sift_min,
#                                                  allDAMs_with_SIFT_Polyphen$sift_pred_any,
#                                                  allDAMs_with_SIFT_Polyphen$polyphen_max,
#                                                  allDAMs_with_SIFT_Polyphen$polyphen_pred_any)
# 
# 
# save(allDAMs_with_SIFT_Polyphen,file=paste(resultPath,'_allDAMs_with_SIFT_PolyPhen_scores.RData',sep=''))
# write.table(allDAMs_with_SIFT_Polyphen,sep="\t",quote=FALSE,file=paste(resultPath,'_allDAMs_with_SIFT_PolyPhen_scores.tsv',sep=''))
# #[END] Uncomment to recompute SIFT and PolyPhen scores


load(paste(resultPath,'_allDAMs_with_SIFT_PolyPhen_scores.RData'))

VEPs<-summary(as.factor(allDAMs_with_SIFT_Polyphen$impact_call))
VEPs<-VEPs[c('High impact','Moderate impact','Possible impact','Low impact','Unknown impact')]

pdf(paste(resultPath,'/_figures_source/DAM_VEP_prediction.pdf',sep=''),11,6)
pie(VEPs,col=c("#800026","#fc4e2a","#feb24c","#ffeda0","grey85"))
dev.off()

paste(sum(VEPs[1:2]),'DAMs have a high/moderate functional impact') 
paste(sum(VEPs[1:3]),'DAMs have a high/moderate/possible functional impact') 


VEPs<-summary(as.factor(allDAMs_with_SIFT_Polyphen$impact_call[
  which(!is.element(allDAMs_with_SIFT_Polyphen$GENE,inTOgen_drivers))]))

VEPs<-VEPs[c('High impact','Moderate impact','Possible impact','Low impact','Unknown impact')]

pdf(paste(resultPath,'/_figures_source/DAM_unreported_VEP_prediction.pdf',sep=''),11,6)
pie(VEPs,col=c("#800026","#fc4e2a","#feb24c","#ffeda0","grey85"))
dev.off()

paste(sum(VEPs[1:2]),'DAMs have a high/moderate functional impact') 
paste(sum(VEPs[1:3]),'DAMs have a high/moderate/possible functional impact') 

load(paste(resultPath,'_all_SAMs.RData',sep=''))

ii<-match(paste(allSAMs$genes,allSAMs$vars),
      paste(allDAMs_with_SIFT_Polyphen$GENE,allDAMs_with_SIFT_Polyphen$var))

allSAMs_with_SIFT_Polyphen<-cbind(allSAMs,allDAMs_with_SIFT_Polyphen[ii,'impact_call'])
colnames(allSAMs_with_SIFT_Polyphen)[ncol(allSAMs_with_SIFT_Polyphen)]<-'impact_call'

save(allSAMs_with_SIFT_Polyphen,file=paste(resultPath,'_allSAMs_with_SIFT_PolyPhen_scores.RData',sep=''))
write.table(allSAMs_with_SIFT_Polyphen,sep="\t",quote=FALSE,file=paste(resultPath,'_allSAMs_with_SIFT_PolyPhen_scores.tsv',sep=''))

samSigs<-paste(allSAMs$cancer_type,allSAMs_with_SIFT_Polyphen$genes,allSAMs_with_SIFT_Polyphen$vars)
usamSigs<-unique(samSigs)
uallSAMs_with_SIFT_Polyphen<-allSAMs_with_SIFT_Polyphen[match(usamSigs,samSigs),]

VEPs<-summary(as.factor(uallSAMs_with_SIFT_Polyphen$impact_call))
VEPs<-VEPs[c('High impact','Moderate impact','Possible impact','Low impact','Unknown impact')]

pdf(paste(resultPath,'/_figures_source/SAM_VEP_prediction.pdf',sep=''),11,6)
pie(VEPs,col=c("#800026","#fc4e2a","#feb24c","#ffeda0","grey85"))
dev.off()



VEPs<-summary(as.factor(uallSAMs_with_SIFT_Polyphen$impact_call[
  which(!is.element(uallSAMs_with_SIFT_Polyphen$genes,inTOgen_drivers))]))
VEPs<-VEPs[c('High impact','Moderate impact','Possible impact','Low impact','Unknown impact')]

VEPs<-VEPs[!is.na(VEPs)]
pdf(paste(resultPath,'/_figures_source/SAM_unreported_VEP_prediction.pdf',sep=''),11,6)
pie(VEPs,col=c("#800026","#fc4e2a","#feb24c"))
dev.off()




