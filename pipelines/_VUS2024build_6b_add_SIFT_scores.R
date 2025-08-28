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


load(paste(resultPath,'_allDAMs.RData',sep=''))

cl_variants <- read_csv(paste(pathdata,'/raw/mutations_all_20241212.csv', sep=""))
### mutations_all_20230202 downloaded from  https://cog.sanger.ac.uk/cmp/download/mutations_all_20241212.zip on 20250221

cl_variants<-cl_variants[which(!is.na(cl_variants$position)),]

variant_signature<-paste(cl_variants$gene_symbol,cl_variants$protein_mutation)
variant_signature1<-paste(allDAMs$GENE,allDAMs$var)

ii<-which(is.element(variant_signature,variant_signature1))
variant_signature<-variant_signature[ii]
cl_variants <- cl_variants[ii,]

uvariant_signature<-unique(variant_signature)

cl_variants<-cl_variants[match(uvariant_signature,variant_signature),]

cl_variants<-cl_variants[order(cl_variants$gene_symbol),]

dt <- as.data.table(cl_variants)  # replace with your object

## 1) seq_region_name
dt[, seq_region_name := str_replace(chromosome, "^chr", "")]

## 2) start / end (SNVs). For indels, adjust end if you wish.
dt[, `:=`(start = as.integer(position),
          end   = as.integer(position))]

## 3) allele_string from REF/ALT
dt[, allele_string := ifelse(!is.na(reference) & !is.na(alternative),
                             paste0(reference, "/", alternative), NA_character_)]

## 4) strand via Ensembl (GRCh38 by default)
ens_ids <- unique(na.omit(ifelse(!is.na(dt$ensembl_gene_id), dt$ensembl_gene_id, dt$gene_id)))

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
strand_map <- getBM(
  attributes = c("ensembl_gene_id","strand"),
  filters    = "ensembl_gene_id",
  values     = ens_ids,
  mart       = mart
)
strand_map <- as.data.table(strand_map)
dt <- strand_map[dt, on = c("ensembl_gene_id" = "ensembl_gene_id")]
setnames(dt, "strand", "strand_numeric")
dt[, strand := as.integer(strand_numeric)][, strand_numeric := NULL]

## 5) final table
out <- dt[, .(gene_symbol, model_id, seq_region_name, start, end, strand, allele_string)]
head(out)

dt<-dt[,-c(4,5,13,14)]

vep_sift_polyphen_safe <- function(df,
                                   batch_size  = 50,
                                   max_retries = 3,
                                   sleep_sec   = 1,
                                   verbose     = TRUE) {
  # ---- Input checks & sanitisation -----------------------------------------
  req_cols <- c("seq_region_name","start","end","reference","alternative")
  stopifnot(all(req_cols %in% names(df)))
  
  x <- as.data.table(df)
  x <- x[!is.na(seq_region_name) & !is.na(start) & !is.na(end) &
           !is.na(reference) & !is.na(alternative)]
  x[, `:=`(start = as.integer(start), end = as.integer(end))]
  x <- x[start > 0 & end >= start]
  
  if (!nrow(x)) {
    warning("No valid variants after basic filtering.")
    return(list(results = data.table(), failed = data.table()))
  }
  
  mk_variant <- function(chr, st, en, ref, alt) sprintf("%s %s %s %s/%s", chr, st, en, ref, alt)
  x[, vep_input := mk_variant(seq_region_name, start, end, reference, alternative)]
  
  batches   <- split(x, ceiling(seq_len(nrow(x))/batch_size))
  n_batches <- length(batches)
  
  if (verbose) {
    message(sprintf("Submitting %d variants in %d batch(es) to Ensembl VEP REST…",
                    nrow(x), n_batches))
    pb <- utils::txtProgressBar(min = 0, max = n_batches, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  
  # ---- Helpers --------------------------------------------------------------
  parse_one <- function(rec) {
    # Return NULL if no transcript consequences
    if (is.null(rec$transcript_consequences)) return(NULL)
    tc <- rec$transcript_consequences
    dt <- rbindlist(lapply(tc, as.data.frame), fill = TRUE)
    
    # Carry over input + core locus fields
    dt[, `:=`(input          = rec$input,
              seq_region_name= rec$seq_region_name,
              start          = rec$start,
              end            = rec$end,
              allele_string  = rec$allele_string)]
    dt
  }
  
  request_batch <- function(dt_batch, depth = 0, batch_id = NA_integer_) {
    if (nrow(dt_batch) == 0) return(list(ok = NULL, bad = NULL))
    
    variants <- as.list(dt_batch$vep_input)
    attempt  <- 0
    
    repeat {
      attempt <- attempt + 1
      if (verbose) {
        msg <- sprintf("Batch %s | size=%d | attempt=%d%s",
                       ifelse(is.na(batch_id), "?", batch_id),
                       nrow(dt_batch), attempt,
                       if (depth > 0) sprintf(" | depth=%d (bisect)", depth) else "")
        message(msg)
      }
      
      res <- try({
        POST(
          url = "https://rest.ensembl.org/vep/human/region",
          add_headers(`Content-Type` = "application/json", `Accept` = "application/json"),
          body = toJSON(list(variants = variants), auto_unbox = TRUE),
          timeout(60)
        )
      }, silent = TRUE)
      
      # Connection error
      if (inherits(res, "try-error")) {
        if (verbose) message("  -> HTTP request failed (connection).")
        if (attempt <= max_retries) { Sys.sleep(sleep_sec * attempt); next } else break
      }
      
      code <- status_code(res)
      
      if (code == 200) {
        ans <- content(res, as = "parsed", type = "application/json")
        ok  <- rbindlist(lapply(ans, parse_one), fill = TRUE)
        return(list(ok = ok, bad = NULL))
      }
      
      if (code %in% c(429, 503)) {
        if (verbose) message(sprintf("  -> Server says %d; backing off…", code))
        if (attempt <= max_retries) { Sys.sleep(sleep_sec * attempt); next } else break
      }
      
      if (code >= 500) {
        if (verbose) message(sprintf("  -> HTTP %d; bisecting batch to isolate offending records…", code))
        if (nrow(dt_batch) == 1 || depth > 10) {
          if (verbose) message("  -> Reached single record or max depth; marking as failed.")
          return(list(ok = NULL, bad = dt_batch))
        }
        mid   <- nrow(dt_batch) %/% 2
        left  <- request_batch(dt_batch[1:mid], depth + 1, batch_id)
        right <- request_batch(dt_batch[(mid+1):nrow(dt_batch)], depth + 1, batch_id)
        return(list(ok = rbindlist(list(left$ok, right$ok), fill = TRUE),
                    bad = rbindlist(list(left$bad, right$bad), fill = TRUE)))
      }
      
      # Other client error (400 etc.) — mark whole batch as failed
      if (verbose) message(sprintf("  -> HTTP %d (client error); marking batch as failed.", code))
      return(list(ok = NULL, bad = dt_batch))
    }
    
    # Exit without success
    list(ok = NULL, bad = dt_batch)
  }
  
  # ---- Process batches ------------------------------------------------------
  ok_list <- list(); bad_list <- list()
  for (i in seq_along(batches)) {
    ans <- request_batch(batches[[i]], batch_id = i)
    if (!is.null(ans$ok))  ok_list[[length(ok_list)+1]]  <- ans$ok
    if (!is.null(ans$bad)) bad_list[[length(bad_list)+1]] <- ans$bad
    if (verbose) utils::setTxtProgressBar(pb, i)
  }
  
  if (verbose) message("\nParsing and summarising transcript-level predictions…")
  resdt <- rbindlist(ok_list, fill = TRUE)
  
  # ---- Robust field harmonisation & summarisation ---------------------------
  if (nrow(resdt)) {
    # Harmonise gene symbol naming
    if (!"gene_symbol" %in% names(resdt) && "hgnc_symbol" %in% names(resdt)) {
      resdt[, gene_symbol := hgnc_symbol]
    }
    # Ensure consequence_terms exists
    if (!"consequence_terms" %in% names(resdt)) {
      if ("consequence_term" %in% names(resdt)) {
        resdt[, consequence_terms := consequence_term]
      } else if ("consequence" %in% names(resdt)) {
        resdt[, consequence_terms := consequence]
      } else {
        resdt[, consequence_terms := NA_character_]
      }
    }
    # Ensure prediction/score columns exist
    for (cn in c("sift_score","sift_prediction","polyphen_score","polyphen_prediction",
                 "gene_symbol","transcript_id")) {
      if (!cn %in% names(resdt)) resdt[[cn]] <- NA
    }
    
    # Build join key and merge back to inputs
    resdt[, join_key := sprintf("%s %s %s %s", seq_region_name, start, end, allele_string)]
    x[, join_key := sprintf("%s %s %s %s/%s", seq_region_name, start, end, reference, alternative)]
    
    keep <- resdt[, .(join_key, gene_symbol, transcript_id, consequence_terms,
                      sift_score, sift_prediction, polyphen_score, polyphen_prediction)]
    
    suppressWarnings({
      keep[, sift_score := as.numeric(sift_score)]
      keep[, polyphen_score := as.numeric(polyphen_score)]
    })
    
    summary <- keep[
      , .(sift_min          = if (all(is.na(sift_score))) NA_real_ else min(sift_score, na.rm = TRUE),
          sift_pred_any     = paste(na.omit(unique(sift_prediction)), collapse = "|"),
          polyphen_max      = if (all(is.na(polyphen_score))) NA_real_ else max(polyphen_score, na.rm = TRUE),
          polyphen_pred_any = paste(na.omit(unique(polyphen_prediction)), collapse = "|")),
      by = join_key
    ]
    
    annotated <- merge(x, summary, by = "join_key", all.x = TRUE)
  } else {
    annotated <- x
    annotated[, c("sift_min","sift_pred_any","polyphen_max","polyphen_pred_any") := NA]
  }
  
  failed <- rbindlist(bad_list, fill = TRUE)
  
  if (verbose) {
    message(sprintf("Done. Annotated: %d | Failed: %d", nrow(annotated), nrow(failed)))
    if (nrow(failed)) {
      message("A few examples of failed inputs:")
      print(head(failed[, .(seq_region_name, start, end, reference, alternative)], 5))
      message("Common causes: wrong build, REF/ALT mismatch, large/complex indels.")
    }
  }
  
  list(
    results = annotated[, !c("join_key","vep_input")],
    failed  = failed
  )
}


annotated <- vep_sift_polyphen_safe(dt)

allDAMs_sigs<-paste(allDAMs$GENE,allDAMs$var)

eff_sigs<-paste(annotated$results$gene_symbol,annotated$results$protein_mutation)

