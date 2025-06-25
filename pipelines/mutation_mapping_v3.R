################################################################################
### Load libraries
################################################################################
library(ensemblVEP)
library(dplyr)
library(reshape2)
library(ggplot2)
library(readr)
library(BiocIO)
library(httr)
library(jsonlite)


################################################################################
### Create annotation function using VEP REST service
################################################################################
f_get_vep_annotation <- function(
    fs_genome_ver,
    fs_vep_id,
    id_type,
    fd_variant_details
) {
  
  
  ################################################################################
  ### Create VEP REST call bases on genome version
  ################################################################################
  fs_vep_input <- {paste0(
    {switch(
      fs_genome_ver,
      GRCh37 = "https://grch37.rest.ensembl.org/vep/human/hgvs/",
      GRCh38 = "https://rest.ensembl.org/vep/human/hgvs/",
    )},
    fs_vep_id,
    "?content-type=application/json"
  )}
  
  
  ################################################################################
  ### Get annotation from VEP
  ################################################################################
  fo_vep_output <- try(fromJSON(rawToChar(GET(fs_vep_input)$content)))
  
  
  ################################################################################
  ### Create result dataframe
  ################################################################################
  if ({
    class(fo_vep_output) != "try-error" &
      {all(
        c(
          "id", "input", "assembly_name",
          "seq_region_name", "start", "end", "strand", "allele_string"
        ) %in%
        names(fo_vep_output)
      )}
  }) {
    fd_vep_output <- {left_join(
      fd_variant_details,
      {data.frame(
        id = fo_vep_output$id,
        genome_ver = fo_vep_output$assembly_name,
        vep_input = fo_vep_output$input,
        seq_region_name = fo_vep_output$seq_region_name,
        start = fo_vep_output$start,
        end = fo_vep_output$end,
        strand = fo_vep_output$strand,
        allele_string = fo_vep_output$allele_string
      ) %>%
          rename(!!id_type := id)
      }
    )}
  } else {
    fd_vep_output <- data.frame()
  }
  
  
  ################################################################################
  ### Add consequences if any match is found
  ################################################################################
  if ({
    !is.null(fd_vep_output) &
      class(fo_vep_output) != "try-error" &
      "transcript_consequences" %in% names(fo_vep_output)
  }) {
    try({fd_vep_output <- left_join(
      fd_vep_output,
      {fo_vep_output$transcript_consequences %>%
          bind_rows() %>%
          filter(!is.na(cds_start) & !is.na(protein_start)) %>%
          mutate(
            cdna_mutation = paste0(
              "c.", 
              ifelse(
                cds_start == cds_end,
                cds_start,
                paste(cds_start, cds_end, sep = "_")
              ),
              substring("tTc/tGc", 2, 2), ">", substring("tTc/tGc", 6, 6)
            )
          )
      }
    )})
  }
  
  
  return(fd_vep_output)
}


################################################################################
### Deine input / output files and folder
################################################################################
# hg19/GRCh37: https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/
# hg38/GRCh38: https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/
# All mutations: https://cog.sanger.ac.uk/cmp/download/mutations_all_20241212.zip
file_input_genomes_list <- {list(
  GRCh37 = "C:\\Data\\VUS\\Homo_sapiens.GRCh37.75.gtf.gz",
  GRCh38 = "C:\\Data\\VUS\\Homo_sapiens.GRCh38.113.gtf.gz"
)}
file_input_variasts <- "C:\\Data\\VUS\\_allDAMs.tsv"
file_input_allmutations <- "C:\\Data\\VUS\\mutations_all_20241212.zip"
file_output_vepinput <- "C:\\Data\\VUS\\vepinput"
setwd("C:\\Data\\VUS")


################################################################################
### Load input files
################################################################################
res_genomes <- {as.list(sapply(
  file_input_genomes_list,
  function(fs) {
    message(paste("Import", fs, Sys.time()))
    return(BiocIO::import(fs) %>% as.data.frame())
  }
))}
res_allmutations <- as.data.frame({read_csv(
  file_input_allmutations,
  col_names = TRUE
)})
res_variants_source <- as.data.frame({read_tsv(
  file_input_variasts,
  col_names = TRUE
)})
gc();gc();gc();gc()


################################################################################
### Add mutations details to variants
################################################################################
res_variants <- {left_join(
  res_variants_source,
  {res_allmutations %>%
    rename(
      sanger_gene_id = gene_id,
      gene_name = gene_symbol,
      gene_id = ensembl_gene_id,
      transcript_CCDS = transcript_id
    )
  }
)}
gc();gc();gc();gc()


################################################################################
### Add Transcripts IDs
################################################################################
res_variants <- {left_join(
  res_variants,
  {lapply(
    names(res_genomes),
    function(fs) {res_genomes[[fs]] %>%
      filter(type == "transcript") %>%
      dplyr::select(any_of(c(
        "gene_name", "gene_id",
        "transcript_id", "transcript_version"
      ))) %>%
      mutate(genome_ver = fs, .before = 1)
    }
  ) %>% bind_rows()}
)}
gc();gc();gc();gc()


################################################################################
### Create VEP id for the mapping
################################################################################
res_variants <- {res_variants %>%
  mutate(
    id_vep_cdna = ifelse(
      genome_ver == "GRCh38",
      paste0(transcript_id, ".", transcript_version, ":", cdna_mutation),
      paste0(transcript_id, ":", cdna_mutation)
    ),
    id_vep_rna = ifelse(
      genome_ver == "GRCh38",
      paste0(transcript_id, ".", transcript_version, ":", rna_mutation),
      paste0(transcript_id, ":", rna_mutation)
    ),
    id_vep_protein = ifelse(
      genome_ver == "GRCh38",
      paste0(transcript_id, ".", transcript_version, ":", protein_mutation),
      paste0(transcript_id, ":", protein_mutation)
    ),
  )
}
gc();gc();gc();gc()


################################################################################
### Annotate Protein mutations using VEP
################################################################################
if (file.exists("vep_annotations.rds")) {
  res_vep <- readRDS("vep_annotations.rds")
} else {
  res_vep <- list()
}
for (current_genome_ver in sort(names(res_genomes), decreasing = TRUE)) {
  if(current_genome_ver %in% names(res_vep)) res_vep[[current_genome_ver]] <- {
    data.frame(genome_ver = character(), id_vep_protein = character())
  }
  for (current_vep_var in {res_variants %>%
    filter(
      genome_ver == current_genome_ver &
      !id_vep_protein %in% res_vep[[current_genome_ver]]$id_vep_protein
    ) %>%
    dplyr::select(id_vep_protein) %>%
    unique() %>%
    unlist()
  }) {
    curent_vep_output <- data.frame()
    
    
    ################################################################################
    ### Run VEP Annotation
    ################################################################################
    for (curent_vep_annotationtype in c(
      "id_vep_protein", "id_vep_rna", "id_vep_cdna" 
    )) {
      if (nrow(curent_vep_output) == 0) {
        print(paste("Annotating ", curent_vep_annotationtype,":", current_genome_ver, current_vep_var, Sys.time()))
        curent_vep_output <- {f_get_vep_annotation(
          fs_genome_ver = current_genome_ver, 
          fs_vep_id = {res_variants %>%
            filter(id_vep_protein == current_vep_var) %>%
            select(all_of(curent_vep_annotationtype)) %>%
            unique()
          },
          id_type = curent_vep_annotationtype,
          fd_variant_details = {res_variants %>%
              filter(
                genome_ver == current_genome_ver &
                  id_vep_protein == current_vep_var
              ) %>%
              dplyr::select(
                genome_ver, gene_name, gene_id, transcript_id,
                chromosome, position,
                cdna_mutation, rna_mutation, protein_mutation, effect,
                id_vep_cdna, id_vep_rna, id_vep_protein
              ) %>%
              unique()
          }
        )}
      }
    }
    
    
    ################################################################################
    ### Merge results
    ################################################################################
    if (nrow(curent_vep_output) > 0) {
      res_vep[[current_genome_ver]] <- {bind_rows(
        res_vep[[current_genome_ver]],
        curent_vep_output
      )}
    }
  }
  saveRDS(res_vep, "vep_annotations.rds")
}


################################################################################
### Report variants
################################################################################
for (current_genome_ver in sort(names(res_genomes), decreasing = TRUE)) {
  {left_join(
    {left_join(
      res_variants_source,
      {res_variants %>%
        dplyr::select(
          gene_name, gene_id,
          chromosome, position,
          cdna_mutation, rna_mutation, protein_mutation, effect
        ) %>%
        unique()
      }
    )},
    {res_vep[[current_genome_ver]] %>%
        filter(
          current_genome_ver == "GRCh37" |
          (
            chromosome == paste0("chr", seq_region_name) &
            position >= start - 1 &
            position <= end + 1
          )
        ) %>%
        dplyr::select(
          gene_name, gene_id,
          chromosome, position,
          cdna_mutation, rna_mutation, protein_mutation, effect,
          seq_region_name, start, end, strand, allele_string 
        ) %>%
        unique()
    }
  ) %>%
    write.table(
      file = paste0("VUS_location_", current_genome_ver,".txt"),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = TRUE
    )
  }
}
