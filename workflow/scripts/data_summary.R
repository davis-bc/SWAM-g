#!/usr/bin/env Rscript

libraries <- c("tidyverse", "writexl")

invisible(lapply(libraries, function(x) {
  suppressMessages(suppressWarnings(library(x, character.only = T)))
  }))


# Input files

checkm2_files        <- snakemake@input[["checkm"]]
mash_file            <- snakemake@input[["mash"]]
mlst_file            <- snakemake@input[["mlst"]]
resfinder_files      <- snakemake@input[["resfinder_files"]]
pf_files             <- snakemake@input[["pf_files"]]
afp_files            <- snakemake@input[["afp_files"]]
seqsero_files        <- snakemake@input[["seqsero_files"]]
sistr_files          <- snakemake@input[["sistr_files"]]
ectyper_files        <- snakemake@input[["ectyper_files"]]
mobtyper_files       <- snakemake@input[["mobtyper_files"]]
coverage_files       <- snakemake@input[["coverage_files"]]
mef_files            <- snakemake@input[["mef_files"]]
contig_files         <- snakemake@input[["contig_files"]]
mobtyper_blast_files <- snakemake@input[["mobtyper_blast_files"]]
txsscan_files        <- snakemake@input[["txsscan_files"]]
txsscan_files2       <- snakemake@input[["txsscan_files2"]]
prodigal_files       <- snakemake@input[["prodigal_files"]]

# Output files

out_file1 <- snakemake@output[[1]]
out_file2 <- snakemake@output[[2]]

# Debug mode
debug_mode <- isTRUE(snakemake@params[["debug"]]) ||
              identical(snakemake@params[["debug"]], "true") ||
              identical(snakemake@params[["debug"]], "True")

if (debug_mode) {
  options(warn = 1)   # Print warnings immediately as they occur
  message("[DEBUG] data_summary.R: debug mode enabled")
  message("[DEBUG] Samples: ", length(checkm2_files))
}

dlog <- function(...) {
  if (debug_mode) message("[DEBUG] ", ...)
}

normalize_empty <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x[x %in% c("", "-", "NA", "N/A", "NULL", "null", "Not available")] <- NA_character_
  x
}

display_serotype <- function(x) {
  x <- normalize_empty(x)
  ifelse(
    is.na(x),
    NA_character_,
    x %>%
      str_replace("^Salmonella\\s+enterica\\s+(subsp\\.?\\s+enterica\\s+)?", "") %>%
      str_replace("^serovar\\s+", "")
  )
}

canonical_serotype <- function(x) {
  x <- display_serotype(x)
  ifelse(
    is.na(x),
    NA_character_,
    x %>%
      str_to_lower() %>%
      str_replace_all("[[:punct:]]", " ") %>%
      str_squish()
  )
}

is_ambiguous_serotype <- function(x) {
  x <- display_serotype(x)
  ifelse(
    is.na(x),
    TRUE,
    str_detect(
      str_to_lower(x),
      "ambiguous|multiple|unknown|unidentified|^\\-+(:\\-+)*$|\\bor\\b|and/or|not\\s+predicted"
    )
  )
}

score_seqsero <- function(identification, serotype, antigenic_profile, status) {
  status <- str_to_upper(coalesce(normalize_empty(status), "REPORTED"))
  identification <- str_to_lower(coalesce(normalize_empty(identification), ""))
  serotype <- display_serotype(serotype)
  antigenic_profile <- normalize_empty(antigenic_profile)

  if (status == "SKIPPED_NOT_SALMONELLA") {
    return(-5)
  }

  score <- 0
  if (str_detect(identification, "salmonella")) score <- score + 1
  if (!is.na(serotype) && !isTRUE(is_ambiguous_serotype(serotype))) {
    score <- score + 3
  } else if (!is.na(serotype)) {
    score <- score + 1
  }
  if (!is.na(antigenic_profile)) score <- score + 1
  if (str_detect(identification, "contamin")) score <- score - 2
  if (str_detect(identification, "unident|incorrect")) score <- score - 2
  score
}

score_sistr <- function(serovar, serovar_antigen, serovar_cgmlst, qc_status, qc_messages, cgmlst_st) {
  qc_status <- str_to_upper(coalesce(normalize_empty(qc_status), ""))
  serovar <- display_serotype(serovar)
  serovar_antigen <- display_serotype(serovar_antigen)
  serovar_cgmlst <- display_serotype(serovar_cgmlst)
  qc_messages <- normalize_empty(qc_messages)
  cgmlst_st <- normalize_empty(cgmlst_st)

  if (qc_status == "SKIPPED_NOT_SALMONELLA") {
    return(-5)
  }

  score <- 0
  if (!is.na(serovar) && !isTRUE(is_ambiguous_serotype(serovar))) {
    score <- score + 3
  } else if (!is.na(serovar)) {
    score <- score + 1
  }
  if (!is.na(serovar_antigen)) score <- score + 1
  if (!is.na(serovar_cgmlst)) score <- score + 1
  if (!is.na(cgmlst_st)) score <- score + 1
  if (qc_status == "PASS") score <- score + 2
  if (qc_status == "FAIL") score <- score - 2
  if (!is.na(qc_messages)) score <- score - 1
  score
}

pick_salmonella_consensus <- function(seqsero_identification, seqsero_serotype, seqsero_antigenic_profile,
                                      seqsero_status, sistr_serovar, sistr_serovar_antigen,
                                      sistr_serovar_cgmlst, sistr_qc_status, sistr_qc_messages, sistr_cgmlst_st) {
  seq_key <- canonical_serotype(seqsero_serotype)
  sistr_key <- canonical_serotype(sistr_serovar)
  seq_score <- score_seqsero(seqsero_identification, seqsero_serotype, seqsero_antigenic_profile, seqsero_status)
  sistr_score <- score_sistr(sistr_serovar, sistr_serovar_antigen, sistr_serovar_cgmlst, sistr_qc_status, sistr_qc_messages, sistr_cgmlst_st)

  agreement <- case_when(
    !is.na(seq_key) && !is.na(sistr_key) && seq_key == sistr_key ~ "agree",
    !is.na(seq_key) && !is.na(sistr_key) ~ "disagree",
    TRUE ~ NA_character_
  )

  if (is.na(seq_key) && is.na(sistr_key)) {
    return(tibble(
      Predicted_Serotype = NA_character_,
      Salmonella_Serotype_Source = NA_character_,
      Serotype_Agreement = agreement,
      SeqSero2_Score = seq_score,
      SISTR_Score = sistr_score
    ))
  }

  if (!is.na(seq_key) && !is.na(sistr_key) && seq_key == sistr_key) {
    return(tibble(
      Predicted_Serotype = coalesce(display_serotype(sistr_serovar), display_serotype(seqsero_serotype)),
      Salmonella_Serotype_Source = "SeqSero2+SISTR",
      Serotype_Agreement = agreement,
      SeqSero2_Score = seq_score,
      SISTR_Score = sistr_score
    ))
  }

  prefer_sistr <- FALSE
  if (is.na(seq_key)) {
    prefer_sistr <- TRUE
  } else if (!is.na(sistr_key) && sistr_score > seq_score) {
    prefer_sistr <- TRUE
  } else if (!is.na(sistr_key) && sistr_score == seq_score &&
             str_to_upper(coalesce(normalize_empty(sistr_qc_status), "")) == "PASS") {
    prefer_sistr <- TRUE
  }

  source <- if (prefer_sistr) "SISTR" else "SeqSero2"
  if (agreement == "disagree") {
    source <- paste0(source, "_preferred")
  }

  tibble(
    Predicted_Serotype = if (prefer_sistr) display_serotype(sistr_serovar) else display_serotype(seqsero_serotype),
    Salmonella_Serotype_Source = source,
    Serotype_Agreement = agreement,
    SeqSero2_Score = seq_score,
    SISTR_Score = sistr_score
  )
}


###############################################
###      Read in pre-generated tables
###############################################

dlog("START: MASH taxonomy")
gtdbtk <- read_tsv(mash_file, col_types = cols())
gtdbtk$Sample <- gsub("\\.chromosome$", "", gtdbtk$user_genome)
gtdbtk$Species <- sub(".*s__", "", gtdbtk$classification)
dlog("DONE:  MASH taxonomy (", nrow(gtdbtk), " rows)")

dlog("START: MLST")
mlst <- read_table(mlst_file, col_names = F, col_types = cols())
names(mlst)[1:3] <- c("Sample", "Scheme", "Sequence_Type")
mlst$Sample <- sub(".*/(.*)\\.fasta", "\\1", mlst$Sample)
dlog("DONE:  MLST (", nrow(mlst), " rows)")

########################################################
###       Gather checkm2 and coverage data, assess assembly quality
########################################################

dlog("START: CheckM2 + coverage")
checkm2 <- do.call(bind_rows, lapply(checkm2_files, function(f) {
  tryCatch(read_tsv(f, col_types = cols(.default = col_character())), error = function(e) NULL)
})) %>%
  mutate(across(c(Completeness, Contamination, Coding_Density, Contig_N50, Average_Gene_Length,
                  Genome_Size, GC_Content, Total_Coding_Sequences, Total_Contigs,
                  Max_Contig_Length, Translation_Table_Used), as.numeric))

coverage <- do.call(rbind, lapply(coverage_files, function(f) 
  read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)))

colnames(coverage) <- c("Sample", "XCoverage")

assembly_stats <- checkm2 %>% 
                   rename(Sample=Name) %>% 
                   left_join(coverage, by="Sample") %>%
                   mutate(QA = ifelse(Contig_N50 > 20000 & Total_Contigs < 500 & XCoverage >= 30, "pass", "fail"))
dlog("DONE:  CheckM2 + coverage (", nrow(assembly_stats), " samples)")



###############################################
###   Gather AMRFinderPlus results
###############################################

dlog("START: AMRFinderPlus")
amrfinderplus <- do.call(bind_rows, lapply(afp_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols(`Contig id` = "c")), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) {return(NULL)}
  x$Sample <- gsub("\\.afp.tsv$", "", basename(f))
  x
}))

if (!is.null(amrfinderplus) && ncol(amrfinderplus) > 0) {
  afp_genotype <- amrfinderplus %>% filter(Type == "AMR") %>% group_by(Sample) %>%
                  summarise(AMRFinder_Genotype = paste(sort(unique(`Element symbol`)), collapse = ", "))
  dlog("DONE:  AMRFinderPlus (", nrow(amrfinderplus), " hits across ", n_distinct(amrfinderplus$Sample), " samples)")
} else {
  amrfinderplus <- tibble()
  afp_genotype  <- tibble(Sample = character(), AMRFinder_Genotype = character())
  dlog("DONE:  AMRFinderPlus (no hits in any sample)")
}


###############################################
###   Gather all Serotyping results
###############################################

dlog("START: SeqSero2 (Salmonella serotyping)")
seqsero <- do.call(bind_rows, lapply(seqsero_files, function(f) {
  x <- tryCatch(
    read_tsv(f, col_types = cols(.default = col_character())),
    error=function(e) NULL
  )
  if (is.null(x) || nrow(x) == 0) return(NULL)
  x %>%
    mutate(
      Sample = coalesce(
        if ("Sample name" %in% names(x)) x[["Sample name"]] else NA_character_,
        if ("Output directory" %in% names(x)) basename(x[["Output directory"]]) else NA_character_,
        basename(dirname(f))
      ),
      SeqSero2_Identification = if ("Predicted identification" %in% names(x)) x[["Predicted identification"]] else NA_character_,
      SeqSero2_Antigenic_Profile = if ("Predicted antigenic profile" %in% names(x)) x[["Predicted antigenic profile"]] else NA_character_,
      SeqSero2_Serotype = if ("Predicted serotype" %in% names(x)) x[["Predicted serotype"]] else NA_character_,
      SeqSero2_Status = if ("Workflow status" %in% names(x)) x[["Workflow status"]] else "REPORTED"
    ) %>%
    select(Sample, SeqSero2_Identification, SeqSero2_Antigenic_Profile, SeqSero2_Serotype, SeqSero2_Status)
}))

if (!is.null(seqsero) && ncol(seqsero) > 0) {
  seqsero <- seqsero %>%
    mutate(across(c(SeqSero2_Identification, SeqSero2_Antigenic_Profile, SeqSero2_Serotype, SeqSero2_Status), normalize_empty))
} else {
  seqsero <- tibble(
    Sample = character(),
    SeqSero2_Identification = character(),
    SeqSero2_Antigenic_Profile = character(),
    SeqSero2_Serotype = character(),
    SeqSero2_Status = character()
  )
}

seqsero_simple <- seqsero %>%
  filter(coalesce(SeqSero2_Status, "REPORTED") != "SKIPPED_NOT_SALMONELLA") %>%
  relocate(Sample)
dlog("DONE:  SeqSero2 (", nrow(seqsero_simple), " Salmonella samples)")

dlog("START: SISTR (Salmonella serotyping)")
sistr <- do.call(bind_rows, lapply(sistr_files, function(f) {
  x <- tryCatch(
    read_tsv(f, col_types = cols(.default = col_character())),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(x) == 0) return(NULL)
  x %>%
    mutate(
      Sample = coalesce(if ("genome" %in% names(x)) x[["genome"]] else NA_character_, basename(dirname(f))),
      SISTR_Serovar = if ("serovar" %in% names(x)) x[["serovar"]] else NA_character_,
      SISTR_Serovar_Antigen = if ("serovar_antigen" %in% names(x)) x[["serovar_antigen"]] else NA_character_,
      SISTR_Serovar_CGMLST = if ("serovar_cgmlst" %in% names(x)) x[["serovar_cgmlst"]] else NA_character_,
      SISTR_cgMLST_ST = if ("cgmlst_ST" %in% names(x)) x[["cgmlst_ST"]] else NA_character_,
      SISTR_Serogroup = if ("serogroup" %in% names(x)) x[["serogroup"]] else NA_character_,
      SISTR_O_Antigen = if ("o_antigen" %in% names(x)) x[["o_antigen"]] else NA_character_,
      SISTR_H1 = if ("h1" %in% names(x)) x[["h1"]] else NA_character_,
      SISTR_H2 = if ("h2" %in% names(x)) x[["h2"]] else NA_character_,
      SISTR_QC_Status = if ("qc_status" %in% names(x)) x[["qc_status"]] else NA_character_,
      SISTR_QC_Messages = if ("qc_messages" %in% names(x)) x[["qc_messages"]] else NA_character_
    ) %>%
    select(Sample, SISTR_Serovar, SISTR_Serovar_Antigen, SISTR_Serovar_CGMLST,
           SISTR_cgMLST_ST, SISTR_Serogroup, SISTR_O_Antigen, SISTR_H1, SISTR_H2,
           SISTR_QC_Status, SISTR_QC_Messages)
}))

if (!is.null(sistr) && ncol(sistr) > 0) {
  sistr <- sistr %>%
    mutate(across(c(SISTR_Serovar, SISTR_Serovar_Antigen, SISTR_Serovar_CGMLST, SISTR_cgMLST_ST,
                    SISTR_Serogroup, SISTR_O_Antigen, SISTR_H1, SISTR_H2, SISTR_QC_Status,
                    SISTR_QC_Messages), normalize_empty))
} else {
  sistr <- tibble(
    Sample = character(),
    SISTR_Serovar = character(),
    SISTR_Serovar_Antigen = character(),
    SISTR_Serovar_CGMLST = character(),
    SISTR_cgMLST_ST = character(),
    SISTR_Serogroup = character(),
    SISTR_O_Antigen = character(),
    SISTR_H1 = character(),
    SISTR_H2 = character(),
    SISTR_QC_Status = character(),
    SISTR_QC_Messages = character()
  )
}

sistr_simple <- sistr %>%
  filter(coalesce(SISTR_QC_Status, "REPORTED") != "SKIPPED_NOT_SALMONELLA") %>%
  relocate(Sample)
dlog("DONE:  SISTR (", nrow(sistr_simple), " Salmonella samples)")

salmonella_serotype <- gtdbtk %>%
  select(Sample, Species) %>%
  left_join(seqsero, by = "Sample") %>%
  left_join(sistr, by = "Sample") %>%
  mutate(
    across(c(SeqSero2_Identification, SeqSero2_Antigenic_Profile, SeqSero2_Serotype, SeqSero2_Status,
             SISTR_Serovar, SISTR_Serovar_Antigen, SISTR_Serovar_CGMLST, SISTR_cgMLST_ST,
             SISTR_Serogroup, SISTR_O_Antigen, SISTR_H1, SISTR_H2, SISTR_QC_Status, SISTR_QC_Messages),
           normalize_empty)
  ) %>%
  {
    bind_cols(
      .,
      pmap_dfr(
        list(.$SeqSero2_Identification, .$SeqSero2_Serotype, .$SeqSero2_Antigenic_Profile, .$SeqSero2_Status,
             .$SISTR_Serovar, .$SISTR_Serovar_Antigen, .$SISTR_Serovar_CGMLST, .$SISTR_QC_Status,
             .$SISTR_QC_Messages, .$SISTR_cgMLST_ST),
        pick_salmonella_consensus
      )
    )
  } %>%
  mutate(
    SeqSero2_Serotype = display_serotype(SeqSero2_Serotype),
    SISTR_Serovar = display_serotype(SISTR_Serovar),
    SISTR_Serovar_Antigen = display_serotype(SISTR_Serovar_Antigen),
    SISTR_Serovar_CGMLST = display_serotype(SISTR_Serovar_CGMLST)
  ) %>%
  filter(
    Species == "Salmonella enterica" |
      !is.na(Predicted_Serotype)
  )

salmonella_serotype_simple <- salmonella_serotype %>%
  relocate(Sample, Species, Predicted_Serotype, Salmonella_Serotype_Source, Serotype_Agreement)

dlog("START: ECTyper (E. coli serotyping)")
ectyper <- do.call(bind_rows, lapply(ectyper_files, function(f) {
  x <- tryCatch(
    read_tsv(f, col_types = cols(.default = col_character())),
    error=function(e) NULL
  )
  if (is.null(x) || nrow(x) == 0) return(NULL)
  x$Sample <- basename(dirname(f))
  x
}))

ectyper_simple <- if (!is.null(ectyper) && ncol(ectyper) > 0) {
  ectyper %>% filter(grepl("Escherichia", Species)) %>% relocate(Sample) %>% select(-Name)
} else {
  tibble()
}
dlog("DONE:  ECTyper (", nrow(ectyper_simple), " E. coli samples)")

###############################################
###       Gather MOB-suite plasmid data
###############################################

dlog("START: MOB-suite mobtyper")
mobtyper <- do.call(bind_rows, lapply(mobtyper_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols(.default = col_character())), error=function(e) NULL)
  if (is.null(x) || nrow(x) == 0) return(NULL)
  x
}))

if (!is.null(mobtyper) && ncol(mobtyper) > 0) {
  mobtyper <- mobtyper %>%
    rename(rep_type        = `rep_type(s)`,
           relaxase_type   = `relaxase_type(s)`,
           orit_type       = `orit_type(s)`,
           associated_pmid = `associated_pmid(s)`) %>%
    separate(sample_id, into = c("Sample", "primary_cluster_id"), sep = ":", extra = "drop") %>%
    mutate(PlasmidID = paste0(Sample, ".plasmid_", primary_cluster_id)) %>%
    relocate(Sample, PlasmidID)

  mobtyper_summary <- mobtyper %>%
                      group_by(Sample) %>%
                      summarise(n_Plasmid_Clusters      = length(primary_cluster_id),
                                Plasmid_Rep_Types        = paste(rep_type,      collapse = " | "),
                                Plasmid_Relaxase_Types   = paste(relaxase_type, collapse = " | "))
  dlog("DONE:  MOB-suite mobtyper (", nrow(mobtyper), " plasmids across ", n_distinct(mobtyper$Sample), " samples)")
} else {
  mobtyper         <- tibble()
  mobtyper_summary <- tibble(Sample = character(), n_Plasmid_Clusters = integer(),
                             Plasmid_Rep_Types = character(), Plasmid_Relaxase_Types = character())
  dlog("DONE:  MOB-suite mobtyper (no plasmids detected)")
}



dlog("START: MOB-suite biomarkers (oriT/oriV)")
mobtyper_blast_raw <- do.call(bind_rows, lapply(mobtyper_blast_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols(.default = col_character())), error=function(e) NULL)
  if (is.null(x) || nrow(x) == 0) return(NULL)
  x$Sample <- basename(dirname(f))
  x
}))

if (!is.null(mobtyper_blast_raw) && ncol(mobtyper_blast_raw) > 0) {
  mobtyper_blast <- mobtyper_blast_raw %>%
    separate(qseqid, into = c("accession", "gene"), sep="\\|", extra = "drop") %>%
    mutate(contig.num = str_split_i(sseqid, " ", 1), sstrand=ifelse(sstrand=="minus", "negative", "positive")) %>%
    select(Sample, contig.num, gene, biomarker, sstart, send, sstrand) %>%
    rename(type=biomarker, start=sstart, end=send, strand=sstrand) %>%
    mutate(start = as.integer(start), end = as.integer(end))
  dlog("DONE:  MOB-suite biomarkers (", nrow(mobtyper_blast), " hits)")
} else {
  mobtyper_blast <- tibble(Sample=character(), contig.num=character(), gene=character(),
                            type=character(), start=integer(), end=integer(), strand=character())
  dlog("DONE:  MOB-suite biomarkers (no hits)")
}

###############################################
###      Parse Resfinder phenotypes
###############################################

dlog("START: ResFinder")
safe_read_resfinder <- function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols(.default = col_character())), error=function(e) NULL)
  sample_name <- basename(dirname(f))
  
  # If empty, or only header, return row with Sample and NA Phenotype
  if (is.null(x) || nrow(x) == 0) {
    return(tibble(Sample = sample_name, Phenotype = NA))
  }
  
  # Coerce Identity if present
  if ("Identity" %in% names(x)) x$Identity <- as.character(x$Identity)
  # Ensure Sample column
  x$Sample <- sample_name
  # If Phenotype missing, add it
  if (!("Phenotype" %in% names(x))) x$Phenotype <- NA
  return(x)
}

resfinder <- do.call(bind_rows, lapply(resfinder_files, safe_read_resfinder))
dlog("DONE:  ResFinder (", nrow(resfinder), " rows)")

dlog("START: ResFinder phenotype + genotype summary")

phenotype <- resfinder %>% group_by(Sample) %>%
             summarise(phenotype = list(unlist(str_split(Phenotype, ",\\s*")))) %>% 
             mutate(Resfinder_Predicted_Phenotype = map_chr(phenotype, ~ paste(sort(unique(trimws(.))), collapse = ", "))) %>%
             select(-phenotype)
             
resfinder_genotype <- resfinder %>% group_by(Sample) %>% 
  summarise(Resfinder_Genotype = paste(sort(unique(`Resistance gene`)), collapse = ", "))
dlog("DONE:  ResFinder phenotype + genotype summary")      
  

###############################################
###        PointFinder + DisinFinder
###############################################

dlog("START: PointFinder")
pf <- do.call(
  bind_rows,
  lapply(pf_files, function(f) {
    x <- tryCatch(read_tsv(f, col_types = cols(.default = col_character())), error = function(e) NULL)
    if (is.null(x) || nrow(x) == 0) return(NULL)
    x$Sample <- basename(dirname(f))
    # Coerce `PMID` (or other problematic columns) to character for consistency
    if ("PMID" %in% colnames(x)) {
      x$PMID <- as.character(x$PMID)
    }
    x
  })
)

if (!is.null(pf) && ncol(pf) > 0) {
  pf <- pf %>%
    group_by(Sample) %>%
    summarise(Pointfinder_Mutation = paste(Mutation, collapse = " | "))
  dlog("DONE:  PointFinder (", nrow(pf), " samples with mutations)")
} else {
  pf <- tibble(Sample = character(), Pointfinder_Mutation = character())
  dlog("DONE:  PointFinder (no mutations)")
}


# df_files <- list.files("resfinder", pattern = "DisinFinder_results_tab.txt", recursive = 1, full.names = T)
# 
# df <- do.call(bind_rows, lapply(df_files, function(f) {
#   x <- tryCatch(read_tsv(f, col_types = cols()), error=function(e) NULL)
#   if (is.null(x)) return(NULL)
#   x$Sample <- basename(dirname(f))
#   x
# }))

###############################################
###           TXSScan results
###############################################

dlog("START: TXSScan + Prodigal coordinate mapping")
process_txsscan_and_prodigal_files <- function(txsscan_files, prodigal_files) {
  # Helper function to parse TXSScan results
  parse_all_systems_file <- function(filepath) {
    data <- read_tsv(filepath, comment = "#", col_types = cols(.default = "c"))
    lines <- readLines(filepath)
    
    # Extract system descriptions from replicon blocks
    system_descriptions <- tibble(line = lines) %>%
      filter(str_starts(line, "#")) %>% 
      mutate(block = cumsum(str_detect(lead(line, default = ""), "^replicon"))) %>% 
      filter(!is.na(block), block > 0) %>% 
      group_by(block) %>% 
      summarise(
        system_name = str_remove(last(str_trim(line)), "^# "),
        .groups = "drop"
      )
    
    # Annotate data with system names and sample info
    data <- data %>%
      mutate(block = cumsum(str_detect(replicon, "^replicon"))) %>%
      left_join(system_descriptions, by = "block") %>%
      select(-block) %>%
      filter(!replicon == "replicon") %>%
      mutate(Sample = sub(".prot", "", replicon)) %>%
      relocate(Sample)
    
    return(data)
  }
  
  # Helper function to parse prodigal headers
  parse_prodigal_headers <- function(filepath) {
    headers <- read_lines(filepath) %>% 
      .[str_starts(., ">")] %>% 
      map_df(~{
        # Parse Prodigal headers string
        info <- str_split(.x, " # ")[[1]]
        orf_id <- sub("^>", "", info[1])  # Extract the gene identifier
        contig_start <- as.integer(info[2])  # Extract the start position
        contig_end <- as.integer(info[3])    # Extract the stop position
        strand <- ifelse(info[4] == "1", "positive", "negative")  # Map orientation
        tibble(orf_id = orf_id, contig_start = contig_start, contig_end = contig_end, strand = strand)
      })
    if (any(is.na(headers$contig_start) | is.na(headers$contig_end))) {
      stop(sprintf("Prodigal file '%s' has improperly parsed numeric values.", filepath))
    }
    return(headers)
  }
  
  # Process all TXSScan files
  txsscan_combined <- txsscan_files %>% 
    lapply(parse_all_systems_file) %>% 
    bind_rows() %>% 
    separate(hit_id, into = c("contig.num", "orf.num"), sep="_", extra = "drop", remove = FALSE) %>%  # Retain hit_id column
    mutate(sys_id2 = sub(".*\\.prot_", "", sys_id)) %>%
    select(Sample, hit_id, contig.num, orf.num, gene_name, sys_id2, hit_begin_match, hit_end_match) %>% 
    mutate(across(c(hit_begin_match, hit_end_match), as.integer)) %>%  # Ensure positions are numeric
    rename(gene = gene_name, type = sys_id2, start = hit_begin_match, end = hit_end_match)
  
  # Add Prodigal data
  processed_results <- txsscan_combined %>%
    group_split(Sample) %>%
    map_df(function(txsscan_data) {
      # Find the matching prodigal file for the given sample
      sample <- unique(txsscan_data$Sample)
      prodigal_file <- prodigal_files[basename(dirname(prodigal_files)) == sample]
      if (length(prodigal_file) != 1) {
        warning(sprintf("No unique prodigal file found for sample '%s'. Skipping.", sample))
        return(txsscan_data)
      }
      
      # Parse prodigal file
      prodigal_mapping <- parse_prodigal_headers(prodigal_file)
      
      # Map Prodigal headers to TXSScan results using hit_id (orf_id)
      txsscan_data %>%
        left_join(prodigal_mapping, by = c("hit_id" = "orf_id")) %>%  # Correctly join using hit_id
        mutate(
          start = ifelse(
            strand == "positive", contig_start + start - 1, contig_end - start + 1
          ),
          end = ifelse(
            strand == "positive", contig_start + end - 1, contig_end - end + 1
          )
        ) %>%
        select(Sample, contig.num, gene, type, start, end, strand)  # Keep relevant columns
    })
  
  return(processed_results)
}


txsscan <- process_txsscan_and_prodigal_files(txsscan_files, prodigal_files)
dlog("DONE:  TXSScan + Prodigal coordinate mapping (", nrow(txsscan), " hits)")


#################################################
###           TXXScan Models
#################################################

dlog("START: TXSScan system models (all_systems.txt)")
process_file <- function(file_path) {
  # Read the file into a character vector
  lines <- readLines(file_path)
  
  # Get the sample name from the sequence-db path
  sequence_db_line <- lines[grep("--sequence-db", lines)]
  sample_name <- str_match(sequence_db_line, ".*unicycler/(.+?)/.*\\.faa")[, 2]
  
  if (is.na(sample_name)) {
    warning(sprintf("Could not parse sample name from '%s'; skipping.", file_path))
    return(NULL)
  }
  
  # Extract systems where wholeness > 0.50
  model_lines <- grep("^model = ", lines, value = TRUE)
  wholeness_lines <- grep("^wholeness = ", lines, value = TRUE)
  
  # Combine models with their wholeness values
  systems_info <- data.frame(
    model = str_match(model_lines, "model = .*/(.+)$")[, 2],  # Extract only the part after the last '/'
    wholeness = as.numeric(str_match(wholeness_lines, "wholeness = ([0-9.]+)")[, 2])
  )
  
  # Filter identified systems based on wholeness > 0.50
  valid_models <- systems_info %>%
    filter(wholeness > 0.50) %>% 
    pull(model)
  
  # Return a dataframe for this sample
  data.frame(
    Sample = sample_name,
    TXSScan_Models_Present = paste(valid_models, collapse = " | ")
  )
}

# Process all_systems.txt files in a list
process_all_files <- function(file_list) {
  bind_rows(Filter(Negate(is.null), map(file_list, process_file)))
}


# Process files and create the final dataframe
txsscan_mod <- process_all_files(txsscan_files2)
dlog("DONE:  TXSScan system models (", nrow(txsscan_mod), " samples with systems)")

#################################################
###        MobileElementFinder Data
#################################################

dlog("START: MobileElementFinder")
mef_raw <- do.call(bind_rows, lapply(mef_files, function(path) {
  x <- tryCatch(
    read_csv(path, skip = 5, col_types = cols(.default = col_character()), show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(x) == 0) return(NULL)
  x$Sample <- basename(dirname(path))
  x
}))

if (!is.null(mef_raw) && ncol(mef_raw) > 0) {
  mef <- mef_raw %>%
    mutate(contig.num = contig, start = as.integer(start), end = as.integer(end)) %>%
    select(Sample, contig.num, name, type, start, end) %>%
    rename(gene=name) %>%
    mutate(strand="positive") %>%
    mutate(type = str_replace_all(type,
                                    c("mite" = "miniature inverted repeat",
                                      "ice"  = "integrative conjugative element",
                                      "ime"  = "integrative mobilizable element")))
  dlog("DONE:  MobileElementFinder (", nrow(mef), " elements)")
} else {
  mef <- tibble(Sample=character(), contig.num=character(), gene=character(),
                type=character(), start=integer(), end=integer(), strand=character())
  dlog("DONE:  MobileElementFinder (no elements found)")
}
               

#################################################
###           Contig mapping
#################################################

dlog("START: Contig map (MOB-recon contig_report)")
contigs <- do.call(bind_rows, lapply(contig_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols(.default = col_character())), error=function(e) NULL)
  if (is.null(x) || nrow(x) == 0) return(NULL)
  x$Sample <- basename(dirname(f))
  x
})) %>% 
  mutate(contig.num = str_split_i(contig_id, " ", 1)) %>% 
  select(Sample, contig_id, circularity_status, molecule_type, primary_cluster_id, contig.num)

if (ncol(mobtyper) > 0) {
  contigs <- contigs %>%
    left_join(mobtyper %>% select(Sample, primary_cluster_id, predicted_mobility, predicted_host_range_overall_name),
              by = c("Sample", "primary_cluster_id"))
}
dlog("DONE:  Contig map (", nrow(contigs), " contigs)")

dlog("START: AMR/MGE contig annotation merge")

if (!is.null(amrfinderplus) && ncol(amrfinderplus) > 0) {
  amrfinderplus_merge <- amrfinderplus %>% select(Sample, `Contig id`, `Element symbol`, Type, Start, Stop, Strand) %>%
    rename(contig.num=`Contig id`,
          gene=`Element symbol`,
          start=Start,
          end=Stop,
          strand=Strand,
          type=Type) %>%
    mutate(strand=ifelse(strand=="+", "positive", "negative"),
           start=as.integer(start), end=as.integer(end))
} else {
  amrfinderplus_merge <- tibble()
}

amr_mge_total <- bind_rows(amrfinderplus_merge, mef, txsscan, mobtyper_blast)

contig_map <- contigs %>% left_join(amr_mge_total, by=c("Sample", "contig.num")) %>% select(-contig.num)
dlog("DONE:  AMR/MGE contig annotation merge (", nrow(contig_map), " rows)")        


###############################################
###       Make plasmid summary table
###############################################

#plasmid_summary <- left_join(mobtyper, amr_plas %>% select(-Sample), by="Name") %>% 
#                    select(Name, Sample, primary_cluster_id, secondary_cluster_id, num_contigs, 
#                            size, gc, `rep_type(s)`, `relaxase_type(s)`, mpf_type, `orit_type(s)`,
#                            predicted_mobility, predicted_host_range_overall_rank, 
#                            predicted_host_range_overall_name, observed_host_range_ncbi_rank, 
#                            observed_host_range_ncbi_name, reported_host_range_lit_rank, 
#                            reported_host_range_lit_name, Start:`Closest reference name`)


###############################################
###         Gather all info
###############################################

dlog("START: Final summary table join")
summary <- gtdbtk %>% select(Sample, Species, closest_genome_reference) %>%
              left_join(mlst %>% select(Sample, Sequence_Type), by="Sample") %>%
              left_join(
                salmonella_serotype %>%
                  select(Sample, Predicted_Serotype, Salmonella_Serotype_Source,
                         SeqSero2_Serotype, SISTR_Serovar, SISTR_QC_Status, Serotype_Agreement),
                by="Sample"
              ) %>%
              left_join(afp_genotype, by="Sample") %>%
              left_join(resfinder_genotype, by="Sample") %>%
              left_join(phenotype, by="Sample") %>%
              left_join(pf %>% select(Sample, Pointfinder_Mutation), by="Sample") %>%
              left_join(mobtyper_summary, by="Sample") %>%
              left_join(txsscan_mod, by="Sample")
dlog("DONE:  Final summary table (", nrow(summary), " samples)")

df_names <- list(
  "summary_out" = summary,
  "AMRFinderPlus" = amrfinderplus,
  "assembly_QA" = assembly_stats,
  "MOBrecon_summary" = mobtyper,
  "salmonella_serotype" = salmonella_serotype_simple,
  "salmonella_seqsero2" = seqsero_simple,
  "salmonella_sistr" = sistr_simple,
  "ecoli_serotype" = ectyper_simple
)


###############################################
###             Write outputs
###############################################

dlog("START: Writing outputs")

write_xlsx(df_names, out_file1)
write.csv(contig_map, out_file2, row.names=FALSE)
dlog("DONE:  Outputs written — data_summary.R complete")
