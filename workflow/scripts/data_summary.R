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


###############################################
###      Read in pre-generated tables
###############################################

gtdbtk <- read_tsv(mash_file, col_types = cols())
gtdbtk$Sample <- gsub("\\.chromosome$", "", gtdbtk$user_genome)
gtdbtk$Species <- sub(".*s__", "", gtdbtk$classification)

mlst <- read_table(mlst_file, col_names = F, col_types = cols())
names(mlst)[1:3] <- c("Sample", "Scheme", "Sequence_Type")
mlst$Sample <- sub(".*/(.*)\\.fasta", "\\1", mlst$Sample)

########################################################
###       Gather checkm2 and coverage data, assess assembly quality
########################################################

checkm2 <- do.call(rbind, lapply(seq_along(checkm2_files), function(i) {
  if (i == 1) {
    read.table(checkm2_files[i], header = TRUE, sep = "\t", 
               stringsAsFactors = FALSE, row.names = NULL)
  } else {
    # Get column names from first file
    col_names <- colnames(read.table(checkm2_files[1], header = TRUE, 
                                     sep = "\t", nrows = 0))
    read.table(checkm2_files[i], header = FALSE, sep = "\t", 
               stringsAsFactors = FALSE, skip = 1, 
               col.names = col_names, row.names = NULL)
  }
}))

coverage <- do.call(rbind, lapply(coverage_files, function(f) 
  read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)))

colnames(coverage) <- c("Sample", "XCoverage")

assembly_stats <- checkm2 %>% 
                   rename(Sample=Name) %>% 
                   left_join(coverage, by="Sample") %>%
                   mutate(QA = ifelse(Contig_N50 > 20000 & Total_Contigs < 500 & XCoverage >= 30, "pass", "fail"))



###############################################
###   Gather AMRFinderPlus results
###############################################

expected_types <- c("AMR", "STRESS", "VIRULENCE")

amrfinderplus <- do.call(bind_rows, lapply(afp_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols(`Contig id` = "c")), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) {return(NULL)}
  x$Sample <- gsub("\\.afp.tsv$", "", basename(f))
  x
}))

afp_genotype <- amrfinderplus %>% filter(Type == "AMR") %>% group_by(Sample) %>%
                summarise(AMRFinder_Genotype = paste(sort(unique(`Element symbol`)), collapse = ", "))


###############################################
###   Gather all Serotyping results
###############################################

seqsero <- do.call(bind_rows, lapply(seqsero_files, function(f) {
  x <- tryCatch(
    read_tsv(f, col_types = cols(.default = col_character())),
    error=function(e) NULL
  )
  if (is.null(x)) return(NULL)
  x
}))
  
seqsero$Sample <- basename(seqsero$`Output directory`)
names(seqsero)[8:9] <- c("Predicted_Antigenic_Profile", "Predicted_Serotype")

seqsero_simple <- seqsero %>% filter(grepl("Salmonella", `Predicted identification`)) %>%
  relocate(Sample) %>% select(-c(`Sample name`, `Output directory`))
  
  
ectyper <- do.call(bind_rows, lapply(ectyper_files, function(f) {
  x <- tryCatch(
    read_tsv(f, col_types = cols(.default = col_character())),
    error=function(e) NULL
  )
  if (is.null(x)) return(NULL)
  x$Sample <- basename(dirname(f))
  x
}))
  
ectyper_simple <- ectyper %>% filter(grepl("Escherichia", Species)) %>% relocate(Sample) %>% select(-Name)

###############################################
###       Gather MOB-suite plasmid data
###############################################

mobtyper <- do.call(bind_rows, lapply(mobtyper_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols(`associated_pmid(s)` = col_character())), error=function(e) NULL)
  if (is.null(x)) return(NULL)
  x
})) %>% 
    separate(sample_id, into=c("Sample", "primary_cluster_id"), sep = ":") %>%
    mutate(PlasmidID = paste0(Sample, ".plasmid_", primary_cluster_id)) %>%
    relocate(Sample, PlasmidID) 

mobtyper_summary <- mobtyper %>% 
                    group_by(Sample) %>% 
                    summarise(n_Plasmid_Clusters = length(primary_cluster_id), 
                                Plasmid_Rep_Types = paste(`rep_type(s)`, collapse = " | "), 
                                Plasmid_Relaxase_Types = paste(`relaxase_type(s)`, collapse = " | "))



mobtyper_blast <- do.call(bind_rows, lapply(mobtyper_blast_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols()), error=function(e) NULL)
  if (is.null(x)) return(NULL)
  x$Sample <- basename(dirname(f))
  x
})) %>% 
      separate(qseqid, into = c("accession", "gene"), sep="\\|") %>%
      mutate(contig.num = str_split_i(sseqid, " ", 1), sstrand=ifelse(sstrand=="minus", "negative", "positive")) %>%
      select(Sample, contig.num, gene, biomarker, sstart, send, sstrand) %>%
      rename(type=biomarker, start=sstart, end=send, strand=sstrand)

###############################################
###      Parse Resfinder phenotypes
###############################################

safe_read_resfinder <- function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols()), error=function(e) NULL)
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

phenotype <- resfinder %>% group_by(Sample) %>%
             summarise(phenotype = list(unlist(str_split(Phenotype, ",\\s*")))) %>% 
             mutate(Resfinder_Predicted_Phenotype = map_chr(phenotype, ~ paste(sort(unique(trimws(.))), collapse = ", "))) %>%
             select(-phenotype)
             
resfinder_genotype <- resfinder %>% group_by(Sample) %>% 
  summarise(Resfinder_Genotype = paste(sort(unique(`Resistance gene`)), collapse = ", "))      
  

###############################################
###        PointFinder + DisinFinder
###############################################

pf <- do.call(
  bind_rows,
  lapply(pf_files, function(f) {
    x <- tryCatch(read_tsv(f, col_types = cols()), error = function(e) NULL)
    if (is.null(x)) return(NULL)
    x$Sample <- basename(dirname(f))
    # Coerce `PMID` (or other problematic columns) to character for consistency
    if ("PMID" %in% colnames(x)) {
      x$PMID <- as.character(x$PMID)
    }
    x
  })
) %>%
group_by(Sample) %>%
summarise(Pointfinder_Mutation = paste(Mutation, collapse = " | "))


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
    separate(hit_id, into = c("contig.num", "orf.num"), sep="_", remove = FALSE) %>%  # Retain hit_id column
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
      prodigal_file <- prodigal_files[grepl(sample, prodigal_files)]
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


#################################################
###           TXXScan Models
#################################################

process_file <- function(file_path) {
  # Read the file into a character vector
  lines <- readLines(file_path)
  
  # Get the sample name from the sequence-db path
  sequence_db_line <- lines[grep("--sequence-db", lines)]
  sample_name <- str_match(sequence_db_line, ".*unicycler/(.+?)/.*\\.faa")[, 2]
  
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
  bind_rows(map(file_list, process_file))
}


# Process files and create the final dataframe
txsscan_mod <- process_all_files(txsscan_files2)

#################################################
###        MobileElementFinder Data
#################################################

mef <- do.call(rbind, lapply(mef_files, function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, skip = 5, header = T)
  x$Sample <- basename(dirname(path))
  x
})) %>% 
    mutate(contig.num = str_split_i(contig, " ", 1)) %>%
    select(Sample, contig.num, name, type, start, end) %>%
    rename(gene=name) %>%
    mutate(strand="positive") %>% 
    mutate(type = str_replace_all(type, 
                                    c("mite" = "miniature inverted repeat", 
                                      "ice"  = "integrative conjugative element",
                                      "ime"  = "integrative mobilizable element")))
               

#################################################
###           Contig mapping
#################################################

contigs <- do.call(bind_rows, lapply(contig_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols()), error=function(e) NULL)
  if (is.null(x)) return(NULL)
  x$Sample <- basename(dirname(f))
  x
})) %>% 
  mutate(contig.num = str_split_i(contig_id, " ", 1)) %>% 
  select(Sample, contig_id, circularity_status, molecule_type, primary_cluster_id, contig.num) %>%
  left_join(mobtyper %>% select(Sample, primary_cluster_id, predicted_mobility, predicted_host_range_overall_name), by=c("Sample", "primary_cluster_id"))

amrfinderplus_merge <- amrfinderplus %>% select(Sample, `Contig id`, `Element symbol`, Type, Start, Stop, Strand) %>%
  rename(contig.num=`Contig id`,
        gene=`Element symbol`,
        start=Start,
        end=Stop,
        strand=Strand,
        type=Type) %>%
  mutate(strand=ifelse(strand=="+", "positive", "negative"))

amr_mge_total <- rbind(amrfinderplus_merge, mef, txsscan, mobtyper_blast)

contig_map <- contigs %>% left_join(amr_mge_total, by=c("Sample", "contig.num")) %>% select(-contig.num)        


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

summary <- gtdbtk %>% select(Sample, Species, closest_genome_reference) %>%
              left_join(mlst %>% select(Sample, Sequence_Type), by="Sample") %>%
              left_join(seqsero %>% select(Sample, Predicted_Serotype), by="Sample") %>%
              left_join(afp_genotype, by="Sample") %>%
              left_join(resfinder_genotype, by="Sample") %>%
              left_join(phenotype, by="Sample") %>%
              left_join(pf %>% select(Sample, Pointfinder_Mutation), by="Sample") %>%
              left_join(mobtyper_summary, by="Sample") %>%
              left_join(txsscan_mod, by="Sample")

df_names <- list(
  "summary_out" = summary,
  "AMRFinderPlus" = amrfinderplus,
  "assembly_QA" = assembly_stats,
  "MOBrecon_summary" = mobtyper,
  "salmonella_serotype" = seqsero_simple,
  "ecoli_serotype" = ectyper_simple
)


###############################################
###             Write outputs
###############################################

write_xlsx(df_names, out_file1)
write.csv(contig_map, out_file2, row.names=FALSE)

