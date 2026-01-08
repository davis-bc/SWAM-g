#!/usr/bin/env Rscript

libraries <- c("tidyverse", "writexl")

invisible(lapply(libraries, function(x) {
  suppressMessages(suppressWarnings(library(x, character.only = T)))
  }))


# Input files

checkm_file          <- snakemake@input[["checkm"]]
checkm_stats         <- snakemake@input[["checkm_stats"]]
gtdbtk_file          <- snakemake@input[["gtdbtk"]]
mlst_file            <- snakemake@input[["mlst"]]
resfinder_files      <- snakemake@input[["resfinder_files"]]
pf_files             <- snakemake@input[["pf_files"]]
amr_cr_files         <- snakemake@input[["amr_cr_files"]]
amr_plas_done        <- snakemake@input[["amr_plas_done"]]
seqsero_files        <- snakemake@input[["seqsero_files"]]
ectyper_files        <- snakemake@input[["ectyper_files"]]
mobtyper_files       <- snakemake@input[["mobtyper_files"]]
coverage_files       <- snakemake@input[["coverage_files"]]
mef_files            <- snakemake@input[["mef_files"]]
contig_files         <- snakemake@input[["contig_files"]]
mobtyper_blast_files <- snakemake@input[["mobtyper_blast_files"]]
txsscan_files        <- snakemake@input[["txsscan_files"]]
txsscan_files2       <- snakemake@input[["txsscan_files2"]]

# Output files

out_file1 <- snakemake@output[[1]]
out_file2 <- snakemake@output[[2]]


###############################################
###      Read in pre-generated tables
###############################################

lines <- readLines(checkm_file)
checkm   <- read_table(lines, skip = 3, n_max = length(lines) - 4, col_names = FALSE)[1:15]
names(checkm) <- c("Sample", "marker", "lineage", "genomes", "markers", "marker.sets", 
                   "0", "1", "2", "3", "4", "5+", 
                   "Completeness", "Contamination", "Strain_Heterogeneity")
checkm$Sample <- gsub("\\.chromosome$", "", checkm$Sample)

gtdbtk <- read_tsv(gtdbtk_file, col_types = cols())
gtdbtk$Sample <- gsub("\\.chromosome$", "", gtdbtk$user_genome)
gtdbtk$Species <- sub(".*s__", "", gtdbtk$classification)

mlst <- read_table(mlst_file, col_names = F, col_types = cols())
names(mlst)[1:3] <- c("Sample", "Scheme", "Sequence_Type")
mlst$Sample <- sub(".*/(.*)\\.chromosome\\.fasta", "\\1", mlst$Sample)

###############################################
###       Gather coverage data
###############################################

coverage <- do.call(rbind, lapply(coverage_files, function(f) 
  read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)))

colnames(coverage) <- c("Sample", "XCoverage")

########################################################
###      Parse ridiculous checkm bin_stats.analyze file
########################################################

# Read lines from file
lines <- read_lines(checkm_stats)

# Function to parse a dictionary string into a named list
parse_dict <- function(dict_string) {
  dict_string <- str_remove_all(dict_string, "^\\{|\\}$")
  kv_pairs <- str_split(dict_string, ",\\s*")[[1]]
  map(kv_pairs, function(pair) {
    parts <- str_split(pair, ":\\s*", n = 2)[[1]]
    key <- str_remove_all(parts[1], "^'|'$")
    value <- str_remove_all(parts[2], "^'|'$")
    # Convert to numeric if possible
    if (str_detect(value, "^[0-9.]+$")) value <- as.numeric(value)
    else if (str_detect(value, "^[0-9]+$")) value <- as.integer(value)
    set_names(list(value), key)
  }) %>% flatten()
}

# Function to parse a line
parse_line <- function(line) {
  filename <- str_trim(str_remove(line, "\\{.*\\}$"))
  dict_str <- str_extract(line, "\\{.*\\}$")
  dict <- parse_dict(dict_str)
  tibble(file = filename, !!!dict)
}

# Parse all lines and bind into a tibble
assembly_stats <- map_df(lines, parse_line)

assembly_stats$Sample <- gsub("\\.chromosome$", "", assembly_stats$file)

assembly_stats <- assembly_stats %>% relocate(Sample) %>% select(-file) %>%
                  left_join(coverage, by="Sample") %>%
                  left_join(checkm %>% select(Sample, Completeness, Contamination, Strain_Heterogeneity), by="Sample") %>% 
                  mutate(QA = ifelse(`N50 (contigs)` > 20000 & `# contigs` < 500 & XCoverage >= 30, "pass", "fail")) %>% 
                  select(Sample, `Genome size`, GC, `# contigs`,  `N50 (contigs)`, `Coding density`, XCoverage, Completeness, Contamination, Strain_Heterogeneity, QA)

###############################################
###   Gather AMRFinderPlus chromosome results
###############################################

expected_types <- c("AMR", "STRESS", "VIRULENCE")

amr_cr <- do.call(bind_rows, lapply(amr_cr_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols(`Contig id` = "c")), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) {return(NULL)}
  x$Sample <- gsub("\\.afp.tsv$", "", basename(f))
  x
}))

amr_cr_summary <- amr_cr %>%
  group_by(Sample, Type) %>%
  summarize(genes = paste(`Element symbol`, collapse = ", "), .groups = 'drop') %>%
  tidyr::pivot_wider(
    names_from = Type,
    values_from = genes,
    names_prefix = "Chromosome_",
    values_fill = list(genes = NA)
  )

# Ensure all expected columns are present
for (col in paste0("Chromosome_", expected_types)) {
  if (!col %in% names(amr_cr_summary)) amr_cr_summary[[col]] <- NA
}

amr_cr_summary <- amr_cr_summary %>%
  select(Sample, Chromosome_AMR, Chromosome_STRESS, Chromosome_VIRULENCE)

###############################################
###   Gather AMRFinderPlus plasmid results
###############################################

amr_plas_files <- unique(unlist(lapply(amr_plas_done, function(marker) {
  dir <- dirname(marker)
  # List all .plasmid*.afp.tsv files for this sample
  list.files(dir, pattern = "*plasmid*", full.names = TRUE)
})))


amr_plas <- do.call(bind_rows, lapply(amr_plas_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols(`Contig id` = "c")), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) {return(NULL)}
  x$Sample <- sub("\\.plasmid.*", "", basename(f))
  x
}))

amr_plas_summary <- amr_plas %>%
  group_by(Sample, Type) %>%
  summarize(genes = paste(`Element symbol`, collapse = ", "), .groups = 'drop') %>%
  tidyr::pivot_wider(
    names_from = Type,
    values_from = genes,
    names_prefix = "Plasmid_",
    values_fill = list(genes = NA)
  )

for (col in paste0("Plasmid_", expected_types)) {
  if (!col %in% names(amr_plas_summary)) amr_plas_summary[[col]] <- NA
}

amr_plas_summary <- amr_plas_summary %>%
  select(Sample, Plasmid_AMR, Plasmid_STRESS, Plasmid_VIRULENCE)
  
  
afp_genotype <- bind_rows(amr_cr, amr_plas) %>% filter(Type == "AMR") %>% group_by(Sample) %>%
                summarise(AMRFinder_Genotype = paste(sort(unique(`Element symbol`)), collapse = ", "))
                
amrfinderplus_all <- rbind(amr_cr, amr_plas) %>% relocate(Sample) %>% select(-Name)

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
  x <- tryCatch(read_tsv(f, col_types = cols()), error=function(e) NULL)
  if (is.null(x)) return(NULL)
  x$Sample <- basename(dirname(f))
  x
})) %>% 
    relocate(Sample) %>% 
    mutate(Name = paste0(Sample, ".plasmid_", str_extract(sample_id, "(?<=assembly:)\\w+"))) %>%
    separate(sample_id, into=c("meh", "primary_cluster_id"), sep = ":")

mobtyper_summary <- mobtyper %>% 
                    group_by(Sample) %>% 
                    summarise(n_Plasmid_Clusters = length(primary_cluster_id), 
                                Plasmid_Rep_Types = paste(`rep_type(s)`, collapse = " | "), 
                                Plasmid_Relaxase_Types = paste(`relaxase_type(s)`, collapse = " | "))

plasmid_summary <- left_join(mobtyper, amr_plas %>% select(-Sample), by="Name") %>% 
                    select(Name, Sample, primary_cluster_id, secondary_cluster_id, num_contigs, 
                            size, gc, `rep_type(s)`, `relaxase_type(s)`, mpf_type, `orit_type(s)`,
                            predicted_mobility, predicted_host_range_overall_rank, 
                            predicted_host_range_overall_name, observed_host_range_ncbi_rank, 
                            observed_host_range_ncbi_name, reported_host_range_lit_rank, 
                            reported_host_range_lit_name, Start:`Closest reference name`)

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


parse_all_systems_file <- function(filepath) {
  # Read the file while skipping comment lines (lines starting with '#')
  data <- read_tsv(filepath, comment = "#", col_types = cols(.default = "c"))
  
  # Extract the corresponding system name for each block
  lines <- readLines(filepath)
  
  # Get the system descriptions from the comments preceding each "replicon" block
  system_descriptions <- tibble(line = lines) %>%
    filter(str_starts(line, "#")) %>%
    mutate(block = cumsum(str_detect(lead(line, default = ""), "^replicon"))) %>%
    filter(!is.na(block), block > 0) %>%
    group_by(block) %>%
    summarise(system_name = str_remove(last(line), "^# "), .groups = "drop")
  
  # Add block-level system_name to data
  data <- data %>%
    mutate(block = cumsum(str_detect(replicon, "^replicon"))) %>%
    left_join(system_descriptions, by = "block") %>%
    select(-block) %>%
    filter(!replicon == "replicon") %>%
    mutate(Sample = sub(".prot", "", replicon)) %>%
    relocate(Sample)
  
  return(data)
}

txsscan <- txsscan_files %>% 
            lapply(parse_all_systems_file) %>% 
            bind_rows() %>% 
            separate(hit_id, into = c("contig.num", "orf.num"), sep="_") %>%
            select(Sample, contig.num, gene_name, sys_id, hit_begin_match, hit_end_match) %>% 
            mutate(strand="unknown", sys_id= sub(".*\\.prot_", "", sys_id)) %>%
            rename(gene=gene_name, type=sys_id, start=hit_begin_match, end=hit_end_match)


#################################################
###           TXXScan Models
#################################################


# Define a function to process a single file
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
    model = str_match(model_lines, "model = .*?/(.+)$")[, 2],
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

afp_total <- rbind(amr_cr %>% select(Sample, `Contig id`, `Element symbol`, Type, Start, Stop, Strand),
                   amr_plas %>% select(Sample, `Contig id`, `Element symbol`, Type, Start, Stop, Strand)) %>%
              rename(contig.num=`Contig id`,
                     gene=`Element symbol`,
                     start=Start,
                     end=Stop,
                     strand=Strand,
                     type=Type) %>%
              mutate(strand=ifelse(strand=="+", "positive", "negative"))

amr_mge_total <- rbind(afp_total, mef, txsscan, mobtyper_blast)

contig_map <- contigs %>% left_join(amr_mge_total, by=c("Sample", "contig.num")) %>% select(-contig.num)         


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
              left_join(amr_cr_summary, by="Sample") %>%
              left_join(mobtyper_summary, by="Sample") %>%
              left_join(amr_plas_summary, by="Sample") %>%
              left_join(txsscan_mod, by="Sample")

df_names <- list(
  "summary_out" = summary,
  "AMRFinderPlus" = amrfinderplus_all,
  "assembly_QA" = assembly_stats,
  "MOBrecon_summary" = plasmid_summary,
  "salmonella_serotype" = seqsero_simple,
  "ecoli_serotype" = ectyper_simple
)


###############################################
###             Write outputs
###############################################

write_xlsx(df_names, out_file1)
write.csv(contig_map, out_file2)

