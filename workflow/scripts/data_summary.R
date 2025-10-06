#!/usr/bin/env Rscript

library(tidyverse)
library(readr)

# Input files
checkm_file  <- snakemake@input[["checkm"]]
checkm_stats  <- snakemake@input[["checkm_stats"]]
gtdbtk_file  <- snakemake@input[["gtdbtk"]]
mlst_file    <- snakemake@input[["mlst"]]
resfinder_files <- snakemake@input[["resfinder_files"]]
amr_cr_files   <- snakemake@input[["amr_cr_files"]]
amr_plas_done <- snakemake@input[["amr_plas_done"]]
seqsero_files  <- snakemake@input[["seqsero_files"]]
mobtyper_files   <- snakemake@input[["mobtyper_files"]]
coverage_files  <- snakemake@input[["coverage_files"]]

# Output files

out_file1 <- snakemake@output[[1]]
out_file2 <- snakemake@output[[2]]
out_file3 <- snakemake@output[[3]]
out_file4 <- snakemake@output[[4]]

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
                    left_join(checkm %>% select(Sample, Completeness, Contamination, Strain_Heterogeneity), by="Sample")
                    
                    
assembly_stats <- assembly_stats %>% mutate(QA = ifelse(`N50 (contigs)` > 20000 & `# contigs` < 500 & XCoverage >= 30, "pass", "fail"))

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


###############################################
### Gather MOB-suite MGE data
###############################################

mobtyper <- do.call(bind_rows, lapply(mobtyper_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols()), error=function(e) NULL)
  if (is.null(x)) return(NULL)
  x$Sample <- str_match(f, "(?<=samples/)(.*?)(?=/mob-suite)")[, 2]
  x
}))

mobtyper <- mobtyper[, c(ncol(mobtyper), 1:(ncol(mobtyper) - 1))]
mobtyper <- mobtyper %>% mutate(Name = paste0(Sample, ".plasmid_", str_extract(sample_id, "(?<=assembly:)\\w+")))

mobtyper_summary <- mobtyper %>% group_by(Sample) %>% summarise(n_Plasmid = length(sample_id), 
                                                                Plasmid_Rep_Types = paste(`rep_type(s)`, collapse = ", "), 
                                                                Plasmid_Relaxase_Types = paste(`relaxase_type(s)`, collapse = ", "))


plasmid_summary <- left_join(mobtyper, amr_plas %>% select(-Sample), by="Name") %>% relocate(Name)

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
###       Make summary table
###############################################

summary <- gtdbtk %>%
  select(Sample, Species) %>%
  left_join(mlst %>% select(Sample, Scheme, Sequence_Type), by="Sample") %>%
  left_join(seqsero %>% select(Sample, Predicted_Serotype), by="Sample") %>%
  left_join(afp_genotype, by="Sample") %>%
  left_join(resfinder_genotype, by="Sample") %>%
  left_join(phenotype, by="Sample") %>%
  left_join(amr_cr_summary, by="Sample") %>%
  left_join(mobtyper_summary, by="Sample") %>%
  left_join(amr_plas_summary, by="Sample") 

###############################################
###       Write outputs
###############################################
 
write_csv(summary, out_file1)
write_csv(seqsero_simple, out_file2)
write_csv(plasmid_summary, out_file3)
write_csv(assembly_stats, out_file4)



