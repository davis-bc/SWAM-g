#!/usr/bin/env Rscript

library(tidyverse)
library(readr)

# Input files
checkm_file  <- snakemake@input[["checkm"]]
gtdbtk_file  <- snakemake@input[["gtdbtk"]]
mlst_file    <- snakemake@input[["mlst"]]
amr_cr_files   <- snakemake@input[["amr_cr_files"]]
amr_plas_done <- snakemake@input[["amr_plas_done"]]
seqsero_files  <- snakemake@input[["seqsero_files"]]
sample_dirs   <- snakemake@input[["sample_dirs"]]

# Output file

out_file <- snakemake@output[[1]]

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
###   Gather AMRFinderPlus chromosome results
###############################################

amr_cr <- do.call(bind_rows, lapply(amr_cr_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols(`Contig id` = "c")), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) {return(NULL)}
  x$Sample <- gsub("\\.afp.tsv$", "", basename(f))
  x
}))

amr_cr_summary <- amr_cr %>% group_by(Sample, Type) %>% summarize(genes = paste(`Element symbol`, collapse = ","), .groups = 'drop') %>% spread(Type, genes)
names(amr_cr_summary) <- c("Sample", "Chromosome_AMR", "Chromosome_Stress", "Chromosome_Virulence")

###############################################
###   Gather AMRFinderPlus plasmid results
###############################################

amr_plas_files <- unlist(lapply(amr_plas_done, function(marker) {
  dir <- dirname(marker)
  # List all .plasmid*.afp.tsv files for this sample
  list.files(dir, pattern = "*plasmid*", full.names = TRUE)
}))

amr_plas <- do.call(bind_rows, lapply(amr_plas_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols(`Contig id` = "c")), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) {return(NULL)}
  x$Sample <- sub("\\.plasmid.*", "", basename(f))
  x
}))

amr_plas_summary <- amr_plas %>% group_by(Sample, Type) %>% summarize(genes = paste(`Element symbol`, collapse = ","), .groups = 'drop') %>% spread(Type, genes)
names(amr_plas_summary) <- c("Sample", "Plasmid_AMR", "Plasmid_Stress", "Plasmid_Virulence")

###############################################
###   Gather all Serotyping results
###############################################

seqsero <- do.call(bind_rows, lapply(seqsero_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols()), error=function(e) NULL)
  if (is.null(x)) return(NULL)
  x
}))

seqsero$Sample <- gsub("\\.chromosome.fasta$", "", seqsero$`Sample name`)
names(seqsero)[8:9] <- c("Predicted_Serotype", "Predicted_Antigenic_Profile")

###############################################
### Gather MOB-suite MGE data
###############################################

mobtyper_files <- list.files(sample_dirs, pattern="mobtyper_results.txt", full.names=TRUE, recursive = 2)

mobtyper <- do.call(bind_rows, lapply(mobtyper_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols()), error=function(e) NULL)
  if (is.null(x)) return(NULL)
  x$Sample <- str_match(f, "(?<=Samples/)(.*?)(?=/mob-suite)")[, 2]
  x
}))

mobtyper_summary <- mobtyper %>% group_by(Sample) %>% summarise(n_Plasmid = length(sample_id), 
                                                                Plasmid_Rep_Types = paste(`rep_type(s)`, collapse = ","), 
                                                                Plasmid_Relaxase_Types = paste(`relaxase_type(s)`, collapse = ","))
###############################################
###       Make summary table
###############################################

summary <- gtdbtk %>% select(Sample, Species) %>%
           left_join(mlst %>% select(Sample, Scheme, Sequence_Type), by="Sample") %>%
           left_join(seqsero %>% select(Sample, Predicted_Serotype, Predicted_Antigenic_Profile), by="Sample") %>%
           left_join(amr_cr_summary, by="Sample") %>% 
           left_join(mobtyper_summary, by="Sample") %>%
           left_join(amr_plas_summary, by="Sample") %>%
           left_join(checkm %>% select(Sample, Completeness, Contamination, Strain_Heterogeneity), by="Sample")

write_tsv(summary, out_file)





