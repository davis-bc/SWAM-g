#!/usr/bin/env Rscript

library(tidyverse)
library(readr)

# Input files
checkm_file  <- snakemake@input[["checkm"]]
gtdbtk_file  <- snakemake@input[["gtdbtk"]]
mlst_file    <- snakemake@input[["mlst"]]
quast_file   <- snakemake@input[["quast"]]
resfinder_files <- snakemake@input[["resfinder_files"]]
amr_cr_files   <- snakemake@input[["amr_cr_files"]]
amr_plas_done <- snakemake@input[["amr_plas_done"]]
seqsero_files  <- snakemake@input[["seqsero_files"]]
sample_dirs   <- snakemake@input[["sample_dirs"]]

# Output files

out_file1 <- snakemake@output[[1]]
out_file2 <- snakemake@output[[2]]
out_file3 <- snakemake@output[[3]]

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

quast <- read.table(quast_file, sep = "\t", header=T)
quast <- quast %>% gather(Sample, Value, 2:length(quast)) %>% 
         pivot_wider(names_from = Assembly, values_from = Value)
names(quast)[9:10] <- c("Genome_Length", "GC")

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


amr_plas_summary <- amr_plas %>% group_by(Sample, Type) %>% summarize(genes = paste(`Element symbol`, collapse = ","), .groups = 'drop') %>% spread(Type, genes)
names(amr_plas_summary) <- c("Sample", "Plasmid_AMR", "Plasmid_Stress", "Plasmid_Virulence")

###############################################
###   Gather all Serotyping results
###############################################

has_seqsero <- FALSE
seqsero <- tibble()

if (length(seqsero_files) > 0) {
  seqsero <- do.call(bind_rows, lapply(seqsero_files, function(f) {
    x <- tryCatch(read_tsv(f, col_types = cols()), error=function(e) NULL)
    if (is.null(x)) return(NULL)
    x
  }))
  if (nrow(seqsero) > 0) {
    has_seqsero <- TRUE
    seqsero$Sample <- gsub("\\.chromosome.fasta$", "", seqsero$`Sample name`)
    names(seqsero)[8:9] <- c("Predicted_Antigenic_Profile", "Predicted_Serotype")
  }
}

seqsero_simple <- seqsero %>% relocate(Sample) %>% select(-c(`Sample name`, `Output directory`))


###############################################
### Gather MOB-suite MGE data
###############################################

mobtyper_files <- list.files(sample_dirs, pattern="mobtyper_results.txt", full.names=TRUE, recursive = 2)

mobtyper <- do.call(bind_rows, lapply(mobtyper_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols()), error=function(e) NULL)
  if (is.null(x)) return(NULL)
  x$Sample <- str_match(f, "(?<=samples/)(.*?)(?=/mob-suite)")[, 2]
  x
}))

mobtyper <- mobtyper[, c(ncol(mobtyper), 1:(ncol(mobtyper) - 1))]
mobtyper <- mobtyper %>% mutate(Name = paste0(Sample, ".plasmid_", str_extract(sample_id, "(?<=assembly:)\\w+")))

mobtyper_summary <- mobtyper %>% group_by(Sample) %>% summarise(n_Plasmid = length(sample_id), 
                                                                Plasmid_Rep_Types = paste(`rep_type(s)`, collapse = ","), 
                                                                Plasmid_Relaxase_Types = paste(`relaxase_type(s)`, collapse = ","))


###############################################
###      Parse Resfinder phenotypes
###############################################

resfinder <- do.call(bind_rows, lapply(resfinder_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols()), error=function(e) NULL)
  if (is.null(x)) return(NULL)
  x$Sample <- basename(dirname(f))
  x
}))

phenotype <- resfinder %>% group_by(Sample) %>%
             summarise(phenotype = list(unlist(str_split(Phenotype, ",\\s*")))) %>% 
             mutate(CGE_Predicted_Phenotype = map_chr(phenotype, ~ paste(sort(unique(trimws(.))), collapse = ", "))) %>%
             select(-phenotype)

###############################################
###       Make summary table
###############################################

summary <- gtdbtk %>%
  select(Sample, Species) %>%
  left_join(mlst %>% select(Sample, Scheme, Sequence_Type), by="Sample") %>%
  left_join(phenotype, by="Sample") 

if (has_seqsero) {
  summary <- summary %>%
    left_join(seqsero %>% select(Sample, Predicted_Serotype, Predicted_Antigenic_Profile), by="Sample")
}

summary <- summary %>%
  left_join(amr_cr_summary, by="Sample") %>%
  left_join(mobtyper_summary, by="Sample") %>%
  left_join(amr_plas_summary, by="Sample") %>%
  left_join(quast %>% select(Sample, Genome_Length, GC, N50), by="Sample") %>%
  left_join(checkm %>% select(Sample, Completeness, Contamination, Strain_Heterogeneity), by="Sample")
  
plasmid_summary <- left_join(mobtyper, amr_plas %>% select(-Sample), by="Name") %>% relocate(Name)

write_csv(summary, out_file1)
write_csv(seqsero_simple, out_file2)
write_csv(plasmid_summary, out_file3)



