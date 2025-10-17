#!/usr/bin/env Rscript

libraries <- c("tidyverse")

invisible(lapply(libraries, function(x) {
  suppressMessages(suppressWarnings(library(x, character.only = T)))
  }))


# Input files
contig_file  <- snakemake@input[["contigs"]]
afp_file     <- snakemake@input[["afp"]]
mef_file     <- snakemake@input[["mef"]]

# Output files

out_file1 <- snakemake@output[[1]]

contigs <- read.table(contig_file, header = T, sep = "\t")
afp     <- read.table(afp_file, header = T, sep = "\t")
mef     <- read.csv(mef_file, header = T, skip = 5)

contigs$Sample <- basename(dirname(contig_file))

contigs <- contigs %>% mutate(contig.num = str_split_i(contig_id, " ", 1)) %>% 
                select(
                    Sample, 
                    contig_id, 
                    molecule_type, 
                    circularity_status, 
                    contig.num
                      )
                       
afp <- afp %>%  select(
                    Contig.id, 
                    Element.symbol,
                    Type,
                    Start, 
                    Stop, 
                    Strand
                      ) %>%
                rename(
                    contig.num=Contig.id, 
                    gene=Element.symbol, 
                    start=Start, 
                    end=Stop, 
                    strand=Strand, 
                    type=Type
                     ) %>%
                mutate(strand=ifelse(strand=="+", "positive", "negative"))
                
mef <- mef %>% mutate(contig.num = str_split_i(contig, " ", 1)) %>%
                select(
                    contig.num,
                    name,
                    type,
                    start,
                    end
                       ) %>%
                rename(gene=name) %>%
                mutate(strand="positive")
                

