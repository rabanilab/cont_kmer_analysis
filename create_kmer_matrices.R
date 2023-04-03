library(kmer)
library(adegenet)
library(dplyr)
library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# test if there are two arguments: if not, return an error
if (length(args) != 2) {
  stop("Please supply fasta file and kmer length", call. = FALSE)
}

if (file.exists(args[1])) {
  fasta.path <- args[1]
} else {
  stop(paste("file", args[1], "doesn't exist"))
}

if (!is.na(as.numeric(args[2]))) {
  k.len <- as.numeric(args[2])
} else {
  stop(paste(args[2], "isn't a valid kmer length"))
}

# fastaPath <- "C:\\Users\\liorf\\Documents\\cluster_all_bg_filtered_longest.fa"
# k.len <- 6

all.obj <- ape::read.dna(fasta.path, format = "fasta")

counts <- kcount(all.obj, k = k.len, residues = "DNA")
counts_df <- data.frame(counts) %>%
  rownames_to_column('ids') %>%
  separate(ids, c('ensembl_gene_id', 'ensembl_gene_id_uniq', 'ensembl_transcript_id', 'ensembl_transcript_id_uniq', 'gene_name'), sep = '\\|')
# took up more memory
# %>%
#   pivot_longer(
#     cols = -c('ensembl_gene_id', 'ensembl_gene_id_uniq', 'ensembl_transcript_id', 'ensembl_transcript_id_uniq', 'gene_name'),
#     names_to = 'kmer',
#     values_to = 'occurrence'
#     )

outpath <- gsub('\\.\\w+$', paste0('_', k.len, "kmer.csv"), fasta.path)
fwrite(counts_df, outpath)
