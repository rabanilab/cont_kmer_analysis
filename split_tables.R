library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Please supply fasta file and kmer length", call. = FALSE)
}

if (file.exists(args[1])) {
  base <- args[1]
} else {
  stop(paste("folder", args[1], "doesn't exist"))
}

base_raw <- paste0(base, '_tmp')

mer5 <- fread(file.path(base_raw, "cluster_all_filtered_longest_5kmer.csv")) %>%
  select(-ensembl_gene_id_uniq, -ensembl_transcript_id_uniq, -ensembl_transcript_id)
mer6 <- fread(file.path(base_raw, "cluster_all_filtered_longest_6kmer.csv")) %>%
  select(-ensembl_gene_id_uniq, -ensembl_transcript_id_uniq, -ensembl_transcript_id)
mer7 <- fread(file.path(base_raw, "cluster_all_filtered_longest_7kmer.csv")) %>%
  select(-ensembl_gene_id_uniq, -ensembl_transcript_id_uniq, -ensembl_transcript_id)
mer8 <- fread(file.path(base_raw, "cluster_all_filtered_longest_8kmer.csv")) %>%
  select(-ensembl_gene_id_uniq, -ensembl_transcript_id_uniq, -ensembl_transcript_id)

fwrite(
  mer5,
  file.path(
    base,
    "cluster_all_filtered_longest_5kmer_p1.csv"
  )
)

fwrite(
  mer6,
  file.path(
    base,
    "cluster_all_filtered_longest_6kmer_p1.csv"
  )
)
starts <- seq(1, 4^7, (4^7) / 4) + 1
ends <- seq(1, 4^7, (4^7) / 4) + ((4^7) / 4)
for (i in 1:4) {
  mer7 %>%
    select(colnames(mer7)[c(1, 2, seq(starts[i]+1,ends[i]+1))]) %>%
    fwrite(file.path(base,
                     paste0("cluster_all_filtered_longest_7kmer_p", i, ".csv")
    ))
}
starts <- seq(1, 4^8, (4^8) / 16) + 1
ends <- seq(1, 4^8, (4^8) / 16) + ((4^8) / 16)
for (i in 1:16) {
  mer8 %>%
    select(colnames(mer8)[c(1, 2, seq(starts[i]+1,ends[i]+1))]) %>%
    fwrite(file.path(base,
                     paste0("cluster_all_filtered_longest_8kmer_p", i, ".csv")
    ))
}
