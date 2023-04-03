rm(list = ls())

library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)

load_data <- function(path) {
  all_kmers <- fread(path) %>%
    distinct() %>%
    pivot_longer(
      cols = -c('ensembl_gene_id', 'gene_name'),
      names_to = 'kmer',
      values_to = 'occurrence'
    ) %>%
    mutate(ingroup = occurrence > 0) %>%
    select(-occurrence)
}

run_ks_tests <- function(all_kmers, params) {
  kmer_test_new_transcription <- all_kmers %>%
    left_join(params) %>%
    filter(!is.na(param_val)) %>%
    group_by(kmer, param_name) %>%
    filter(max(ingroup) == 1) %>%
    filter(min(ingroup) == 0) %>%
    summarise(
      # n = n(),
      ks_pval_gr = ks.test(param_val[ingroup], param_val[!ingroup], alternative = 'greater')$p.value,
      ks_pval_ts = ks.test(param_val[ingroup], param_val[!ingroup], alternative = 'two.sided')$p.value,
      ks_pval_ls = ks.test(param_val[ingroup], param_val[!ingroup], alternative = 'less')$p.value,
      tot_mean = mean(param_val, na.rm = T),
      kmer_mean = mean(param_val[ingroup], na.rm = T),
      bg_mean = mean(param_val[!ingroup], na.rm = T),
      tot_median = median(param_val, na.rm = T),
      kmer_median = median(param_val[ingroup], na.rm = T),
      bg_median = median(param_val[!ingroup], na.rm = T),
      tot_sd = sd(param_val, na.rm = T),
      kmer_sd = sd(param_val[ingroup], na.rm = T),
      bg_sd = sd(param_val[!ingroup], na.rm = T),
      kmer_size = sum(ingroup, na.rm = T),
      bg_size = sum(!ingroup, na.rm = T)
    )
}

get_ks_res <- function(path, params) {
  gc()
  data <- load_data(path)
  gc()
  return(run_ks_tests(data, params))
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) <= 3) {
  stop("Please supply path to folder with kmer tables, parameter table and output path", call. = FALSE)
}

if (file.exists(args[1])) {
  kmer.base <- args[1]
} else {
  stop(paste("folder", args[1], "doesn't exist"))
}

if (file.exists(args[2])) {
  params.path <- args[2]
} else {
  stop(paste("file", args[2], "doesn't exist"))
}

if (file.exists(args[3])) {
  out.base <- args[3]
} else {
  stop(paste("folder", args[3], "doesn't exist"))
}

params_df <- fread(params.path)

if (length(args) == 4) {
  filter.term <- args[4]

  params_df <- params_df %>%
    mutate(param_name = paste(category, param_name)) %>%
    filter(category == filter.term) %>%
    select(-category)
}

file_paths <- list.files(path = kmer.base)

t0 <- Sys.time()
ks_list <- lapply(file.path(kmer.base, file_paths), get_ks_res, params = params_df)
t1 <- Sys.time()

all_ks <- rbindlist(ks_list)

all_ks_3 <- all_ks %>%
  pivot_longer(
    cols = c(ks_pval_gr, ks_pval_ts, ks_pval_ls),
    names_to = 'alternative',
    values_to = "ks_pval"
  ) %>%
  mutate(
    alternative = factor(alternative,
                         levels = c('ks_pval_gr', 'ks_pval_ts', 'ks_pval_ls'),
                         labels = c('greater', 'two.sided', 'less')
    )
  )

all_ks_plot <- all_ks_3 %>%
  mutate(
    # ks_pval = ifelse(ks_pval == 0, 4.563133e-20, ks_pval),
    ks_pval_adj = p.adjust(ks_pval, "fdr"),
    pval_fdr_log = -log10(ks_pval_adj),
    mean_diff = (kmer_mean - bg_mean),
    effect_size = mean_diff / tot_sd,
  ) %>%
  select(kmer, param_name, alternative, pval_fdr_log, effect_size)

fwrite(all_ks, file.path(out.base, paste(filter.term, "ks_raw_with_stats.tsv", sep = '_')), sep = '\t')
# fwrite(all_ks_plot, file.path(out.base, paste(filter.term, "ks_adjusted_effect.tsv", sep = '_')), sep = '\t')

