rm(list = ls())

library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(viridis)
library(patchwork)
library(ggrepel)

es <- 0.2
esm <- 1.5
pv <- -log10(0.01)

# calculate density. https://slowkow.com/notes/ggplot2-color-by-density/:
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

correct_p_values <- function(all_ks) {
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
    filter(alternative != 'two.sided') %>%
    group_by(param_name) %>%
    mutate(ks_pval_adj = p.adjust(ks_pval, "fdr")) %>%
    ungroup() %>%
    mutate(
      pval_fdr_log = -log10(ks_pval_adj),
      mean_diff = (kmer_mean - bg_mean),
      effect_size = mean_diff / tot_sd,
    ) %>%
    select(kmer, param_name, alternative, pval_fdr_log, effect_size) #%>%
  all_ks_plot_max <- all_ks_plot %>%
    filter(alternative != 'two.sided') %>%
    group_by(kmer, param_name, effect_size) %>%
    summarise(pval_fdr_log = max(pval_fdr_log))

  return(all_ks_plot_max)
}

plot_parameter_volcano <- function(all_data, param) {
  filtered <- all_data %>%
    filter(param_name == param) %>%
    mutate(regulated =
             case_when(
               pval_fdr_log > pv & effect_size > es ~ 'ur',
               pval_fdr_log > pv & effect_size < -es ~ 'dr',
               TRUE ~ 'not'
             )
    ) %>%
    mutate(effect_size = ifelse(effect_size < -esm, -esm, effect_size),
           effect_size = ifelse(effect_size > esm, esm, effect_size))
  filtered$density <- get_density(filtered$effect_size, filtered$pval_fdr_log, h = c(1, 2), n = 300)

  ggplot(filtered, aes(effect_size, pval_fdr_log, label = gsub("T", "U", kmer))) +
    geom_point(alpha = 1, aes(color = density)) +
    geom_point(data = filtered %>% filter(regulated == 'ur'), color = 'blue') +
    geom_point(data = filtered %>% filter(regulated == 'dr'), color = 'red') +
    geom_text_repel(data = filtered %>% filter(regulated != 'not'), fontface ='bold', size = 8, max.overlaps = 50) +
    geom_hline(yintercept = pv, linetype = 'dashed') +
    geom_vline(xintercept = c(-es, es), linetype = 'dashed') +
    facet_wrap(param_name ~ ., scales = 'free_x') +
    scale_color_viridis(direction = 1) +
    theme_classic(base_size = 20) +
    xlim(-esm, esm) +
    xlab('effect size (SMD)') +
    ylab("-log10(p-value)") +
    theme(
      legend.position = "none"
    )
}


args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Please supply path to file with ks test results and output path", call. = FALSE)
}

if (file.exists(args[1])) {
  all.ks.path <- args[1]
} else {
  stop(paste("file", args[1], "doesn't exist"))
}

if (file.exists(args[2])) {
  out.base <- args[2]
} else {
  stop(paste("folder", args[2], "doesn't exist"))
}

all_ks <- fread(all.ks.path)
all_ks_plot_max <- correct_p_values(all_ks)

for (p in unique(all_ks_plot_max$param_name)) {
  plot_parameter_volcano(all_ks_plot_max, p)
  ggsave(file.path(out.base,
                   paste0(gsub(" |-", "_", p), "_volcano_one_sided.png")),
         height = 13,
         width = 15,
         dpi = 1000
  )
}

