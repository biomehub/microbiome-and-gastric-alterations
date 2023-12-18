# Load packages ----------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)
library(DESeq2)
library(corncob)
library(limma)
library(edgeR)  # BiocManager::install("edgeR")


# Basic setup -----------------------------------------------------
source("R/utils.R")
theme_set(theme_bw(base_size = 13))
outdir <- "output/bioinfo864"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# Load data ---------------------------------------------------------------------
.otu <- read_tsv("data/220326-195253_otu_table.tsv",
                 show_col_types = FALSE) %>% 
  mutate(sample = str_extract(sample, "\\d+")) %>% 
  column_to_rownames("sample") %>% 
  as.matrix()
.tax <- read_tsv("data/220326-195253_tax_table.tsv",
                 show_col_types = FALSE) %>% 
  select(-`...1`) %>% 
  column_to_rownames("otu") %>% 
  as.matrix()
.meta <- read_tsv("data/metadata.tsv",
                  show_col_types = FALSE) %>% 
  mutate(
    endoscopia = as.character(endoscopia),
    y = factor(
      ifelse(endoscopia == "normal", "normal", "alterada"),
      levels = c("normal", "alterada")
    ),
    endoscopia = case_when(
      endoscopia == "normal" ~ "normal",
      endoscopia == "doença péptica gastroduodenal" ~ "DPG",
      endoscopia == "esofagite erosiva" ~ "EE",
      endoscopia == "esofagite erosiva, doença péptica gastroduodenal" ~ "DPG+EE",
      TRUE ~ NA_character_
    ) %>% 
      factor(
        levels = c(
          "normal",
          "DPG",
          "EE",
          "DPG+EE"
        )
      ),
    coleta = factor(
      as.character(coleta),
      levels = c("coleta-1", "coleta-2")
    )
  ) %>% 
  column_to_rownames('amostra')

phylo <- phyloseq(
  otu_table(.otu, taxa_are_rows = FALSE),
  tax_table(.tax),
  sample_data(.meta)
)

zeroed_samples <- sample_names(phylo)[
  sample_sums(phylo) == 0
]
write_lines(zeroed_samples, "output/zeroed-samples.tsv")
phylo <- prune_samples(
  sample_sums(phylo) > 0, phylo
)
phylo <- prune_taxa(
  taxa_sums(phylo) > 0,
  phylo
)
phylo.prop <- transform_sample_counts(
  phylo, 
  \(x) x/sum(x)
)

included_samples <- sample_names(phylo)
write_lines(included_samples, "output/included-samples.tsv")

## colors ----

.cols <- list(
  endoscopia = RColorBrewer::brewer.pal(4, "Dark2") %>% 
    set_names(nm = levels(phylo@sam_data$endoscopia)),
  y = RColorBrewer::brewer.pal(3, "Set1")[1:2] %>% 
    set_names(nm = levels(phylo@sam_data$y)),
  coleta = RColorBrewer::brewer.pal(3, "Set2")[1:2] %>% 
    set_names(nm = levels(phylo@sam_data$coleta))
)


# Diversity analysis ----

## alpha-diversity ----

df_alpha <- left_join(
  estimate_richness(phylo, measures = c("Observed", "Shannon", "InvSimpson")) %>% 
    rownames_to_column("sample_id") %>% 
    mutate(sample_id = str_extract(sample_id, "\\d+")),
  data.frame(sample_data(phylo)) %>% 
    rownames_to_column("sample_id"),
  by = "sample_id"
) %>% 
  dplyr::rename(
    Richness := Observed
  ) %>% 
  pivot_longer(cols = c(Richness, Shannon, InvSimpson))

for (.metric in c("Richness", "Shannon", "InvSimpson")) {
  
  .df_alpha <- df_alpha %>% 
    filter(name == .metric)
  
  pvals <- tribble(
    ~name, ~pval,
    "coleta", wilcox.test(value ~ coleta, data = .df_alpha)$p.value,
    "y", wilcox.test(value ~ y, data = .df_alpha)$p.value,
    "endoscopia", kruskal.test(value ~ endoscopia, data = .df_alpha)$p.value
  ) %>% 
    mutate(padj = p.adjust(pval)) %>% 
    group_split(name) %>% 
    set_names(map(., ~ .x$name))
  
  p1 <- .df_alpha %>% 
    ggplot(aes(endoscopia, value)) +
    geom_quasirandom(aes(fill = endoscopia),
                     shape = 21, size = 3, color = "gray20") +
    geom_boxplot(alpha = 0, lwd = 1) +
    scale_fill_manual(values = .cols$endoscopia) +
    labs(
      x = "Endoscopia", y = .metric,
      subtitle = paste0(
        "Kruskal-Wallis Rank Sum Test p=", 
        round(pvals$endoscopia$padj, 3)
      )
    ) +
    guides(fill = "none")
  
  p2 <- .df_alpha %>% 
    ggplot(aes(y, value)) +
    geom_quasirandom(aes(fill = y),
                     shape = 21, size = 3, color = "gray20") +
    geom_boxplot(alpha = 0, lwd = 1) +
    scale_fill_manual(values = .cols$y) +
    labs(
      x = "Endoscopia", y = .metric,
      subtitle = paste0(
        "Wilcoxon Rank Sum Test p=", 
        round(pvals$y$padj, 3)
      )
    ) +
    guides(fill = "none")
  
  p3 <- .df_alpha %>% 
    ggplot(aes(coleta, value)) +
    geom_quasirandom(aes(fill = coleta),
                     shape = 21, size = 3, color = "gray20") +
    geom_boxplot(alpha = 0, lwd = 1) +
    scale_fill_manual(values = .cols$coleta) +
    labs(
      x = "Coleta", y = .metric,
      subtitle = paste0(
        "Wilcoxon Rank Sum Test p=", 
        round(pvals$coleta$padj, 3)
      )
    ) +
    guides(fill = "none")
  
  .metric <- str_to_lower(.metric)
  
  ggsave(
    f("{outdir}/alpha-diversity/alpha-diversity-{.metric}-endoscopy.png"),
    p1, 
    width = 6, height = 4.5
  )
  ggsave(
    f("{outdir}/alpha-diversity/alpha-diversity-{.metric}-endoscopy-binary.png"),
    p2, 
    width = 6, height = 4.5
  )
  ggsave(
    f("{outdir}/alpha-diversity/alpha-diversity-{.metric}-coleta.png"),
    p3, 
    width = 6, height = 4.5
  )
}



## beta-diversity ----

.distance <- phyloseq::distance(
  phylo.prop, method = 'bray'
)

.ord <- ordinate(phylo.prop, method = "PCoA", distance = .distance)

pcoa1 <- plot_ordination(phylo.prop, .ord, axes = c(1,2))
pcoa2 <- plot_ordination(phylo.prop, .ord, axes = c(1,3))

get_pmvn <- \(x, .names) x %>% 
  as_tibble(rownames = "term") %>% 
  filter(term %in% .names) %>% 
  dplyr::select(term, R2, pval := `Pr(>F)`)

pmnv <- vegan::adonis2(
  .distance ~ endoscopia + coleta, 
  data = data.frame(sample_data(phylo.prop)),
  by = "margin",
  permutations = 3000
) %>% 
  get_pmvn(.names = c("coleta", "endoscopia"))

pmnv_binary <- vegan::adonis2(
  .distance ~ y + coleta, 
  data = data.frame(sample_data(phylo.prop)),
  by = "margin",
  permutations = 3000
) %>% 
  get_pmvn(.names = c("y"))

pmnv <- bind_rows(
  pmnv, pmnv_binary
) %>% 
  mutate(padj = p.adjust(pval, "BH")) %>% 
  group_split(term) %>% 
  set_names(map(., ~ .x$term))

### Endoscopia ----

p1 <- pcoa1$data %>% 
  ggplot(aes(Axis.1, Axis.2)) +
  geom_point(aes(fill = endoscopia),
             shape = 21, size = 3, color = "gray20") +
  scale_fill_manual(values = .cols$endoscopia) +
  labs(
    x = pcoa1$labels$x,
    y = pcoa1$labels$y,
    fill = NULL
  )
p2 <- pcoa2$data %>% 
  ggplot(aes(Axis.1, Axis.3)) +
  geom_point(aes(fill = endoscopia),
             shape = 21, size = 3, color = "gray20") +
  scale_fill_manual(values = .cols$endoscopia) +
  labs(
    x = pcoa2$labels$x,
    y = pcoa2$labels$y,
    fill = NULL
  )


.p <- (p1 | p2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    paste0(
      "PERMANOVA p=",
      round(pmnv$endoscopia$padj, 3),
      ", R2=",
      round(pmnv$endoscopia$R2, 3)
    )
  )

ggsave(
  f("{outdir}/beta-diversity/beta-diversity-endoscopy.png"),
  .p, 
  width = 12, height = 4.5,bg = 'white'
)


### Endoscopia (binário) ----


p3 <- pcoa1$data %>% 
  ggplot(aes(Axis.1, Axis.2)) +
  geom_point(aes(fill = y),
             shape = 21, size = 3, color = "gray20") +
  scale_fill_manual(values = .cols$y) +
  labs(
    x = pcoa1$labels$x,
    y = pcoa1$labels$y,
    fill = NULL
  )
p4 <- pcoa2$data %>% 
  ggplot(aes(Axis.1, Axis.3)) +
  geom_point(aes(fill = y),
             shape = 21, size = 3, color = "gray20") +
  scale_fill_manual(values = .cols$y) +
  labs(
    x = p2$labels$x,
    y = p2$labels$y,
    fill = NULL
  )

.p <- (p3 | p4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    paste0(
      "PERMANOVA p=",
      round(pmnv$y$padj, 3),
      ", R2=",
      round(pmnv$y$R2, 3)
    )
  )

ggsave(
  f("{outdir}/beta-diversity/beta-diversity-endoscopy-binary.png"),
  .p, 
  width = 12, height = 4.5,bg = 'white'
)

### Coleta ----

p5 <- pcoa1$data %>% 
  ggplot(aes(Axis.1, Axis.2)) +
  geom_point(aes(fill = coleta),
             shape = 21, size = 3, color = "gray20") +
  scale_fill_manual(values = .cols$coleta) +
  labs(
    x = pcoa1$labels$x,
    y = pcoa1$labels$y,
    fill = NULL
  )
p6 <- pcoa2$data %>% 
  ggplot(aes(Axis.1, Axis.3)) +
  geom_point(aes(fill = coleta),
             shape = 21, size = 3, color = "gray20") +
  scale_fill_manual(values = .cols$coleta) +
  labs(
    x = p2$labels$x,
    y = p2$labels$y,
    fill = NULL
  )


.p <- (p5 | p6) +
  plot_layout(guides = "collect") +
  plot_annotation(
    paste0(
      "PERMANOVA p=",
      round(pmnv$coleta$padj, 3),
      ", R2=",
      round(pmnv$coleta$R2, 3)
    )
  )

ggsave(
  f("{outdir}/beta-diversity/beta-diversity-coleta.png"),
  .p, 
  width = 12, height = 4.5,bg = 'white'
)


# Differential Abundance ----

## DESeq2 ----

rank_list <- c(
  "lowest",
  "species",
  "genus",
  "family",
  "phylum"
) %>% 
  set_names(.)

min_prevalence <- 0.05

da_deseq2 <- map_df(rank_list, \(.rank) {
  print(str_glue("\n\nRank: {.rank}\n"))
  phylo.rank <- tax_glom2(phylo, .rank = .rank)
  phylo.rank <- prune_taxa(
    # taxa must be present in at least `min_prevalence` of samples
    taxa_prevalence(phylo.rank) >= min_prevalence,
    phylo.rank
  )
  phylo.rank@sam_data$coleta <- factor(
    str_replace(phylo.rank@sam_data$coleta, "-", "_"),
    levels = c("coleta_1", "coleta_2")
  )
  phylo.rank@sam_data$endoscopia <- factor(
    str_replace(phylo.rank@sam_data$endoscopia, "\\+", "_"),
    levels = c("normal", "DPG", "EE", "DPG_EE")
  )
  
  lrt_pvals <- phylo.rank %>%
    phyloseq_to_deseq2(~ coleta + endoscopia) %>%
    DESeq(sfType = 'poscounts', test = "LRT", reduced = ~ coleta, quiet = T) %>%
    results(tidy = T) %>% 
    select(bac := row, 
           lrt_pvalue := pvalue, 
           lrt_padj := padj)
  poshoc_ds <- phylo.rank %>%
    phyloseq_to_deseq2(~ endoscopia + coleta) %>%
    DESeq(sfType = 'poscounts', betaPrior = TRUE, quiet = TRUE)
  beta_estimates <- inner_join(
    results(poshoc_ds, 
            contrast = c("endoscopia", "DPG", "normal"), tidy = T) %>% 
      select(bac := row, 
             DPG_vs_normal := log2FoldChange, 
             DPG_vs_normal_SE := lfcSE, 
             DPG_vs_normal_pvalue := pvalue, 
             DPG_vs_normal_padj := padj),
    results(poshoc_ds, 
            contrast = c("endoscopia", "EE", "normal"), tidy = T) %>% 
      select(bac := row, 
             EE_vs_normal := log2FoldChange, 
             EE_vs_normal_SE := lfcSE, 
             EE_vs_normal_pvalue := pvalue, 
             EE_vs_normal_padj := padj),
    by = "bac"
  ) %>% 
    inner_join(
      results(poshoc_ds, 
              contrast = c("endoscopia", "DPG_EE", "normal"), tidy = T) %>% 
        select(bac := row, 
               DPG_EE_vs_normal := log2FoldChange, 
               DPG_EE_vs_normal_SE := lfcSE, 
               DPG_EE_vs_normal_pvalue := pvalue, 
               DPG_EE_vs_normal_padj := padj),
      by = "bac"
    )
  
  res <- inner_join(
    lrt_pvals,
    beta_estimates,
    by = "bac"
  )
  
  return(res)
}, .id = "tax_rank") %>% 
  mutate(tool = "DESeq2")

### no differential abundant taxa:
### (decreasing order of "correctness")

da_deseq2 %>% 
  filter(lrt_padj < 0.1) %>% 
  nrow()
da_deseq2 %>% 
  filter(
    DPG_EE_vs_normal_padj < .1 |
      DPG_vs_normal_padj < .1 |
      EE_vs_normal_padj < .1
  )

# just so we have some fig from this

.p <- da_deseq2 %>% 
  select(ends_with('vs_normal')) %>%
  pivot_longer(
    everything()
  ) %>% 
  mutate(
    name = case_when(
      str_detect(name, "DPG_EE_vs") ~ "DPG+EE",
      str_detect(name, "DPG_vs") ~ "DPG",
      str_detect(name, "EE_vs") ~ "EE"
    ) %>% 
      factor(
        levels = c("DPG", "EE", "DPG+EE")
      )
  ) %>% 
  ggplot(aes(x = value, y = name, fill = name)) +
  ggridges::geom_density_ridges(
    jittered_points = TRUE,
    position = ggridges::position_points_jitter(width = 0, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
  ) +
  scale_fill_manual(values = .cols$endoscopia) +
  guides(fill = "none") +
  theme(
    plot.caption.position = "plot",
    plot.caption = element_text(
      size = 8, hjust = 0
    ),
    plot.subtitle = element_text(
      size = 10, hjust = 0
    )
  ) +
  labs(
    x = "log2FC", y = NULL,
    caption = "Comparisons against control participants.\nLog2FC for lowest taxonomic rank identified.",
    subtitle = "No effects were considered statistically significant at FDR < 10%"
  ) +
  xlim(-3.5, 3.5)

ggsave(
  f("output/bioinfo864/differential-abundance/diff-abund.png"),
  .p, width = 6, height = 5
)
## corcob ----


da_cc <- map_df(rank_list, \(.rank) {
  
  print(str_glue("\n\nRank: {.rank}\n"))
  phylo.rank <- tax_glom2(phylo, .rank = .rank)
  phylo.rank <- prune_taxa(
    # taxa must be present in at least `min_prevalence` of samples
    taxa_prevalence(phylo.rank) >= min_prevalence,
    phylo.rank
  )
  phylo.rank@sam_data$coleta <- factor(
    str_replace(phylo.rank@sam_data$coleta, "-", "_"),
    levels = c("coleta_1", "coleta_2")
  )
  phylo.rank@sam_data$endoscopia <- factor(
    str_replace(phylo.rank@sam_data$endoscopia, "\\+", "_"),
    levels = c("normal", "DPG", "EE", "DPG_EE")
  )
  
  .contrasts <- list(
    "endoscopiaDPG",
    "endoscopiaEE",
    "endoscopiaDPG_EE"
  ) %>% set_names(.)
  d <- contrastsTest(formula = ~ coleta + endoscopia,
                     phi.formula = ~ coleta,
                     contrasts_DA = .contrasts,
                     data = phylo.rank,
                     test = "Wald", 
                     fdr_cutoff = 1)
  res <- map_df(seq_along(.contrasts), ~{
    tibble(
      bac = taxa_names(d$data),
      term = .contrasts[[.x]],
      padj = d$p_fdr[, .x]
    )
  })

  return(res)
}, .id = "tax_rank") %>% 
  mutate(tool = "corncob")

### no differential abundant taxa:
### (decreasing order of "correctness")

sigs_corncob <- da_cc %>% 
  filter(padj < 0.1)
sigs_corncob

### probably just noise
phylo.prop %>% 
  tax_glom2("genus") %>% 
  get_taxon_couns('Mycoplasma') %>% 
  ggplot(aes(endoscopia, value)) +
  geom_boxplot() +
  geom_quasirandom() +
  stat_compare_means()

## limma + voom ----

da_lv <- map_df(rank_list, \(.rank) {
  print(str_glue("\n\nRank: {.rank}\n"))
  phylo.rank <- tax_glom2(phylo, .rank = .rank)
  phylo.rank <- prune_taxa(
    # taxa must be present in at least `min_prevalence` of samples
    taxa_prevalence(phylo.rank) >= min_prevalence,
    phylo.rank
  )
  phylo.rank@sam_data$coleta <- factor(
    str_replace(phylo.rank@sam_data$coleta, "-", "_"),
    levels = c("coleta_1", "coleta_2")
  )
  phylo.rank@sam_data$endoscopia <- factor(
    str_replace(phylo.rank@sam_data$endoscopia, "\\+", "_"),
    levels = c("normal", "DPG", "EE", "DPG_EE")
  )
  
  counts_table <- t(phylo.rank@otu_table@.Data)
  rownames(counts_table) <- str_replace(rownames(counts_table), " ", "_")
  
  dge_list <- edgeR::DGEList(counts_table) %>% 
    edgeR::calcNormFactors(method="TMM")
  
  design_matrix <- model.matrix(~ coleta + endoscopia, 
                                data.frame(phylo.rank@sam_data))
  voomvoom <- voom(dge_list, design_matrix, plot = FALSE)
  fit <- lmFit(voomvoom, design_matrix) %>% 
    eBayes()
  res <- c(
    "endoscopiaDPG",
    "endoscopiaEE",
    "endoscopiaDPG_EE"
  ) %>% 
    set_names(.) %>% 
    map_df(
      ~{
        fit %>% 
          topTable(
            coef = .x, 
            number = Inf, 
            adjust.method = "BH", 
            confint = TRUE
          ) %>% 
          as_tibble(rownames = "bac") %>% 
          mutate(
            term = .x
          ) %>% 
          dplyr::select(
            bac, 
            term,
            estimate := logFC,
            .lower := CI.L,
            .upper := CI.R,
            pvalue := P.Value,
            padj := adj.P.Val
          )
      }
    )
  
  return(res)
}, .id = "tax_rank") %>% 
  mutate(tool = "limma+voom")

sigs_lv <- da_lv %>% 
  filter(padj < 0.1)

## diff presence ----

dp_res <- map_df(rank_list, \(.rank) {
  print(str_glue("\n\nRank: {.rank}\n"))
  phylo.rank <- tax_glom2(phylo, .rank = .rank)
  phylo.rank <- prune_taxa(
    # taxa must be present in at least `min_prevalence` of samples
    taxa_prevalence(phylo.rank) >= min_prevalence,
    phylo.rank
  )
  phylo.rank@sam_data$coleta <- factor(
    str_replace(phylo.rank@sam_data$coleta, "-", "_"),
    levels = c("coleta_1", "coleta_2")
  )
  phylo.rank@sam_data$endoscopia <- factor(
    str_replace(phylo.rank@sam_data$endoscopia, "\\+", "_"),
    levels = c("normal", "DPG", "EE", "DPG_EE")
  )
  taxa_list <- taxa_names(phylo.rank) %>% 
    set_names(.)
  rank_res <- map_df(taxa_list, \(.taxname) {
    .df <- get_taxon_couns(.phyloseq = phylo.rank, taxon = .taxname) %>% 
      mutate(value = as.numeric(value > 0))
    fit <- logistf::logistf(
      value ~ coleta + endoscopia, data = .df
    )
    summary_logistf(fit, .print = FALSE) %>% 
      as_tibble(rownames = "term") %>% 
      filter(str_detect(term, "endoscopia")) %>% 
      dplyr::select(
        term, estimate := coef,
        .lower := `lower 0.95`,
        .upper := `upper 0.95`,
        pval := p
      )
  }, .id = "bac") %>% 
    mutate(
      padj = p.adjust(pval, method = "BH")
    )
}, .id = "tax_rank")

sigs_dp <- dp_res %>% 
  filter(padj < .1)

sigs_freq <- list(sigs_corncob, sigs_dp, sigs_lv) %>% 
  map(
    ~ {
      .x %>% 
        mutate(x = paste0(bac, "-", term, "-", tax_rank)) %>% 
        pull(x)
    }
  ) %>% 
  unlist() %>% 
  table() %>% 
  sort()
actually_sigs <- sigs_freq[sigs_freq > 1]
# they are all the same ASVs
as.data.frame(tax_table(phylo.prop)) %>% 
  filter(genus == "Mycoplasma")
df <- get_taxon_couns(phylo.prop, "Mycoplasma", "genus")


prev_data <- df %>% 
  group_by(endoscopia) %>% 
  summarise(
    x = binom::binom.wilson(sum(value>0), n())
  ) %>%
  unnest(x) %>% 
  mutate(
    .lab = paste0(
      "Prev. ",
      scales::percent(mean, accuracy = 1),
      "\n[",
      scales::percent(lower, accuracy = 1),
      " \u2012 ",
      scales::percent(upper, accuracy = 1),
      "]"
    ),
    value = 0.015
  ) %>% 
  dplyr::select(endoscopia, value, .lab)

df %>% 
  ggplot(aes(endoscopia, value)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(
    color = "gray40"
  ) +
  stat_summary(fun.data = mean_ci, col='red') +
  geom_text(
    data = prev_data,
    aes(label = .lab)
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    y = "Relative abundance",
    x = NULL,
    subtitle = "Mycoplasma (genus)"
  )

# Overall Variation ----

phylo_prop_list <- map(rank_list, ~{
  phylo.prop %>% 
    tax_glom2(.x)
})

for (.rank in rank_list) {
  
  phylo.prop.rank <- phylo_prop_list[[.rank]]
  
  bacs_ix <- colMeans(otu_table(phylo.prop.rank) > 0) > min_prevalence
  counts_matrix <- phylo.prop.rank@otu_table@.Data[,bacs_ix]
  counts_matrix <- counts_matrix[, order(colMeans(counts_matrix), decreasing = T)]
  
  file_path <- f('output/bioinfo864/relative-abundances/{.rank}.png')
  png(
    file_path,
    res = 600, 
    width = ifelse(sum(bacs_ix) > 15, 1e4, 4500), 
    height = ifelse(sum(bacs_ix) > 15, 9e3, 3500)
  )
  pheatmap::pheatmap(
    100*t(counts_matrix),
    annotation_col = phylo.prop.rank@sam_data %>% 
      as_tibble(rownames = "sample_id") %>% 
      mutate(alteração = ifelse(endoscopia == "normal", "não", "sim"),
             coleta = str_replace(coleta, "-", " ")) %>% 
      column_to_rownames("sample_id") %>% 
      dplyr::select(alteração, endoscopia, coleta) ,
    show_colnames = F,
    cluster_rows = FALSE,
    scale = 'none',
    color = viridis::viridis(100), 
    legend_breaks = seq(0, 100, 25),
    legend_labels = paste0(seq(0, 100, 25), "%")
  )
  dev.off()
}
