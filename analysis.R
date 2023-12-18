# Load packages ----------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)
library(DESeq2)
library(corncob)
library(limma)
library(MicrobiomeStat) # install.packages('MicrobiomeStat')
library(mice)  # install.packages('mice')



# Basic setup -----------------------------------------------------
source("R/utils.R")
theme_set(theme_bw(base_size = 13))
outdir <- "output/analysis"
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
      ifelse(endoscopia == "normal", "no", "yes"),
      levels = c("no", "yes")
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
      ifelse(coleta == "coleta-1", "Center 1", "Center 2"),
      levels = c("Center 1", "Center 2")
    ),
    sintomas = factor(
      ifelse(sintomas == "no", "no", "yes"),
      levels = c("no", "yes")
    ),
    ibp = factor(
      as.character(ibp),
      levels = c("no", "yes")
    ),
    tabagismo = factor(
      as.character(tabagismo),
      levels = c("no", "yes")
    ),
    etilismo = factor(
      as.character(etilismo),
      levels = c("no", "yes")
    )
  ) %>%
  dplyr::select(-atb) %>% # too few positives
  column_to_rownames('amostra')
.meta$complete_case <- apply(.meta, 1, \(j) !any(is.na(j)))

# multiple imputation of metadata values

pick_mode <- function(z) {
  as.character(z[which.max(table(na.omit(z)))])
}
imp <- .meta %>% 
  dplyr::select(-endoscopia, -y, -complete_case) %>% 
  mice::mice(maxit = 10, m = 10, seed = 12345, printFlag = FALSE) %>% 
  mice::complete('long', inc = TRUE) %>% 
  filter(.imp > 0) %>% 
  group_by(.id) %>% 
  summarise(
    ibp = pick_mode(ibp),
    tabagismo = pick_mode(tabagismo),
    etilismo = pick_mode(etilismo),
    coleta = pick_mode(coleta),
    sintomas = pick_mode(sintomas)
  )

.meta$ibp <- factor(imp$ibp, levels = c("no", "yes"))
.meta$tabagismo <- factor(imp$tabagismo, levels = c("no", "yes"))
.meta$etilismo <- factor(imp$etilismo, levels = c("no", "yes"))
.meta$sintomas <- factor(imp$sintomas, levels = c("no", "yes"))

any(is.na(.meta))

phylo <- phyloseq(
  otu_table(.otu, taxa_are_rows = FALSE),
  tax_table(.tax),
  sample_data(.meta)
)

zeroed_samples <- sample_names(phylo)[
  sample_sums(phylo) == 0
]
write_lines(zeroed_samples, "output/analysis/zeroed-samples.tsv")
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
write_lines(included_samples, "output/analysis/included-samples.tsv")

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
    # same results with lm(log(value) ~ coleta + y (or endoscopia) + ibp + ..., data = .df_value)
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
      x = "Endoscopy", y = .metric,
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
      x = "Any endoscopy alteration", y = .metric,
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
      x = NULL, y = .metric,
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


.p1 <- (p1 | p2) +
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
  .p1, 
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
    fill = "Any endoscopy\nalteration"
  )
p4 <- pcoa2$data %>% 
  ggplot(aes(Axis.1, Axis.3)) +
  geom_point(aes(fill = y),
             shape = 21, size = 3, color = "gray20") +
  scale_fill_manual(values = .cols$y) +
  labs(
    x = p2$labels$x,
    y = p2$labels$y,
    fill = "Any endoscopy\nalteration"
  )

.p2 <- (p3 | p4) +
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
  .p2, 
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


.p3 <- (p5 | p6) +
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
  .p3, 
  width = 12, height = 4.5, bg = 'white'
)

### all together ----

ggsave(
  f("{outdir}/beta-diversity/beta-diversity-all.png"),
  ggarrange(.p2, .p1, .p3, nrow = 3, ncol = 1),
  width = 12, height = 4.5*2.5, bg = 'white'
)


# Differential Abundance ----

## DESeq2 ----

rank_list <- c(
  "lowest",    # ASV level
  "species",
  "genus",
  "family",
  "phylum"
) %>% 
  set_names(.)

phylo_list <- map(rank_list, ~ tax_glom2(phylo, .rank = .x))

min_prevalence <- 0.05

da_deseq2 <- map_df(rank_list, \(.rank) {
  print(str_glue("\n\nRank: {.rank}\n"))
  phylo.rank <- phylo_list[[.rank]]
  phylo.rank <- prune_taxa(
    # taxa must be present in at least `min_prevalence` of samples
    taxa_prevalence(phylo.rank) >= min_prevalence,
    phylo.rank
  )
  phylo.rank@sam_data$coleta <- factor(
    str_replace(phylo.rank@sam_data$coleta, " ", "_"),
    levels = c("Center_1", "Center_2")
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
  f("output/analysis/differential-abundance/diff-abund.png"),
  .p, width = 6, height = 5
)
## corcob ----


da_cc <- map_df(rank_list, \(.rank) {
  
  print(str_glue("\n\nRank: {.rank}\n"))
  phylo.rank <- phylo_list[[.rank]]
  phylo.rank <- prune_taxa(
    # taxa must be present in at least `min_prevalence` of samples
    taxa_prevalence(phylo.rank) >= min_prevalence,
    phylo.rank
  )
  phylo.rank@sam_data$coleta <- factor(
    str_replace(phylo.rank@sam_data$coleta, " ", "_"),
    levels = c("Center_1", "Center_2")
  )
  phylo.rank@sam_data$endoscopia <- factor(
    str_replace(phylo.rank@sam_data$endoscopia, "\\+", "_"),
    levels = c("normal", "DPG", "EE", "DPG_EE")
  )
  
  d_lrt <- corncob::differentialTest(
    formula = ~ coleta + endoscopia,
    formula_null = ~ coleta,
    phi.formula = ~ coleta,
    phi.formula_null = ~ coleta,
    data = phylo.rank,
    test = "LRT", 
    fdr_cutoff = 1
  )
  
  .level <- ifelse(.rank == "lowest", "lowest_taxonomy", .rank)
  estimates <- corncob:::plot.differentialTest(d_lrt, level = .level, data_only = T) %>%
    mutate(
      contrast = paste0(
        str_extract(variable, "DPG_EE|DPG|EE"), 
        "_vs_normal"
      )
    ) %>% 
    select(bac := taxa, x, contrast) %>% 
    pivot_wider(
      names_from = contrast,
      values_from = x
    )
  
  res_lrt <- tibble(
    bac = names(d_lrt$p_fdr),
    lrt_pvalue = d_lrt$p,
    lrt_padj = d_lrt$p_fdr
  ) %>% 
    inner_join(
      estimates,
      by = 'bac'
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
      pvalue = d$p[, .x],
      padj = d$p_fdr[, .x]
    )
  }) %>%
    pivot_wider(
      names_from = term,
      values_from = c(pvalue, padj),
      names_repair = function(x) ifelse(
        x == "bac", x,
        paste0(
          str_extract(x, "DPG_EE|DPG|EE"),
          "_vs_normal_",
          str_extract(x, "pvalue|padj")
        )
      )
    )
  
  res <- inner_join(res_lrt, res, by = 'bac')

  return(res)
}, .id = "tax_rank") %>% 
  mutate(tool = "corncob")


da_cc %>% 
  filter(lrt_padj < 0.1)

## limma + voom ----

da_lv <- map_df(rank_list, \(.rank) {
  print(str_glue("\n\nRank: {.rank}\n"))
  phylo.rank <- phylo_list[[.rank]]
  phylo.rank <- prune_taxa(
    # taxa must be present in at least `min_prevalence` of samples
    taxa_prevalence(phylo.rank) >= min_prevalence,
    phylo.rank
  )
  phylo.rank@sam_data$coleta <- factor(
    str_replace(phylo.rank@sam_data$coleta, " ", "_"),
    levels = c("Center_1", "Center_2")
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
            pvalue := P.Value,
            padj := adj.P.Val
          )
      }
    ) %>% 
    pivot_wider(
      names_from = term,
      values_from = c(estimate, pvalue, padj),
      names_repair = function(x) ifelse(
        x == "bac",
        x,
        paste0(
          str_extract(x, "DPG_EE|DPG|EE"),
          "_vs_normal",
          ifelse(
            str_detect(x, "estimate"),
            "",
            paste0("_", str_extract(x, "pvalue|padj"))
          )
        )
      )
    )
  res_f <- topTable(fit, coef = c(3, 4, 5), number = Inf) %>% 
    as_tibble(rownames = "bac") %>% 
    select(bac, lrt_pvalue := P.Value, lrt_padj := adj.P.Val)
  
  res <- inner_join(res, res_f, by = "bac")
  
  return(res)
}, .id = "tax_rank") %>% 
  mutate(tool = "limma+voom")

da_lv %>% 
  filter(lrt_padj < 0.1)

## diff presence ----

dp_res <- map_df(rank_list, \(.rank) {
  print(str_glue("\n\nRank: {.rank}\n"))
  phylo.rank <- phylo_list[[.rank]]
  phylo.rank <- prune_taxa(
    # taxa must be present in at least `min_prevalence` of samples
    taxa_prevalence(phylo.rank) >= min_prevalence,
    phylo.rank
  )
  phylo.rank@sam_data$coleta <- factor(
    str_replace(phylo.rank@sam_data$coleta, " ", "_"),
    levels = c("Center_1", "Center_2")
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
    fitnull <- logistf::logistf(
      value ~ coleta, data = .df,
      control = logistf::logistf.control(maxit = 100, maxstep = 10)
    )
    fit <- logistf::logistf(
      value ~ coleta + endoscopia, data = .df,
      control = logistf::logistf.control(maxit = 100, maxstep = 10)
    )
    
    summary_logistf(fit) %>% 
      as_tibble(rownames = "term") %>% 
      filter(str_detect(term, "endoscopia")) %>% 
      dplyr::select(
        term, estimate := coef,
        pval := p
      ) %>% 
      pivot_wider(
        names_from = term,
        values_from = c(estimate, pval),
        names_repair = function(x) paste0(
          str_extract(x, "DPG_EE|DPG|EE"),
          "_vs_normal",
          ifelse(
            str_detect(x, "estimate"),
            "",
            paste0("_", str_extract(x, "pval"))
          )
        )
      ) %>% 
      mutate(
        lrt_pval = anova(fit, fitnull, method = "PLR")$pval[[1]]
      )
  }, .id = "bac") %>%
    rename_with(~str_replace(., 'pval', 'pvalue')) %>% 
    mutate(
      lrt_padj = p.adjust(lrt_pvalue, method = "BH"),
      DPG_vs_normal_padj = p.adjust(DPG_vs_normal_pvalue, method = "BH"),
      DPG_EE_vs_normal_padj = p.adjust(DPG_EE_vs_normal_pvalue, method = "BH"),
      EE_vs_normal_padj = p.adjust(EE_vs_normal_pvalue, method = "BH")
    )
}, .id = "tax_rank") %>% 
  mutate(tool = "Firth")

dp_res %>% 
  filter(lrt_padj < .1)

## LinDA / MicrobiomeStat ----

da_linda <- map_df(rank_list, \(.rank) {
  print(str_glue("\n\nRank: {.rank}\n"))
  phylo.rank <- phylo_list[[.rank]]
  phylo.rank <- prune_taxa(
    # taxa must be present in at least `min_prevalence` of samples
    taxa_prevalence(phylo.rank) >= min_prevalence,
    phylo.rank
  )
  phylo.rank@sam_data$coleta <- factor(
    str_replace(phylo.rank@sam_data$coleta, " ", "_"),
    levels = c("Center_1", "Center_2")
  )
  phylo.rank@sam_data$endoscopia <- factor(
    str_replace(phylo.rank@sam_data$endoscopia, "\\+", "_"),
    levels = c("normal", "DPG", "EE", "DPG_EE")
  )

  .meta_data <- data.frame(phylo.rank@sam_data)
  .otu_data <- t(phylo.rank@otu_table@.Data)[, rownames(.meta_data)]

  # no overall F test - annoying
  res_overall <- MicrobiomeStat::linda(
    feature.dat = .otu_data, 
    meta.dat = .meta_data, 
    formula = "~ coleta + y",
    feature.dat.type = 'count',
    prev.filter = min_prevalence,
    zero.handling = "imputation",
    n.cores = 4
  )$output$yyes %>%
    rownames_to_column("bac") %>% 
    select(bac, lrt_pvalue := pvalue, lrt_padj := padj)
  
  res <- MicrobiomeStat::linda(
    feature.dat = .otu_data, 
    meta.dat = .meta_data, 
    formula = "~ coleta + endoscopia",
    feature.dat.type = 'count',
    prev.filter = min_prevalence,
    zero.handling = "imputation",
    n.cores = 4
  )$output %>% 
    map_df(~ rownames_to_column(.x, "bac"), .id = "term") %>%
    filter(str_detect(term, "endoscopia")) %>% 
    select(term, bac, estimate := log2FoldChange, pvalue, padj) %>% 
    pivot_wider(
      names_from = term,
      values_from = c(estimate, pvalue, padj),
      names_repair = function(x) ifelse(
        x == "bac",
        x,
        paste0(
          str_extract(x, "DPG_EE|DPG|EE"),
          "_vs_normal",
          ifelse(
            str_detect(x, "estimate"),
            "",
            paste0("_", str_extract(x, "pvalue|padj"))
          )
        )
      )
    )

  res <- inner_join(res, res_overall, by = "bac")
  
  return(res)
}, .id = "tax_rank") %>% 
  mutate(tool = "linda")

da_linda %>% 
  filter(lrt_padj < 0.1)

# combine DA results ----

all_res <- bind_rows(da_deseq2, da_cc, da_lv, dp_res, da_linda) %>% 
  mutate(
    bac = str_replace_all(bac, "_", " ")
  ) %>% 
  select(-contains("_SE"))

all_res %>% 
  filter(lrt_padj < .1) %>% 
  group_by(tax_rank, bac) %>% 
  summarise(
    n = n(),
    tools = paste0(unique(tool), collapse = "; ")
  ) %>%
  arrange(desc(n))
all_res %>%
  group_by(tool) %>% 
  mutate(lrt_padj2 = p.adjust(lrt_pvalue, method = 'BH')) %>% 
  select(bac, tax_rank, tool, contains("lrt")) %>% 
  filter(lrt_padj2<.1)

df <- all_res %>% 
  select(-contains("_SE")) %>%
  mutate(tool = ifelse(tool == "Firth", "logistf", tool)) %>% 
  filter(tax_rank != "lowest")

## volcano ----
p1 <- df %>%
  mutate(
    sig = ifelse(
      DPG_vs_normal_padj< 0.1 & lrt_padj < .1,
      "Significant at FDR 10%",
      "Not significant at FDR 10%"
    )
  ) %>% 
  na.omit() %>% 
  ggplot(aes(DPG_vs_normal, -log10(DPG_vs_normal_pvalue))) +
  ggrepel::geom_text_repel(
    aes(label = str_replace_all(bac, "_", " ")),
    data = . %>% filter(!str_detect(sig, "Not"))
  ) +
  geom_point(aes(color = sig)) +
  labs(
    color = NULL,
    x = NULL,
    y = "-log10(p-value)",
    title = "DPG vs normal"
  ) +
  facet_wrap(~tool, ncol = 5) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 5)) +
  guides(color = "none")
p2 <- df %>%
  mutate(
    sig = ifelse(
      EE_vs_normal_padj< 0.1 & lrt_padj < .1,
      "Significant at FDR 10%",
      "Not significant at FDR 10%"
    )
  ) %>% 
  na.omit() %>% 
  ggplot(aes(EE_vs_normal, -log10(EE_vs_normal_pvalue))) +
  ggrepel::geom_text_repel(
    aes(label = str_replace_all(bac, "_", " ")),
    data = . %>% filter(!str_detect(sig, "Not")),
    max.time = 10, max.iter = 3e4
  ) +
  geom_point(aes(color = sig)) +
  labs(
    color = NULL,
    x = NULL,
    y = "-log10(p-value)",
    title = "EE vs normal"
  ) +
  facet_wrap(~tool, ncol = 5) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 5))
p3 <- df %>%
  mutate(
    sig = ifelse(
      DPG_EE_vs_normal_padj< 0.1 & lrt_padj < .1,
      "Significant at FDR 10%",
      "Not significant at FDR 10%"
    )
  ) %>% 
  na.omit() %>% 
  ggplot(aes(DPG_EE_vs_normal, -log10(DPG_EE_vs_normal_pvalue))) +
  ggrepel::geom_text_repel(
    aes(label = str_replace_all(bac, "_", " ")),
    data = . %>% filter(!str_detect(sig, "Not")),
    max.time = 5, nudge_x = -1
  ) +
  geom_point(aes(color = sig)) +
  labs(
    color = NULL,
    x = "log2FC",
    y = "-log10(p-value)",
    title = "DPG+EE vs normal"
  ) +
  facet_wrap(~tool, ncol = 5) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 5))

.p <- guide_area() + (p1/p2/p3) +
  plot_layout(guides = "collect",
              nrow = 2,
              heights = c(1, 20)) & theme(legend.position = 'top') 
ggsave(
  f("{outdir}/differential-abundance/volcano-all.png"),
  .p, width = 14, height = 10
)

# accont for tools
all_res %>% 
  group_by(tax_rank) %>% 
  mutate(lrt_padj2 = p.adjust(lrt_pvalue, method = 'BH')) %>% 
  select(bac, tax_rank, tool, contains("lrt")) %>% 
  filter(lrt_padj2<.1)

# account for tax ranks
all_res %>%
  group_by(tool) %>% 
  mutate(lrt_padj2 = p.adjust(lrt_pvalue, method = 'BH')) %>% 
  select(bac, tax_rank, tool, contains("lrt")) %>% 
  filter(lrt_padj2<.1)

# account for both
all_res %>%
  mutate(lrt_padj2 = p.adjust(lrt_pvalue, method = 'BH')) %>% 
  select(bac, tax_rank, tool, contains("lrt")) %>% 
  filter(lrt_padj2<.1)

## potential-ish hit ----

.bac <- "Tenericutes"
.tax_rank <- "phylum"
all_res %>% 
  filter(str_detect(bac, .bac), tax_rank == .tax_rank) %>% 
  view

df <- get_taxon_couns(phylo.prop, .bac, .tax_rank)


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
      "]\n(",
      x, "/", n, ")"
    ),
    value = 0.015
  ) %>% 
  dplyr::select(endoscopia, value, .lab)

df %>% 
  ggplot(aes(endoscopia, value)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(
    color = "gray40", pch = 21, size = 3,
    aes(fill = coleta),
    width = .3
  ) +
  stat_summary(fun.data = mean_ci, col='red') +
  geom_text(
    data = prev_data,
    aes(label = .lab)
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  scale_fill_manual(
    values = .cols$coleta
  ) +
  labs(
    y = "Relative abundance",
    x = NULL, fill = NULL,
    subtitle = str_glue("{.bac} ({.tax_rank})")
  )

ggsave(
  f("{outdir}/differential-abundance/{.bac}-{.tax_rank}.png"),
  width = 10, height = 4.75
)


# Overall Variation ----

phylo_prop_list <- map(rank_list, ~{
  phylo.prop %>% 
    tax_glom2(.x)
})

for (.rank in rank_list) {
  cat(.rank, sep = "\n")
  
  phylo.prop.rank <- phylo_prop_list[[.rank]]
  
  bacs_ix <- colMeans(otu_table(phylo.prop.rank) > 0) >= min_prevalence
  others <- phylo.prop.rank@otu_table@.Data[, !bacs_ix]
  others <- ifelse(
    is.vector(others), others, rowSums(others)
  )
  counts_matrix <- cbind(
    phylo.prop.rank@otu_table@.Data[,bacs_ix],
    "Others (<5% prevalence)" = others
  )
  
  taxa_order <- order(colMeans(counts_matrix), decreasing = T)
  counts_matrix <- counts_matrix[, taxa_order]
  
  file_path <- f('{outdir}/relative-abundances/{.rank}.png')
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
      mutate(`Any alteration` = ifelse(endoscopia == "normal", "no", "yes"),
             Center = str_replace(coleta, "-", " ")) %>%
      select(sample_id, `Any alteration`, Endoscopy := endoscopia, Center) %>% 
      column_to_rownames("sample_id"),
    show_colnames = F,
    cluster_rows = FALSE,
    scale = 'none',
    color = viridis::viridis(100), 
    legend_breaks = seq(0, 100, 25),
    legend_labels = paste0(seq(0, 100, 25), "%")
  )
  dev.off()
  
  top_taxa <- as.vector(na.omit(colnames(counts_matrix)[1:20]))
  .p <- cbind(
    phylo.prop.rank@sam_data, counts_matrix[, top_taxa]
  ) %>% 
    as_tibble(rownames = "sample_id") %>% 
    pivot_longer(cols = all_of(unname(top_taxa))) %>% 
    ggplot(aes(fct_relevel(name, top_taxa), value, color = coleta)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_color_manual(values = .cols$coleta) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
          legend.position = c(.8, .8)) +
    labs(x = NULL, color = NULL, y = "Abundance (%)")
  ggsave(
    f("{outdir}/relative-abundances/top-taxa-{.rank}.png"),
    .p, width = 12, height = 6
  )
}