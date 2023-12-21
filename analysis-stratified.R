library(tidyverse)
library(phyloseq)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)
library(DESeq2)
library(corncob)
library(limma)
library(edgeR)  # BiocManager::install("edgeR")
library(here)


# Basic setup -----------------------------------------------------
source(here("R/utils.R"))
theme_set(theme_bw(base_size = 13))
outdir <- here("output/stratified-analysis")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# Load data ---------------------------------------------------------------------
.otu <- read_tsv(here("data/220326-195253_otu_table.tsv"),
                 show_col_types = FALSE) %>% 
  mutate(sample = str_extract(sample, "\\d+")) %>% 
  column_to_rownames("sample") %>% 
  as.matrix()
.tax <- read_tsv(here("data/220326-195253_tax_table.tsv"),
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
      endoscopia == "doença péptica gastroduodenal" ~ "GPD",
      endoscopia == "esofagite erosiva" ~ "EE",
      endoscopia == "esofagite erosiva, doença péptica gastroduodenal" ~ "GPD+EE",
      TRUE ~ NA_character_
    ) %>% 
      factor(
        levels = c(
          "normal",
          "GPD",
          "EE",
          "GPD+EE"
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

coletas <- .meta$coleta %>% 
  as.character() %>% 
  unique() %>% 
  set_names(.)

phylo_list <- map(
  coletas, \(.coleta) {
  .phyl <- phyloseq(
    otu_table(.otu, taxa_are_rows = FALSE),
    tax_table(.tax),
    .meta %>% 
      filter(coleta == .coleta) %>% 
      sample_data()
  )
  .phyl <- prune_taxa(
    taxa_sums(.phyl) > 0,
    .phyl
  )
  .phyl <- prune_samples(
    sample_sums(.phyl) > 0,
    .phyl
  )
  return(.phyl)
})


## colors ----

.cols <- list(
  endoscopia = RColorBrewer::brewer.pal(4, "Dark2") %>% 
    set_names(nm = levels(.meta$endoscopia)),
  y = RColorBrewer::brewer.pal(3, "Set1")[1:2] %>% 
    set_names(nm = levels(.meta$y))
)

## alpha-diversity ----

alpha_metrics <- c("Richness", "Shannon", "InvSimpson") %>% 
  set_names(.)
alpha_results <- map(coletas, \(.coleta) {
  .phyl <- phylo_list[[.coleta]]
  .df <- left_join(
    # alpha div measures
    estimate_richness(.phyl, measures = c("Observed", alpha_metrics[-1])) %>% 
      rownames_to_column("sample_id") %>% 
      mutate(sample_id = str_extract(sample_id, "\\d+")) %>% 
      dplyr::rename(
        Richness := Observed
      ),
    # metadata
    data.frame(sample_data(.phyl)) %>% 
      rownames_to_column("sample_id"),
    
    by = "sample_id"
  ) %>% 
    pivot_longer(cols = all_of(unname(alpha_metrics)))
  
  walk(alpha_metrics, \(.metric) {
    .df_alpha <- .df %>% 
      dplyr::filter(name == .metric)
    
    pvals <- tribble(
      ~name, ~pval,
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
        x = "Endoscopy", y = .metric,
        subtitle = paste0(
          "Wilcoxon Rank Sum Test p=", 
          round(pvals$y$padj, 3)
        )
      ) +
      guides(fill = "none")
    
    .metric <- str_to_lower(.metric)
    
    ggsave(
      str_glue("{outdir}/alpha-diversity/alpha-diversity-{.coleta}-{.metric}-endoscopy.png"),
      p1, 
      width = 6, height = 4.5
    )
    ggsave(
      str_glue("{outdir}/alpha-diversity/alpha-diversity-{.coleta}-{.metric}-endoscopy-binary.png"),
      p2, 
      width = 6, height = 4.5
    )
  })
})

## beta-diversity ----

walk(coletas, \(.coleta) {
  phylo.prop <- transform_sample_counts(phylo_list[[.coleta]], function(i) i/sum(i))
  
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
    .distance ~ endoscopia, 
    data = data.frame(sample_data(phylo.prop)),
    by = "margin",
    permutations = 3000
  ) %>% 
    get_pmvn(.names = c("endoscopia"))
  
  pmnv_binary <- vegan::adonis2(
    .distance ~ y, 
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
    str_glue("{outdir}/beta-diversity/beta-diversity-endoscopy-{.coleta}.png"),
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
    str_glue("{outdir}/beta-diversity/beta-diversity-endoscopy-binary-{.coleta}.png"),
    .p, 
    width = 12, height = 4.5,bg = 'white'
  )
})


# Differential Abundance ----

rank_list <- c(
  "lowest",
  "species",
  "genus",
  "family",
  "phylum"
) %>% 
  set_names(.)

min_prevalence <- 0.05
da_res <- map(coletas, \(.coleta) {
  output <- vector("list", length = 4)
  cat(cli::bg_red(str_glue("Running DA analysis for {.coleta}")))
  map_df(
    rank_list,
    function(.rank) {
      print(str_glue("\n\nRank: {.rank}\n"))
      phylo.rank <- tax_glom2(phylo_list[[.coleta]], .rank = .rank)
      phylo.rank <- prune_taxa(
        # taxa must be present in at least `min_prevalence` of samples
        taxa_prevalence(phylo.rank) >= min_prevalence,
        phylo.rank
      )
      phylo.rank@sam_data$endoscopia <- fct_relevel(
        str_replace(phylo.rank@sam_data$endoscopia, "\\+", "_"),
        "normal", "GPD", "EE", "GPD_EE"
      )
      
      # DESeq2 
      res_deseq2 <- phylo.rank %>%
        phyloseq_to_deseq2(~ endoscopia) %>%
        DESeq(sfType = 'poscounts', test = "LRT", reduced = ~ 1, quiet = T) %>%
        results(tidy = T) %>% 
        select(bac := row, 
               padj) %>% 
        mutate(tool = "DESeq2")
      
      # Corncob
      .contrasts <- list(
        "endoscopiaGPD",
        "endoscopiaEE"
      ) %>% set_names(.)
      
      if (.coleta %in% c("coleta-2", "Center 2")) {
        .contrasts[["endoscopiaGPD_EE"]] <- "endoscopiaGPD_EE"
      }
      
      d <- contrastsTest(formula = ~ coleta + endoscopia,
                         phi.formula = ~ coleta,
                         contrasts_DA = .contrasts,
                         data = phylo.rank,
                         test = "Wald", 
                         fdr_cutoff = 1)
      res_cc <- map_df(
        seq_along(.contrasts), ~{
        tibble(
          bac = taxa_names(d$data),
          term = .contrasts[[.x]],
          padj = d$p_fdr[, .x],
          tool = "corncob"
        )
      })
      
      ## limma + voom
      
      counts_table <- t(phylo.rank@otu_table@.Data)
      rownames(counts_table) <- str_replace(rownames(counts_table), " ", "_")
      
      dge_list <- edgeR::DGEList(counts_table) %>% 
        edgeR::calcNormFactors(method = "TMM")
      
      design_matrix <- model.matrix(~ endoscopia, 
                                    data.frame(phylo.rank@sam_data))
      voomvoom <- voom(dge_list, design_matrix, plot = FALSE)
      fit <- lmFit(voomvoom, design_matrix) %>% 
        eBayes()
      
      res_limma <- .contrasts %>% 
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
              ) %>% 
              mutate(
                tool = "limma+voom"
              )
          }
        )
      
      taxa_list <- taxa_names(phylo.rank) %>% 
        set_names(.)
      res_logistf <- map_df(taxa_list, \(.taxname) {
        .df <- get_taxon_couns(.phyloseq = phylo.rank, taxon = .taxname) %>% 
          mutate(value = as.numeric(value > 0))
        fit <- logistf::logistf(
          value ~ endoscopia, data = .df
        )
        summary_logistf(fit, .print = FALSE) %>% 
          as_tibble(rownames = "term") %>% 
          filter(str_detect(term, "endoscopia")) %>% 
          dplyr::select(
            term, 
            estimate := coef,
            .lower := `lower 0.95`,
            .upper := `upper 0.95`,
            pval := p
          )
      }, .id = "bac") %>% 
        mutate(
          padj = p.adjust(pval, method = "BH"),
          tool = "logistf"
        )
      
      all_res <- bind_rows(
        res_deseq2 %>% dplyr::filter(padj < 0.1),
        res_cc %>% dplyr::filter(padj < 0.1),
        res_limma %>% dplyr::filter(padj < 0.1),
        res_logistf %>% dplyr::filter(padj < 0.1)
      )
      return(all_res)
    }, 
    .id = "tax_rank"
  )
  
})

da_res %>% 
  map_df(
    ~ .x %>% 
      group_by(tax_rank, bac) %>% 
      summarise(n_tools = n_distinct(tool),
                tools = paste0(unique(tool), collapse = "; "),
                .groups = "drop"),
    .id = "coleta"
  ) %>% 
  arrange(-n_tools)

# no taxa identified by more than one tool
# in both samples. Only the Tenericutes phylum
# identified by more than one tool in first sample
# (limma+voom and logistf),
# though it was only identified by a single tool in
# the second sample. Moreover, the association
# in the first sample was with the EE group, whereas
# the association in the second sample was with the
# GPD+EE group. All reads from this phylum came from the
# Mycoplasma genus, though associations at lower levels likely
# diseapeared after FDR adjustment. Given the inconsistent
# results across 4 DA tools and between the two samples,
# this association is unlikely to replicate.

# likely just noise
phylo_list %>% 
  map_df(
    ~ transform_sample_counts(.x, \(j) j / sum(j)) %>% 
      get_taxon_couns("Tenericutes", .rank = "phylum")
  ) %>% 
  ggplot(aes(endoscopia, value)) +
  geom_boxplot() +
  geom_quasirandom() +
  facet_wrap(~coleta, scales = 'free') +
  labs(
    subtitle = paste0(
      "Tenericutes association\nEE vs normal (coleta 1, FDR 10%)\nGPD+EE vs normal (coleta 2, FDR 10%)"
    )
  )

