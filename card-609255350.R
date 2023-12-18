if (!require(Maaslin2)) {
  BiocManager::install("Maaslin2")
  require(Maaslin2)
}
library(tidyverse)

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

## colors ----

.cols <- list(
  endoscopia = RColorBrewer::brewer.pal(4, "Dark2") %>% 
    set_names(nm = levels(phylo@sam_data$endoscopia)),
  y = RColorBrewer::brewer.pal(3, "Set1")[1:2] %>% 
    set_names(nm = levels(phylo@sam_data$y)),
  coleta = RColorBrewer::brewer.pal(3, "Set2")[1:2] %>% 
    set_names(nm = levels(phylo@sam_data$coleta))
)
