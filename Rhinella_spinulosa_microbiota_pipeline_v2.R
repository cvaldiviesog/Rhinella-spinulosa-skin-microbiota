# ==============================================================================
# Rhinella spinulosa skin microbiota
# Analysis pipeline — version 2
# ==============================================================================
#
# Purpose
# -------
# This script contains the analytical and statistical workflow for the
# Rhinella spinulosa skin microbiota study. Figure generation is kept in a
# separate script.
#
# Reproducibility notes
# ---------------------
# 1. The Farellones taxonomy is read from data/Rhinella2_taxonomy.tsv.
# 2. Network metrics reported in Table 1 were estimated in Cytoscape and are
#    therefore not recalculated here as substitutes.
# 3. The n = 4 sensitivity analysis uses five independent random selections of
#    four T and four J replicates; the A dataset already contains four replicates.
# 4. Topological states are assigned independently within each group-specific
#    network. No inferential comparisons between networks are performed.
# 5. This script has been structurally refactored from the original analysis
#    script. It has not been executed in the current environment because R is
#    unavailable here.
#
# Expected data files in data/
# ----------------------------
# Rhinella2_feature-table.biom
# sample_data.txt
# Rhinella2_taxonomy.tsv
# representative-sequences_FAR.fasta
# representative-sequences_Longo.fasta
# data_feature-table2.biom
# sample_data_longo.txt
# taxonomy_Longo.tsv
#
# ==============================================================================

# ---- 0. Setup ----------------------------------------------------------------

options(stringsAsFactors = FALSE)

DATA_DIR <- "data"
RESULTS_DIR <- "results"
FIGURES_DIR <- "figures"

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
if (!dir.exists(FIGURES_DIR)) dir.create(FIGURES_DIR, recursive = TRUE)

BIOM_FILE <- file.path(DATA_DIR, "Rhinella2_feature-table.biom")
METADATA_FILE <- file.path(DATA_DIR, "sample_data.txt")
TAXONOMY_FILE <- file.path(DATA_DIR, "Rhinella2_taxonomy.tsv")
FAR_SEQS_FILE <- file.path(DATA_DIR, "representative-sequences_FAR.fasta")
LONGO_SEQS_FILE <- file.path(DATA_DIR, "representative-sequences_Longo.fasta")
LONGO_BIOM_FILE <- file.path(DATA_DIR, "data_feature-table2.biom")
LONGO_METADATA_FILE <- file.path(DATA_DIR, "sample_data_longo.txt")
LONGO_TAXONOMY_FILE <- file.path(DATA_DIR, "taxonomy_Longo.tsv")

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
  library(biomformat)
  library(tidyverse)
  library(ComplexUpset)
  library(SpiecEasi)
  library(igraph)
  library(RCy3)
  library(stringr)
  library(rstatix)
  library(ggpubr)
  library(patchwork)
  library(ggraph)
  library(graphlayouts)
  library(tidytext)
  library(cowplot)
  library(ANCOMBC)
  library(ggbeeswarm)
  library(rcompanion)
  library(boot)
  library(ggalluvial)
})

# ---- 1. Import and preprocess Farellones data -------------------------------

# Import input files
biom <- read_biom(BIOM_FILE)

otu_mat <- as.matrix(biom_data(biom))

OTU <- otu_table(
  otu_mat,
  taxa_are_rows = TRUE
)
meta <- read.delim(
  METADATA_FILE,
  sep="\t",
  stringsAsFactors = FALSE
)

rownames(meta) <- meta$sample.id
meta <- meta[colnames(otu_mat), ]

SAM <- sample_data(meta)

# Taxonomy table used for the Farellones dataset.
# The expected input is a tab-delimited taxonomy table with one row per ASV.
# If the file contains a single semicolon-delimited Taxon column, convert it
# to rank columns before constructing the phyloseq tax_table object.
# ---- Taxonomy ----------------------------------------------------------------

tax_raw <- read.delim(
  TAXONOMY_FILE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

tax_split <- strsplit(
  as.character(tax_raw$Taxon),
  ";"
)

tax_ranks <- c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

tax_mat <- matrix(
  NA_character_,
  nrow = nrow(tax_raw),
  ncol = length(tax_ranks),
  dimnames = list(
    tax_raw$`Feature ID`,
    tax_ranks
  )
)

for (i in seq_along(tax_split)) {
  
  x <- tax_split[[i]]
  
  n <- min(length(x), length(tax_ranks))
  
  if (n > 0) {
    tax_mat[i, seq_len(n)] <- x[seq_len(n)]
  }
}

# Remove taxonomic prefixes
tax_mat <- apply(
  tax_mat,
  2,
  function(x) sub("^[a-z]__", "", x)
)

rownames(tax_mat) <- tax_raw$`Feature ID`

# Ensure exact correspondence and order with OTU table
stopifnot(
  setequal(taxa_names(OTU), rownames(tax_mat))
)

tax_mat <- tax_mat[
  taxa_names(OTU),
  ,
  drop = FALSE
]

TAX <- tax_table(tax_mat)

ps <- phyloseq(
  
  OTU,
  TAX,
  SAM
  
)


# ---- 2. Taxonomic filtering ---------------------------------------------------

ps <- subset_taxa(
  ps,
  Kingdom == "Bacteria"
)

ps <- subset_taxa(
  ps,
  Family != "Mitochondria"
)

ps <- subset_taxa(
  ps,
  Order != "Chloroplast"
)

ps <- prune_taxa(
  taxa_sums(ps) > 0,
  ps
)

ps_far <- subset_samples(
  ps,
  environment == "skin"
)

ps_water <- subset_samples(
  ps,
  environment == "water"
)

ps_far <- prune_taxa(
  taxa_sums(ps_far) > 0,
  ps_far
)

ps_water <- prune_taxa(
  taxa_sums(ps_water) > 0,
  ps_water
)


# ---- 3. Abundance matrices ----------------------------------------------------

otu_skin  <- as(otu_table(ps_far), "matrix")
otu_water <- as(otu_table(ps_water), "matrix")

# Relative-abundance matrices used throughout the host-filtering analyses.
rel_skin <- sweep(otu_skin, 2, colSums(otu_skin), "/")
rel_water <- sweep(otu_water, 2, colSums(otu_water), "/")

## Ensure that ASVs are arranged in rows
if(!taxa_are_rows(ps_far)){
  otu_skin <- t(otu_skin)
}

if(!taxa_are_rows(ps_water)){
  otu_water <- t(otu_water)
}

# Calculate summary statistics by environment
summary_skin <- data.frame(
  
  FeatureID = rownames(otu_skin),
  
  Prev_skin = rowSums(otu_skin > 0),
  
  Prev_skin_prop = rowSums(otu_skin > 0) / ncol(otu_skin),
  
  MeanCount_skin = rowMeans(otu_skin),
  
  TotalCount_skin = rowSums(otu_skin)
  
)

summary_water <- data.frame(
  
  FeatureID = rownames(otu_water),
  
  Prev_water = rowSums(otu_water > 0),
  
  Prev_water_prop = rowSums(otu_water > 0) / ncol(otu_water),
  
  MeanCount_water = rowMeans(otu_water),
  
  TotalCount_water = rowSums(otu_water)
  
)

# Combine the two summary tables
asv_master <- full_join(
  
  summary_skin,
  summary_water,
  
  by="FeatureID"
  
)

asv_master[is.na(asv_master)] <- 0

tax_df <- as.data.frame(
  tax_table(ps)
)

tax_df$FeatureID <- rownames(tax_df)

asv_master <- left_join(
  
  asv_master,
  tax_df,
  
  by="FeatureID"
  
)


# ---- 4. Host-filtering categories ---------------------------------------------

## Maximum observed relative abundance per ASV
max_skin <- apply(rel_skin, 1, max)
max_water <- apply(rel_water, 1, max)

max_skin_df <- data.frame(
  FeatureID = names(max_skin),
  MaxRel_skin = max_skin
)

max_water_df <- data.frame(
  FeatureID = names(max_water),
  MaxRel_water = max_water
)

sens_master <- asv_master %>%
  left_join(max_skin_df, by = "FeatureID") %>%
  left_join(max_water_df, by = "FeatureID") %>%
  mutate(
    MaxRel_skin = replace_na(MaxRel_skin, 0),
    MaxRel_water = replace_na(MaxRel_water, 0)
  )

thresholds <- c(
  "0%"    = 0,
  "0.01%" = 0.0001,
  "0.1%"  = 0.001
)

# Core ASVs
core_skin <- taxa_names(
  prune_taxa(
    rowSums(otu_skin > 0) >= 4,
    ps_far
  )
)

# Add core status to ASV master table
asv_master$core_skin <-
  asv_master$FeatureID %in% core_skin

sensitivity_water <- bind_rows(
  
  lapply(names(thresholds), function(th){
    
    x <- thresholds[th]
    
    sens_master %>%
      left_join(
        asv_master %>%
          select(FeatureID, core_skin),
        by = "FeatureID"
      ) %>%
      mutate(
        
        ## Skin occurrence and core membership remain fixed
        skin_detected = Prev_skin > 0,
        
        ## Only environmental detection threshold changes
        water_detected = MaxRel_water > x,
        
        Category_sens = case_when(
          
          !skin_detected & water_detected ~ "Water_only",
          
          skin_detected & water_detected & !core_skin ~
            "Shared_non_core",
          
          skin_detected & !water_detected & !core_skin ~
            "Skin_only_non_core",
          
          skin_detected & water_detected & core_skin ~
            "Core_detected_in_water",
          
          skin_detected & !water_detected & core_skin ~
            "Core_not_detected_in_water",
          
          TRUE ~ "Not_detected"
        ),
        
        Threshold = th
      )
  })
)





## Focus only on ASVs occurring on tadpole skin

sensitivity_skin_ASVs <- sensitivity_water %>%
  filter(Prev_skin > 0) %>%
  count(Threshold, Category_sens) %>%
  group_by(Threshold) %>%
  mutate(
    Percent = 100 * n / sum(n)
  ) %>%
  ungroup()

sensitivity_skin_ASVs

## Number of skin ASVs considered environmentally shared
## at each threshold

sensitivity_shared <- sensitivity_water %>%
  filter(Prev_skin > 0) %>%
  group_by(Threshold) %>%
  summarise(
    Total_skin_ASVs = n(),
    Shared = sum(
      Category_sens %in%
        c("Shared_non_core",
          "Core_detected_in_water")
    ),
    Skin_restricted = sum(
      Category_sens %in%
        c("Skin_only_non_core",
          "Core_not_detected_in_water")
    ),
    Percent_shared = 100 * Shared / Total_skin_ASVs,
    Percent_skin_restricted =
      100 * Skin_restricted / Total_skin_ASVs,
    .groups = "drop"
  )

sensitivity_shared

## ==========================================
## Sensitivity analysis: detection threshold
## ==========================================

thresholds <- c(
  "0%"    = 0,
  "0.01%" = 0.0001,
  "0.1%"  = 0.001
)

classify_threshold <- function(threshold){
  
  skin_detected <- rowSums(rel_skin > threshold) > 0
  water_detected <- rowSums(rel_water > threshold) > 0
  
  tibble(
    FeatureID = rownames(rel_skin),
    Detection_threshold = case_when(
      skin_detected & !water_detected ~ "Skin_only",
      !skin_detected & water_detected ~ "Water_only",
      skin_detected & water_detected  ~ "Shared",
      TRUE                            ~ "Not_detected"
    )
  )
}

sensitivity_detection <- bind_rows(
  lapply(
    names(thresholds),
    function(x){
      classify_threshold(thresholds[x]) %>%
        mutate(Threshold = x)
    }
  )
)

## Number of ASVs in each category at each threshold
sensitivity_summary <- sensitivity_detection %>%
  count(Threshold, Detection_threshold)

sensitivity_summary

## ==========================================
## Sensitivity of ecological classification
## ==========================================

sensitivity_categories <- sensitivity_detection %>%
  left_join(
    asv_master %>%
      select(FeatureID, core_skin),
    by = "FeatureID"
  ) %>%
  mutate(
    Category = case_when(
      Detection_threshold == "Water_only" ~ "Water_only",
      Detection_threshold == "Shared" & !core_skin ~ "Shared_non_core",
      Detection_threshold == "Skin_only" & !core_skin ~ "Skin_only_non_core",
      Detection_threshold == "Shared" & core_skin ~ "Core_detected_in_water",
      Detection_threshold == "Skin_only" & core_skin ~ "Core_not_detected_in_water",
      TRUE ~ "Not_detected"
    )
  )

sensitivity_category_summary <- sensitivity_categories %>%
  count(Threshold, Category) %>%
  group_by(Threshold) %>%
  mutate(
    Percent = 100 * n / sum(n)
  ) %>%
  ungroup()

sensitivity_category_summary


## Compare classifications with the original 0% threshold

baseline <- sensitivity_categories %>%
  filter(Threshold == "0%") %>%
  select(FeatureID, Baseline = Category)

sensitivity_changes <- sensitivity_categories %>%
  left_join(baseline, by = "FeatureID") %>%
  mutate(
    Changed = Category != Baseline
  ) %>%
  group_by(Threshold) %>%
  summarise(
    Total = n(),
    Changed = sum(Changed),
    Unchanged = sum(!Changed),
    Percent_changed = 100 * Changed / Total,
    .groups = "drop"
  )

sensitivity_changes


# ---- 5. Core microbiota -------------------------------------------------------

rel_skin_df <- data.frame(
  FeatureID = rownames(rel_skin),
  MeanRel_skin = rowMeans(rel_skin)
)

rel_water_df <- data.frame(
  FeatureID = rownames(rel_water),
  MeanRel_water = rowMeans(rel_water)
)

asv_master <- left_join(
  asv_master,
  rel_skin_df,
  by = "FeatureID"
)

asv_master <- left_join(
  asv_master,
  rel_water_df,
  by = "FeatureID"
)

asv_master$MeanRel_skin[is.na(asv_master$MeanRel_skin)] <- 0
asv_master$MeanRel_water[is.na(asv_master$MeanRel_water)] <- 0

asv_master <- asv_master %>%
  mutate(
    PrevClass = case_when(
      Prev_skin == 5 ~ "100%",
      Prev_skin == 4 ~ "80%",
      Prev_skin == 3 ~ "60%",
      Prev_skin == 2 ~ "40%",
      Prev_skin == 1 ~ "20%",
      TRUE           ~ "Absent"
    )
  )

# Original detection classification
asv_master <- asv_master %>%
  mutate(
    Detection = case_when(
      Prev_skin > 0 & Prev_water == 0 ~ "Skin_only",
      Prev_skin == 0 & Prev_water > 0 ~ "Water_only",
      Prev_skin > 0 & Prev_water > 0 ~ "Shared"
    )
  )

# Ecological classification incorporating core status
asv_master <- asv_master %>%
  mutate(
    Category = case_when(
      Detection == "Water_only" ~ "Water_only",
      Detection == "Shared" & !core_skin ~ "Shared_non_core",
      Detection == "Skin_only" & !core_skin ~ "Skin_only_non_core",
      Detection == "Shared" & core_skin ~ "Core_detected_in_water",
      Detection == "Skin_only" & core_skin ~ "Core_not_detected_in_water"
    )
  )

##


## ===========================
## Presence/absence sets
## ===========================

asv_sets <- asv_master %>%
  dplyr::transmute(
    
    `Environmental pool` = Prev_water > 0,
    
    `Skin microbiota` = Prev_skin > 0,
    
    `Core microbiota` = core_skin
    
  )

skin_col  <- "#2C7FB8"
core_col  <- "#FDAE61"
water_col <- "#41AE76"


# ---- 6. Differential abundance (ANCOM-BC2) ----------------------------------

ancom <- ancombc2(
  
  data = ps,
  
  assay_name = "counts",
  
  tax_level = NULL,
  
  fix_formula = "environment",
  
  rand_formula = NULL,
  
  p_adj_method = "BH",
  
  prv_cut = 0,
  
  lib_cut = 0,
  
  group = "environment",
  
  struc_zero = TRUE,
  
  neg_lb = TRUE,
  
  alpha = 0.05,
  
  global = FALSE,
  
  pairwise = FALSE,
  
  dunnet = FALSE,
  
  trend = FALSE,
  
  iter_control = list(
    tol = 1e-2,
    max_iter = 20,
    verbose = TRUE
  ),
  
  em_control = list(
    tol = 1e-5,
    max_iter = 100
  ),
  
  lme_control = lme4::lmerControl(),
  
  mdfdr_control = list(
    fwer_ctrl_method = "holm",
    B = 100
  )
  
)


# Add ANCOM-BC2 results to the master table

ancom_res <- ancom$res %>%
  rename(FeatureID = taxon)

asv_master <- left_join(
  asv_master,
  ancom_res,
  by = "FeatureID"
)


# Summarize ANCOM-BC2 results by category

ancom_summary <- asv_master %>%
  group_by(Category) %>%
  summarise(
    Total = n(),
    Tested = sum(!is.na(diff_environmentwater)),
    Significant = sum(diff_environmentwater, na.rm = TRUE),
    Percent_tested = round(100 * Tested / Total, 1),
    Percent_significant = round(100 * Significant / Tested, 1)
  )


### Plots de prevalencia

# Are core ASVs more abundant?

kw_data <- asv_master %>%
  filter(
    MeanRel_skin > 0
  ) %>%
  mutate(
    Category = recode(
      Category,
      Water_only = "Water only",
      Shared_non_core = "Shared (non-core)",
      Skin_only_non_core = "Skin only (non-core)",
      Core_detected_in_water = "Core (shared)",
      Core_not_detected_in_water = "Core (skin-restricted)"
    ),
    Category = factor(
      Category,
      levels = c(
        "Water only",
        "Shared (non-core)",
        "Skin only (non-core)",
        "Core (shared)",
        "Core (skin-restricted)"
      )
    )
  )

kruskal.test(
  MeanRel_skin ~ Category,
  data = kw_data
)

dunn_res <- kw_data %>%
  dunn_test(
    MeanRel_skin ~ Category,
    p.adjust.method = "BH"
  )


# Do core ASVs show higher prevalence?

kw_prev <- asv_master %>%
  filter(
    Category != "Water_only"
  ) %>%
  mutate(
    Category = factor(
      Category,
      levels = c(
        "Shared_non_core",
        "Skin_only_non_core",
        "Core_detected_in_water",
        "Core_not_detected_in_water"
      ),
      labels = c(
        "Shared (non-core)",
        "Skin only (non-core)",
        "Core (shared)",
        "Core (skin-restricted)"
      )
    )
  )

kruskal.test(
  Prev_skin ~ Category,
  data = kw_prev
)

dunn_prev <- kw_prev %>%
  dunn_test(
    Prev_skin ~ Category,
    p.adjust.method = "BH"
  )

prev_summary <- kw_prev %>%
  group_by(Category) %>%
  summarise(
    Median = median(Prev_skin) * 20,
    Q1 = quantile(Prev_skin, 0.25) * 20,
    Q3 = quantile(Prev_skin, 0.75) * 20,
    .groups = "drop"
  )

prev_summary$Letters <- c(
  "c",
  "b",
  "a",
  "a"
)


# Is category associated with ANCOM-BC2 significance?
# Category versus significance (chi-square or Fisher's exact test)

ancom_tested <- asv_master %>%
  filter(
    StructuralZero == FALSE,
    !is.na(diff_robust_environmentwater)
  ) %>%
  mutate(
    ANCOM = ifelse(
      diff_robust_environmentwater,
      "Significant",
      "Not significant"
    )
  )

tab <- table(
  ancom_tested$Category,
  ancom_tested$ANCOM
)

chisq.test(tab)
cramerV(tab)


# Is prevalence associated with abundance?

cor.test(
  asv_master$Prev_skin,
  asv_master$MeanRel_skin,
  method = "spearman",
  exact = FALSE
)

plot_cor <- asv_master %>%
  filter(
    Category != "Water_only"
  ) %>%
  mutate(
    Category = factor(
      Category,
      levels = c(
        "Shared_non_core",
        "Skin_only_non_core",
        "Core_detected_in_water",
        "Core_not_detected_in_water"
      ),
      labels = c(
        "Shared (non-core)",
        "Skin only (non-core)",
        "Core (shared)",
        "Core (skin-restricted)"
      )
    )
  )

spearman <- cor.test(
  plot_cor$Prev_skin,
  plot_cor$MeanRel_skin,
  method = "spearman",
  exact = FALSE
)

label_cor <- paste0(
  "Spearman \u03C1 = ",
  round(spearman$estimate, 2),
  "\nP < 0.001"
)


skin_only_asvs <- asv_master %>%
  filter(
    Detection == "Skin_only",
    !core_skin
  ) %>%
  pull(FeatureID)
# ---- 7. Tadpole association network ------------------------------------------

ps_skin_network <- subset_taxa(
  ps,
  taxa_names(ps) %in% c(core_skin, skin_only_asvs))

ps_core_network <- prune_taxa(
  core_skin,
  ps)

ps_skin_network <- filter_taxa(
  ps_skin_network,
  function(x) sum(x > 0) >= 2,
  TRUE
)

se_skin <- spiec.easi(
  ps_skin_network,
  method = "mb",
  lambda.min.ratio = 1e-2,
  nlambda = 20,
  pulsar.params = list(
    thresh = 0.05
  ))

g_skin <- adj2igraph(
  getRefit(se_skin),
  vertex.attr = list(
    name = taxa_names(ps_skin_network)
  )
)

vertex_df <- asv_master %>%
  filter(
    FeatureID %in% V(g_skin)$name
  )

idx <- match(
  V(g_skin)$name,
  vertex_df$FeatureID
)

V(g_skin)$Kingdom <- vertex_df$Kingdom[idx]
V(g_skin)$Phylum  <- vertex_df$Phylum[idx]
V(g_skin)$Class   <- vertex_df$Class[idx]
V(g_skin)$Order   <- vertex_df$Order[idx]
V(g_skin)$Family  <- vertex_df$Family[idx]
V(g_skin)$Genus   <- vertex_df$Genus[idx]
V(g_skin)$Species <- vertex_df$Species[idx]

createNetworkFromIgraph(
  g_skin,
  title = "Skin_microbiome",
  collection = "Rhinella"
)



### Degree/betweenness correlations versus abundance
FAR_mods <- read.table(
  "Larval_skin_network default node.txt",
  header = TRUE,
  sep = ";",
  quote = "\"",
  stringsAsFactors = FALSE
)

num_cols <- c(
  "AverageShortestPathLength",
  "BetweennessCentrality",
  "ClosenessCentrality",
  "ClusteringCoefficient",
  "MeanCount_skin",
  "MeanCount_water",
  "MeanRel_skin",
  "MeanRel_water",
  "NeighborhoodConnectivity",
  "Prev_skin_prop",
  "Prev_water_prop",
  "Radiality",
  "TopologicalCoefficient"
)

FAR_mods[num_cols] <- lapply(
  FAR_mods[num_cols],
  function(x) as.numeric(gsub(",", ".", x))
)

network_df <- FAR_mods %>%
  
  filter(
    MeanRel_skin > 0
  ) %>%
  
  mutate(
    
    log10MeanRel = log10(MeanRel_skin),
    
    Category = recode(
      Category,
      Shared_non_core = "Shared (non-core)",
      Skin_only_non_core = "Skin only (non-core)",
      Core_detected_in_water = "Core (shared)",
      Core_not_detected_in_water = "Core (skin-restricted)"
    )
    
  )

spearman_degree <- cor.test(
  network_df$log10MeanRel,
  network_df$Degree,
  method = "spearman",
  exact = F
)

pearson_degree <- cor.test(
  network_df$log10MeanRel,
  network_df$Degree,
  method = "pearson"
)

lm_degree <- lm(
  Degree ~ log10MeanRel,
  data = network_df
)

summary(lm_degree)

R2_degree <- summary(lm_degree)$r.squared

spearman_between <- cor.test(
  
  network_df$log10MeanRel,
  
  network_df$BetweennessCentrality,
  
  method = "spearman",
  exact = F
  
)

pearson_between <- cor.test(
  
  network_df$log10MeanRel,
  
  network_df$BetweennessCentrality,
  
  method = "pearson"
  
)

lm_between <- lm(
  
  BetweennessCentrality ~ log10MeanRel,
  
  data = network_df
  
)

R2_between <- summary(lm_between)$r.squared

  

### Degree/betweenness correlations versus ecological category
kw_net <- FAR_mods %>%
  
  filter(
    Category != "Water_only"
  ) %>%
  
  mutate(
    
    Category = factor(
      Category,
      levels = c(
        "Shared_non_core",
        "Skin_only_non_core",
        "Core_detected_in_water",
        "Core_not_detected_in_water"
      ),
      labels = c(
        "Shared (non-core)",
        "Skin only (non-core)",
        "Core (shared)",
        "Core (skin-restricted)"
      )
    )
    
  )

kruskal.test(
  
  Degree ~ Category,
  
  data = kw_net
  
)  

dunn_degree <-
  
  kw_net %>%
  
  dunn_test(
    
    Degree ~ Category,
    
    p.adjust.method = "BH"
    
  )

degree_summary <-
  
  kw_net %>%
  
  group_by(Category) %>%
  
  summarise(
    
    Median = median(Degree),
    
    Q1 = quantile(Degree,0.25),
    
    Q3 = quantile(Degree,0.75),
    
    .groups="drop"
    
  )

degree_summary$Letters <- c(
  "b",
  "a",
  "a",
  "a"
)


kruskal.test(
  
  BetweennessCentrality ~ Category,
  
  data=kw_net
  
)

dunn_between <-
  
  kw_net %>%
  
  dunn_test(
    
    BetweennessCentrality ~ Category,
    
    p.adjust.method="BH"
    
  )

between_summary <-
  
  kw_net %>%
  
  group_by(Category) %>%
  
  summarise(
    
    Median=median(BetweennessCentrality),
    
    Q1=quantile(BetweennessCentrality,.25),
    
    Q3=quantile(BetweennessCentrality,.75),
    
    .groups="drop"
    
  )
between_summary$Letters <- c(
  "c",
  "b",
  "a",
  "a"
)



####Bootstrap de analisis estadisticos
spearman_boot <- function(data, indices){
  
  d <- data[indices, ]
  
  cor(
    d$log10MeanRel,
    d$BetweennessCentrality,
    method = "spearman"
  )
  
}

set.seed(123)

boot_spear <- boot(
  
  data = network_df,
  
  statistic = spearman_boot,
  
  R = 1000
  
)

boot_ci <- boot.ci(
  
  boot_spear,
  
  type = c("perc","bca")
  
)

mean(boot_spear$t)

median(boot_spear$t)

sd(boot_spear$t)

quantile(
  
  boot_spear$t,
  
  c(.025,.975)
  
)

boot_data <- FAR_mods %>%
  
  filter(Category != "Water_only") %>%
  
  mutate(
    
    Core = ifelse(
      
      Category %in% c(
        "Core_detected_in_water",
        "Core_not_detected_in_water"
      ),
      
      "Core",
      
      "Non-core"
      
    )
    
  )

boot_core_between <- function(data, indices){
  
  d <- data[indices, ]
  
  med_core <-
    
    median(
      
      d$BetweennessCentrality[
        d$Core == "Core"
      ]
      
    )
  
  med_non <-
    
    median(
      
      d$BetweennessCentrality[
        d$Core == "Non-core"
      ]
      
    )
  
  med_core - med_non
  
}

set.seed(123)

boot_core <- boot(
  
  data = boot_data,
  
  statistic = boot_core_between,
  
  R = 1000
  
)

boot.ci(
  
  boot_core,
  
  type = c("perc","bca")
  
)


mean(boot_core$t)

median(boot_core$t)

sd(boot_core$t)

quantile(
  
  boot_core$t,
  
  c(.025,.975)
  
)


###Subsampling de test host-filtering
subsample_spearman <- function(data, prop_remove = 0.1){
  
  n_remove <- floor(nrow(data) * prop_remove)
  
  keep <- sample(
    seq_len(nrow(data)),
    size = nrow(data) - n_remove,
    replace = FALSE
  )
  
  cor(
    data$log10MeanRel[keep],
    data$BetweennessCentrality[keep],
    method = "spearman"
  )
  
}


set.seed(123)

sub10 <- replicate(
  1000,
  subsample_spearman(network_df,0.10)
)

sub20 <- replicate(
  1000,
  subsample_spearman(network_df,0.20)
)

sub30 <- replicate(
  1000,
  subsample_spearman(network_df,0.30)
)

summarise_sub <- function(x){
  
  tibble(
    
    Mean = mean(x),
    
    SD = sd(x),
    
    Median = median(x),
    
    CI_low = quantile(x,.025),
    
    CI_high = quantile(x,.975)
    
  )
  
}

summarise_sub(sub10)

summarise_sub(sub20)

summarise_sub(sub30)

subsample_core <- function(data, prop_remove = 0.1){
  
  n_remove <- floor(nrow(data) * prop_remove)
  
  keep <- sample(
    seq_len(nrow(data)),
    size = nrow(data)-n_remove,
    replace = FALSE
  )
  
  d <- data[keep,]
  
  median(
    d$BetweennessCentrality[
      d$Core=="Core"
    ]
  ) -
    
    median(
      d$BetweennessCentrality[
        d$Core=="Non-core"
      ]
    )
  
}

set.seed(123)

core10 <- replicate(
  1000,
  subsample_core(boot_data,0.10)
)

core20 <- replicate(
  1000,
  subsample_core(boot_data,0.20)
)

core30 <- replicate(
  1000,
  subsample_core(boot_data,0.30)
)

summarise_sub(core10)

summarise_sub(core20)

summarise_sub(core30)

subsampling_summary <-
  
  bind_rows(
    
    summarise_sub(sub10) %>%
      mutate(
        Analysis="Spearman",
        Removal="10%"
      ),
    
    summarise_sub(sub20) %>%
      mutate(
        Analysis="Spearman",
        Removal="20%"
      ),
    
    summarise_sub(sub30) %>%
      mutate(
        Analysis="Spearman",
        Removal="30%"
      ),
    
    summarise_sub(core10) %>%
      mutate(
        Analysis="Core vs Non-core",
        Removal="10%"
      ),
    
    summarise_sub(core20) %>%
      mutate(
        Analysis="Core vs Non-core",
        Removal="20%"
      ),
    
    summarise_sub(core30) %>%
      mutate(
        Analysis="Core vs Non-core",
        Removal="30%"
      )
    
  )

subsampling_summary




# ---- 8. Network modules and hubs ----------------------------------------------

module_summary <-
  
  FAR_mods %>%
  
  filter(Category != "Water_only") %>%
  
  mutate(
    
    Group = case_when(
      
      Category %in% c(
        "Core_detected_in_water",
        "Core_not_detected_in_water"
      ) ~ "Core",
      
      Category == "Shared_non_core" ~ "Shared",
      
      Category == "Skin_only_non_core" ~ "Skin only"
      
    )
    
  ) %>%
  
  count(X__glayCluster, Group) %>%
  
  group_by(X__glayCluster) %>%
  
  mutate(
    
    Nodes = sum(n),
    
    Percent = round(100 * n / Nodes, 1)
    
  ) %>%
  
  ungroup() %>%
  
  select(
    X__glayCluster,
    Nodes,
    Group,
    Percent
  ) %>%
  
  pivot_wider(
    
    names_from = Group,
    
    values_from = Percent,
    
    values_fill = 0
    
  ) %>%
  
  arrange(desc(Nodes)) %>%
  
  rename(
    
    Module = X__glayCluster,
    
    `% Core` = Core,
    
    `% Shared` = Shared,
    
    `% Skin only` = `Skin only`
    
  )

dominant_genus <-
  
  FAR_mods %>%
  
  filter(Category != "Water_only") %>%
  
  count(X__glayCluster, Genus) %>%
  
  group_by(X__glayCluster) %>%
  
  slice_max(
    n,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup() %>%
  
  rename(
    Module = X__glayCluster,
    Dominant_genus = Genus,
    Genus_nodes = n
  )

module_summary <-
  
  module_summary %>%
  
  left_join(
    dominant_genus,
    by = "Module"
  )

FAR_mods <- FAR_mods %>%
  
  filter(Category != "Water_only") %>%
  
  mutate(
    
    Group = case_when(
      
      Category %in% c(
        "Core_detected_in_water",
        "Core_not_detected_in_water"
      ) ~ "Core",
      
      Category == "Shared_non_core" ~ "Shared",
      
      Category == "Skin_only_non_core" ~ "Skin only"
      
    )
    
  )

tab_modules <-
  
  table(
    FAR_mods$X__glayCluster,
    FAR_mods$Group
  )

chisq_mod <- chisq.test(tab_modules)

chisq_mod$expected

chisq_mod
cramerV(tab_modules)

round(chisq_mod$stdres,2)

## Hub identity

# Degree, BetweennessCentrality y ClosenessCentrality
# were imported from the Cytoscape node table.

FAR_mods$Degree_norm <-
  FAR_mods$Degree/(nrow(FAR_mods)-1)

deg95 <- quantile(
  FAR_mods$Degree_norm,
  0.95
)

bet95 <- quantile(
  FAR_mods$BetweennessCentrality,
  0.95
)

FAR_mods <- FAR_mods %>%
  mutate(
    Hub_degree = Degree_norm >= deg95,
    Hub_between = BetweennessCentrality >= bet95,
    Hub = Hub_degree | Hub_between
  )

degree_set <- FAR_mods$name[FAR_mods$Hub_degree]

between_set <- FAR_mods$name[FAR_mods$Hub_between]

length(intersect(degree_set, between_set)) /
  length(union(degree_set, between_set))

# Hub enrichment
FAR_mods <- FAR_mods %>%
  mutate(
    CoreGroup = ifelse(
      Category %in% c(
        "Core_detected_in_water",
        "Core_not_detected_in_water"
      ),
      "Core",
      "Non-core"
    )
  )

tab_deg <- table(
  FAR_mods$Hub_degree,
  FAR_mods$CoreGroup
)

fisher.test(tab_deg)

tab_btw <- table(
  FAR_mods$Hub_between,
  FAR_mods$CoreGroup
)

fisher.test(tab_btw)

hub_prop <- FAR_mods %>%
  group_by(Category) %>%
  summarise(
    PropHub = mean(Hub_degree) * 100
  )

tab_hub <- table(
  FAR_mods$Category,
  FAR_mods$Hub_degree
)

tab_hub

chisq_hub <- chisq.test(tab_hub)

chisq_hub

chisq_hub$expected
cramerV(tab_hub)
round(chisq_hub$stdres, 2)

hub_prop <-
  
  FAR_mods %>%
  
  group_by(Category) %>%
  
  summarise(
    
    PropHub = mean(Hub_degree) * 100,
    
    Hubs = sum(Hub_degree),
    
    Total = n(),
    
    .groups = "drop"
    
  ) %>%
  
  mutate(
    
    Category = factor(
      
      Category,
      
      levels = c(
        
        "Shared_non_core",
        
        "Skin_only_non_core",
        
        "Core_not_detected_in_water",
        
        "Core_detected_in_water"
        
      ),
      
      labels = c(
        
        "Shared (non-core)",
        
        "Skin only (non-core)",
        
        "Core (skin-restricted)",
        
        "Core (shared)"
        
      )
      
    )
    
  )

hub_prop$Label <-
  
  paste0(
    
    round(hub_prop$PropHub,1),
    
    "%\n(n=",
    
    hub_prop$Hubs,
    
    ")"
    
  )





## Core network

core_asvs <-
  
  asv_master %>%
  
  filter(
    
    core_skin
    
  ) %>%
  
  pull(FeatureID)

ps_core <-
  
  prune_taxa(
    
    core_asvs,
    
    ps_far
    
  )

se_core <-
  
  spiec.easi(
    
    ps_core,
    
    method = "mb",
    
    lambda.min.ratio = 1e-2,
    
    nlambda = 20,
    
    pulsar.params = list(
      
      rep.num = 50
      
    )
    
  )

ig_core <-
  
  adj2igraph(
    
    getRefit(se_core),
    
    vertex.attr = list(
      
      name = taxa_names(ps_core)
      
    )
    
  )

vertex_df <-
  
  asv_master %>%
  
  filter(
    
    FeatureID %in% V(ig_core)$name
    
  )

idx <-
  
  match(
    
    V(ig_core)$name,
    
    vertex_df$FeatureID
    
  )

V(ig_core)$Category <-
  
  vertex_df$Category[idx]

V(ig_core)$Genus <-
  
  vertex_df$Genus[idx]

V(ig_core)$Family <-
  
  vertex_df$Family[idx]

V(ig_core)$Order <-
  
  vertex_df$Order[idx]



createNetworkFromIgraph(
  
  ig_core,
  
  title = "Core network",
  
  collection = "Rhinella skin"
  
)

FAR_mods <- FAR_mods %>%
  mutate(
    Hub_type = case_when(
      Hub_degree & Hub_between ~ "Both",
      Hub_degree               ~ "Degree",
      Hub_between              ~ "Betweenness",
      TRUE                     ~ "None"
    )
  )
hub_nodes <-
  
  FAR_mods %>%
  
  select(
    
    name,
    Hub_type,
    Hub_degree,
    Hub_between,
    Degree,
    BetweennessCentrality,
    Category,
    Genus,
    X__glayCluster
    
  )

degree_hubs <-
  
  FAR_mods %>%
  
  filter(Hub_degree) %>%
  
  pull(name)

between_hubs <-
  
  FAR_mods %>%
  
  filter(Hub_between) %>%
  
  pull(name)

both_hubs <-
  
  FAR_mods %>%
  
  filter(Hub_degree & Hub_between) %>%
  
  pull(name)


V(ig_core)$MeanRel_skin <-
  
  vertex_df$MeanRel_skin[idx]

V(ig_core)$Prev_skin <-
  
  vertex_df$Prev_skin[idx]

V(ig_core)$Degree <-
  
  degree(ig_core)

V(ig_core)$Betweenness <-
  
  betweenness(
    
    ig_core,
    
    normalized = TRUE
    
  )

V(ig_core)$Closeness <-
  
  closeness(
    
    ig_core,
    
    normalized = TRUE
    
  )

V(ig_core)$Degree_norm <-
  
  degree(ig_core)/(vcount(ig_core)-1)

CORE_mods <- data.frame(
  name = V(ig_core)$name,
  MeanRel_skin = V(ig_core)$MeanRel_skin,
  Prev_skin = V(ig_core)$Prev_skin,
  Degree = V(ig_core)$Degree,
  Degree_norm = V(ig_core)$Degree_norm,
  Betweenness = V(ig_core)$Betweenness,
  Closeness = V(ig_core)$Closeness,
  stringsAsFactors = FALSE
)

deg95_core <- quantile(
  CORE_mods$Degree_norm,
  0.95
)

bet95_core <- quantile(
  CORE_mods$Betweenness,
  0.95
)

CORE_mods <- CORE_mods %>%
  mutate(
    Hub_degree = Degree_norm >= deg95_core,
    Hub_between = Betweenness >= bet95_core,
    Hub = Hub_degree | Hub_between
  )

# How many hubs from the full network are also present in the core network?
degree_full <- FAR_mods %>%
  filter(Hub_degree)

between_full <- FAR_mods %>%
  filter(Hub_between)

# Degree hubs also present in the core network
sum(degree_full$name %in% V(ig_core)$name)



# ---- 9. Group-specific network analyses --------------------------------------

####Longo dataset
biom_longo <- read_biom(
  LONGO_BIOM_FILE
)

otu_mat_longo <- as.matrix(
  biom_data(biom_longo)
)

dim(otu_mat_longo)

OTU_longo <- otu_table(
  otu_mat_longo,
  taxa_are_rows = TRUE
)

tax_longo <- read.delim(
  LONGO_TAXONOMY_FILE,
  sep = "\t",
  header = TRUE
)

tax_split_longo <- tax_longo %>%
  separate(
    Taxon,
    into = c(
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    ),
    sep = ";",
    fill = "right"
  )

tax_split_longo <- tax_split_longo %>%
  mutate(
    across(
      Kingdom:Species,
      ~gsub("^[a-z]__", "", .)
    )
  )

rownames(tax_split_longo) <-
  tax_split_longo$Feature.ID

tax_mat_longo <- as.matrix(
  tax_split_longo[,c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  )]
)

TAX_longo <- tax_table(
  tax_mat_longo
)

meta_longo <- data.frame(
  SampleID = colnames(
    otu_mat_longo
  )
)

meta_longo$Stage <- ifelse(
  grepl("^YL", meta_longo$SampleID),
  "Juvenile",
  "Adult"
)

rownames(meta_longo) <-
  meta_longo$SampleID

META_longo <- sample_data(
  meta_longo
)

### Function to process count matrices

build_phyloseq_dataset <- function(
    biom_file,
    metadata_file,
    taxonomy_file,
    group_column,
    remove_unclassified = TRUE){
  
  library(phyloseq)
  library(biomformat)
  library(dplyr)
  library(tidyr)
  
  ##----------------------------------------------------------
  ## Leer BIOM
  ##----------------------------------------------------------
  
  biom <- read_biom(biom_file)
  
  otu_mat <- as.matrix(
    biom_data(biom)
  )
  
  OTU <- otu_table(
    otu_mat,
    taxa_are_rows = TRUE
  )
  
  ##----------------------------------------------------------
  ## Read taxonomy
  ##----------------------------------------------------------
  
  tax <- read.delim(
    taxonomy_file,
    sep="\t",
    stringsAsFactors = FALSE
  )
  
  tax_split <- tax %>%
    
    separate(
      Taxon,
      into=c(
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      ),
      sep=";",
      fill="right"
    ) %>%
    
    mutate(
      across(
        Kingdom:Species,
        ~gsub("^[a-z]__", "", .)
      )
    )
  
  rownames(tax_split) <-
    tax_split$Feature.ID
  
  tax_mat <-
    
    as.matrix(
      
      tax_split[,
                c(
                  "Kingdom",
                  "Phylum",
                  "Class",
                  "Order",
                  "Family",
                  "Genus",
                  "Species"
                )]
      
    )
  
  TAX <- tax_table(tax_mat)
  
  ##----------------------------------------------------------
  ## Metadata
  ##----------------------------------------------------------
  
  meta <- read.delim(
    metadata_file,
    sep="\t",
    stringsAsFactors = FALSE
  )
  
  rownames(meta) <- meta[,1]
  
  meta <- meta[
    colnames(otu_mat),
  ]
  
  SAM <- sample_data(meta)
  
  ##----------------------------------------------------------
  ## Build the phyloseq object
  ##----------------------------------------------------------
  
  ps <- phyloseq(
    OTU,
    TAX,
    SAM
  )
  
  ##----------------------------------------------------------
  ## Filtrado
  ##----------------------------------------------------------
  
  ps <- subset_taxa(
    ps,
    Kingdom=="Bacteria"
  )
  
  ps <- subset_taxa(
    ps,
    Family!="Mitochondria"
  )
  
  ps <- subset_taxa(
    ps,
    Order!="Chloroplast"
  )
  
  if(remove_unclassified){
    
    ps <- subset_taxa(
      ps,
      !is.na(Genus)
    )
    
    ps <- subset_taxa(
      ps,
      Genus!=""
    )
    
    ps <- subset_taxa(
      ps,
      Genus!="unclassified"
    )
    
  }
  
  ps <- prune_taxa(
    taxa_sums(ps)>0,
    ps
  )
  
  ##----------------------------------------------------------
  ## Metadata limpio
  ##----------------------------------------------------------
  
  meta_df <-
    
    data.frame(
      sample_data(ps)
    )
  
  if(!(group_column %in% colnames(meta_df))){
    
    stop(
      paste(
        "La columna",
        group_column,
        "no existe."
      )
    )
    
  }
  
  groups <-
    
    unique(
      meta_df[[group_column]]
    )
  
  ##----------------------------------------------------------
  ## Objetos por grupo
  ##----------------------------------------------------------
  
  ps_groups <- list()
  
  otu_groups <- list()
  
  relab_groups <- list()
  
  summary_list <- list()
  
  for(g in groups){
    
    samples_keep <-
      
      rownames(meta_df)[
        meta_df[[group_column]]==g
      ]
    
    ps_tmp <-
      
      prune_samples(
        samples_keep,
        ps
      )
    
    ps_tmp <-
      
      prune_taxa(
        taxa_sums(ps_tmp)>0,
        ps_tmp
      )
    
    otu <-
      
      as(
        otu_table(ps_tmp),
        "matrix"
      )
    
    if(!taxa_are_rows(ps_tmp))
      otu <- t(otu)
    
    relab <-
      
      sweep(
        otu,
        2,
        colSums(otu),
        "/"
      )
    
    ps_groups[[as.character(g)]] <- ps_tmp
    
    otu_groups[[as.character(g)]] <- otu
    
    relab_groups[[as.character(g)]] <- relab
    
    tmp <-
      
      data.frame(
        
        FeatureID = rownames(otu),
        
        Prev = rowSums(
          otu>0
        ),
        
        Prev_prop = rowMeans(
          otu>0
        ),
        
        MeanCount = rowMeans(otu),
        
        TotalCount = rowSums(otu),
        
        MeanRel = rowMeans(relab),
        
        stringsAsFactors = FALSE
        
      )
    
    colnames(tmp)[-1] <-
      
      paste0(
        colnames(tmp)[-1],
        "_",
        g
      )
    
    summary_list[[as.character(g)]] <- tmp
    
  }
  
  ##----------------------------------------------------------
  ## Build the ASV master table
  ##----------------------------------------------------------
  
  asv_master <- summary_list[[1]]
  
  if(length(summary_list)>1){
    
    for(i in 2:length(summary_list)){
      
      asv_master <-
        
        full_join(
          asv_master,
          summary_list[[i]],
          by="FeatureID"
        )
      
    }
    
  }
  
  asv_master[
    is.na(asv_master)
  ] <- 0
  
  tax_df <-
    
    as.data.frame(
      tax_table(ps)
    )
  
  tax_df$FeatureID <-
    rownames(tax_df)
  
  asv_master <-
    
    left_join(
      asv_master,
      tax_df,
      by="FeatureID"
    )
  
  ##----------------------------------------------------------
  ## Salida
  ##----------------------------------------------------------
  
  datos <- list(
    
    ps = ps,
    
    ps_groups = ps_groups,
    
    otu_groups = otu_groups,
    
    relab_groups = relab_groups,
    
    metadata = meta_df,
    
    taxonomy = tax_df,
    
    asv_master = asv_master
    
  )
  
  return(datos)
  
}

##Procesar datos de Longo
datos = build_phyloseq_dataset(biom_file = LONGO_BIOM_FILE,
                      metadata_file =  LONGO_METADATA_FILE,
                      taxonomy_file = LONGO_TAXONOMY_FILE,
                      group_column = "Stage")


##SPIEC-EASI
ps_JUV <- datos$ps_groups$Juvenile

ps_ADU <- datos$ps_groups$Adult

set.seed(123)

se_JUV <- spiec.easi(
  ps_JUV,
  method = "mb",
  lambda.min.ratio = 1e-2,
  nlambda = 20,
  pulsar.params = list(
    rep.num = 50
  )
)

set.seed(123)

se_ADU <- spiec.easi(
  ps_ADU,
  method = "mb",
  lambda.min.ratio = 1e-2,
  nlambda = 20,
  pulsar.params = list(
    rep.num = 50
  )
)

#Matrices de adyacencia
adj_JUV <- getRefit(se_JUV)

adj_ADU <- getRefit(se_ADU)

# Build networks
g_JUV <- graph_from_adjacency_matrix(
  adj_JUV,
  mode = "undirected",
  diag = FALSE
)

g_ADU <- graph_from_adjacency_matrix(
  adj_ADU,
  mode = "undirected",
  diag = FALSE
)

rownames(adj_JUV) <- taxa_names(ps_JUV)
colnames(adj_JUV) <- taxa_names(ps_JUV)

rownames(adj_ADU) <- taxa_names(ps_ADU)
colnames(adj_ADU) <- taxa_names(ps_ADU)

g_JUV <- graph_from_adjacency_matrix(
  adj_JUV,
  mode = "undirected",
  diag = FALSE
)

g_ADU <- graph_from_adjacency_matrix(
  adj_ADU,
  mode = "undirected",
  diag = FALSE
)

#Dataframes de atributos
vertex_JUV <-
  
  datos$asv_master %>%
  
  dplyr::select(
    
    FeatureID,
    
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    
    MeanRel_Juvenile,
    Prev_Juvenile,
    MeanCount_Juvenile,
    TotalCount_Juvenile
    
  )

vertex_ADU <-
  
  datos$asv_master %>%
  
  dplyr::select(
    
    FeatureID,
    
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    
    MeanRel_Adult,
    Prev_Adult,
    MeanCount_Adult,
    TotalCount_Adult
    
  )

idx <- match(
  V(g_JUV)$name,
  vertex_JUV$FeatureID
)

idx2 <- match(
  V(g_ADU)$name,
  vertex_ADU$FeatureID
)

V(g_JUV)$Kingdom <- vertex_JUV$Kingdom[idx]
V(g_JUV)$Phylum  <- vertex_JUV$Phylum[idx]
V(g_JUV)$Class   <- vertex_JUV$Class[idx]
V(g_JUV)$Order   <- vertex_JUV$Order[idx]
V(g_JUV)$Family  <- vertex_JUV$Family[idx]
V(g_JUV)$Genus   <- vertex_JUV$Genus[idx]
V(g_JUV)$Species <- vertex_JUV$Species[idx]

V(g_ADU)$Kingdom <- vertex_ADU$Kingdom[idx2]
V(g_ADU)$Phylum  <- vertex_ADU$Phylum[idx2]
V(g_ADU)$Class   <- vertex_ADU$Class[idx2]
V(g_ADU)$Order   <- vertex_ADU$Order[idx2]
V(g_ADU)$Family  <- vertex_ADU$Family[idx2]
V(g_ADU)$Genus   <- vertex_ADU$Genus[idx2]
V(g_ADU)$Species <- vertex_ADU$Species[idx2]

V(g_JUV)$MeanRel <-
  vertex_JUV$MeanRel_Juvenile[idx]

V(g_JUV)$Prev <-
  vertex_JUV$Prev_Juvenile[idx]

V(g_JUV)$MeanCount <-
  vertex_JUV$MeanCount_Juvenile[idx]

V(g_JUV)$TotalCount <-
  vertex_JUV$TotalCount_Juvenile[idx]

V(g_ADU)$MeanRel <-
  vertex_ADU$MeanRel_Adult[idx2]

V(g_ADU)$Prev <-
  vertex_ADU$Prev_Adult[idx2]

V(g_ADU)$MeanCount <-
  vertex_ADU$MeanCount_Adult[idx2]

V(g_ADU)$TotalCount <-
  vertex_ADU$TotalCount_Adult[idx2]



# Export networks
createNetworkFromIgraph(
  g_JUV,
  title = "Juvenile_network",
  collection = "Ontogeny"
)

createNetworkFromIgraph(
  g_ADU,
  title = "Adult_network",
  collection = "Ontogeny"
)

## Import node tables
YL_mods <- read.table(
  "Juvenile_network_module_assignation.csv",
  header = TRUE,
  sep = ";",
  quote = "\"",
  stringsAsFactors = FALSE
)

PAR_mods <- read.table(
  "Adult_network default node.csv",
  header = TRUE,
  sep = ";",
  quote = "\"",
  stringsAsFactors = FALSE
)

# Calculate topological metrics
YL_mods <- YL_mods %>%
  mutate(
    Degree_norm = Degree/(n()-1)
  )
deg95_YL <- quantile(
  YL_mods$Degree_norm,
  .95,
  na.rm = TRUE
)

YL_mods$BetweennessCentrality <-
  
  as.numeric(
    
    gsub(
      ",",
      ".",
      YL_mods$BetweennessCentrality
    )
    
  )

# Define thresholds
bet95_YL <- quantile(
  YL_mods$BetweennessCentrality,
  .95,
  na.rm = TRUE
)

# Define hubs
YL_mods <- YL_mods %>%
  mutate(
    
    Hub_degree =
      Degree_norm >= deg95_YL,
    
    Hub_between =
      BetweennessCentrality >= bet95_YL,
    
    Hub =
      Hub_degree | Hub_between
    
  )
  YL_mods$MeanRel <-
  
  as.numeric(
    
    gsub(
      ",",
      ".",
      YL_mods$MeanRel
    )
    
  )
  
  YL_mods$log10MeanRel <- log10(YL_mods$MeanRel + 1e-6)

cor.test(
  YL_mods$log10MeanRel,
  YL_mods$Degree,
  method = "spearman"
)

cor.test(
  YL_mods$log10MeanRel,
  YL_mods$BetweennessCentrality,
  method = "spearman"
)

degree_set <-
  YL_mods$name[
    YL_mods$Hub_degree
  ]

between_set <-
  YL_mods$name[
    YL_mods$Hub_between
  ]

length(intersect(
  degree_set,
  between_set
))

length(union(
  degree_set,
  between_set
))

length(intersect(
  degree_set,
  between_set
))/
  
  length(union(
    degree_set,
    between_set
  ))

tab_deg <- table(
  YL_mods$X__glayCluster,
  YL_mods$Hub_degree
)

chisq_deg_mc <- chisq.test(
  tab_deg,
  simulate.p.value = TRUE,
  B = 10000
)

chisq_deg_mc

round(chisq_deg_mc$stdres,2)

cramerV(tab_deg)

tab_btw <- table(
  YL_mods$X__glayCluster,
  YL_mods$Hub_between
)

chisq_btw_mc <- chisq.test(
  tab_btw,
  simulate.p.value = TRUE,
  B = 10000
)

chisq_btw_mc

round(chisq_btw$stdres, 2)

cramerV(tab_btw)

boot_spearman <- function(data, indices, x, y){
  
  d <- data[indices, ]
  
  cor(
    d[[x]],
    d[[y]],
    method = "spearman",
    use = "complete.obs"
  )
  
}

set.seed(123)

boot_deg <- boot(
  
  data = YL_mods,
  
  statistic = function(data, indices){
    
    boot_spearman(
      data,
      indices,
      x = "log10MeanRel",
      y = "Degree"
    )
    
  },
  
  R = 1000
  
)

mean(boot_deg$t)

median(boot_deg$t)

sd(boot_deg$t)

quantile(
  boot_deg$t,
  c(.025,.975)
)

set.seed(123)

boot_btw <- boot(
  
  data = YL_mods,
  
  statistic = function(data, indices){
    
    boot_spearman(
      data,
      indices,
      x = "log10MeanRel",
      y = "BetweennessCentrality"
    )
    
  },
  
  R = 1000
  
)

mean(boot_btw$t)

median(boot_btw$t)

sd(boot_btw$t)

quantile(
  boot_btw$t,
  c(.025,.975)
)

subsample_spearman <- function(
    data,
    x,
    y,
    remove = .10,
    n_iter = 1000){
  
  n <- nrow(data)
  
  rho <- numeric(n_iter)
  
  for(i in seq_len(n_iter)){
    
    keep <- sample(
      
      n,
      
      size = round(n*(1-remove)),
      
      replace = FALSE
      
    )
    
    rho[i] <- cor(
      
      data[[x]][keep],
      
      data[[y]][keep],
      
      method = "spearman",
      
      use = "complete.obs"
      
    )
    
  }
  
  rho
  
}

set.seed(123)

sub_deg10 <- subsample_spearman(
  
  data = YL_mods,
  
  x = "log10MeanRel",
  
  y = "Degree",
  
  remove = .10
  
)

sub_deg20 <- subsample_spearman(
  
  YL_mods,
  
  "log10MeanRel",
  
  "Degree",
  
  remove = .20
  
)

sub_deg30 <- subsample_spearman(
  
  YL_mods,
  
  "log10MeanRel",
  
  "Degree",
  
  remove = .30
  
)

set.seed(123)

sub_btw10 <- subsample_spearman(
  
  YL_mods,
  
  "log10MeanRel",
  
  "BetweennessCentrality",
  
  remove = .10
  
)

sub_btw20 <- subsample_spearman(
  
  YL_mods,
  
  "log10MeanRel",
  
  "BetweennessCentrality",
  
  remove = .20
  
)

sub_btw30 <- subsample_spearman(
  
  YL_mods,
  
  "log10MeanRel",
  
  "BetweennessCentrality",
  
  remove = .30
  
)

summarise_sub <- function(x){
  
  tibble(
    
    Mean = mean(x),
    
    SD = sd(x),
    
    Median = median(x),
    
    CI_low = quantile(x,.025),
    
    CI_high = quantile(x,.975)
    
  )
  
}

summarise_sub(sub_deg10)

summarise_sub(sub_deg20)

summarise_sub(sub_deg30)

summarise_sub(sub_btw10)

summarise_sub(sub_btw20)

summarise_sub(sub_btw30)

######## Adults
PAR_mods <- PAR_mods %>%
  mutate(
    Degree_norm = Degree/(n()-1)
  )

deg95 <- quantile(
  PAR_mods$Degree_norm,
  .95,
  na.rm = TRUE
)

PAR_mods$BetweennessCentrality <-
  
  as.numeric(
    
    gsub(
      ",",
      ".",
      PAR_mods$BetweennessCentrality
    )
    
  )

# Define thresholds
bet95 <- quantile(
  PAR_mods$BetweennessCentrality,
  .95,
  na.rm = TRUE
)

# Define hubs
PAR_mods <- PAR_mods %>%
  mutate(
    
    Hub_degree =
      Degree_norm >= deg95,
    
    Hub_between =
      BetweennessCentrality >= bet95,
    
    Hub =
      Hub_degree | Hub_between
    
  )
PAR_mods$MeanRel <-
  
  as.numeric(
    
    gsub(
      ",",
      ".",
      PAR_mods$MeanRel
    )
    
  )

PAR_mods$log10MeanRel <- log10(PAR_mods$MeanRel + 1e-6)

cor.test(
  PAR_mods$log10MeanRel,
  PAR_mods$Degree,
  method = "spearman"
)

cor.test(
  PAR_mods$log10MeanRel,
  PAR_mods$BetweennessCentrality,
  method = "spearman"
)

degree_set <-
  PAR_mods$name[
    PAR_mods$Hub_degree
  ]

between_set <-
  PAR_mods$name[
    PAR_mods$Hub_between
  ]

length(intersect(
  degree_set,
  between_set
))

length(union(
  degree_set,
  between_set
))

length(intersect(
  degree_set,
  between_set
))/
  
  length(union(
    degree_set,
    between_set
  ))

tab_deg <- table(
  PAR_mods$X__glayCluster,
  PAR_mods$Hub_degree
)

chisq_deg_mc <- chisq.test(
  tab_deg,
  simulate.p.value = TRUE,
  B = 10000
)

chisq_deg_mc

round(chisq_deg_mc$stdres,2)

cramerV(tab_deg)

tab_btw <- table(
  PAR_mods$X__glayCluster,
  PAR_mods$Hub_between
)

chisq_btw_mc <- chisq.test(
  tab_btw,
  simulate.p.value = TRUE,
  B = 10000
)

chisq_btw_mc

round(chisq_btw$stdres, 2)

cramerV(tab_btw)

boot_spearman <- function(data, indices, x, y){
  
  d <- data[indices, ]
  
  cor(
    d[[x]],
    d[[y]],
    method = "spearman",
    use = "complete.obs"
  )
  
}

set.seed(123)

boot_deg <- boot(
  
  data = PAR_mods,
  
  statistic = function(data, indices){
    
    boot_spearman(
      data,
      indices,
      x = "log10MeanRel",
      y = "Degree"
    )
    
  },
  
  R = 1000
  
)

mean(boot_deg$t)

median(boot_deg$t)

sd(boot_deg$t)

quantile(
  boot_deg$t,
  c(.025,.975)
)

set.seed(123)

boot_btw <- boot(
  
  data = PAR_mods,
  
  statistic = function(data, indices){
    
    boot_spearman(
      data,
      indices,
      x = "log10MeanRel",
      y = "BetweennessCentrality"
    )
    
  },
  
  R = 1000
  
)

mean(boot_btw$t)

median(boot_btw$t)

sd(boot_btw$t)

quantile(
  boot_btw$t,
  c(.025,.975)
)

subsample_spearman <- function(
    data,
    x,
    y,
    remove = .10,
    n_iter = 1000){
  
  n <- nrow(data)
  
  rho <- numeric(n_iter)
  
  for(i in seq_len(n_iter)){
    
    keep <- sample(
      
      n,
      
      size = round(n*(1-remove)),
      
      replace = FALSE
      
    )
    
    rho[i] <- cor(
      
      data[[x]][keep],
      
      data[[y]][keep],
      
      method = "spearman",
      
      use = "complete.obs"
      
    )
    
  }
  
  rho
  
}

set.seed(123)

sub_deg10 <- subsample_spearman(
  
  data = PAR_mods,
  
  x = "log10MeanRel",
  
  y = "Degree",
  
  remove = .10
  
)

sub_deg20 <- subsample_spearman(
  
  PAR_mods,
  
  "log10MeanRel",
  
  "Degree",
  
  remove = .20
  
)

sub_deg30 <- subsample_spearman(
  
  PAR_mods,
  
  "log10MeanRel",
  
  "Degree",
  
  remove = .30
  
)

set.seed(123)

sub_btw10 <- subsample_spearman(
  
  PAR_mods,
  
  "log10MeanRel",
  
  "BetweennessCentrality",
  
  remove = .10
  
)

sub_btw20 <- subsample_spearman(
  
  PAR_mods,
  
  "log10MeanRel",
  
  "BetweennessCentrality",
  
  remove = .20
  
)

sub_btw30 <- subsample_spearman(
  
  PAR_mods,
  
  "log10MeanRel",
  
  "BetweennessCentrality",
  
  remove = .30
  
)

summarise_sub <- function(x){
  
  tibble(
    
    Mean = mean(x),
    
    SD = sd(x),
    
    Median = median(x),
    
    CI_low = quantile(x,.025),
    
    CI_high = quantile(x,.975)
    
  )
  
}

summarise_sub(sub_deg10)

summarise_sub(sub_deg20)

summarise_sub(sub_deg30)

summarise_sub(sub_btw10)

summarise_sub(sub_btw20)

summarise_sub(sub_btw30)


# ---- 10. Sequence-equivalent ASVs and topological configurations --------------

# Sequence equivalence is defined by identical nucleotide sequences across the
# shared 199-bp region. This is a sequence-level correspondence, not evidence
# that ASVs represent the same strain or a continuously tracked lineage.
#
# Group labels used throughout this section:
# T = tadpole, J = juvenile, A = adult.
#
# The T topological annotations are taken from the 875-node Cytoscape export
# (FAR_mods), not from g_skin, which is a separate core-restricted network.

# ---- 10.1 Import representative sequences ------------------------------------

seq_T <- readDNAStringSet(FAR_SEQS_FILE)
seq_post <- readDNAStringSet(LONGO_SEQS_FILE)

# Retain the shared 199-bp region for sequence matching.
seq_T_199 <- DNAStringSet(
  substr(as.character(seq_T), 1, 199)
)
names(seq_T_199) <- names(seq_T)

seq_post_199 <- DNAStringSet(
  substr(as.character(seq_post), 1, 199)
)
names(seq_post_199) <- names(seq_post)

sequence_length_audit <- tibble(
  Dataset = c("T", "J/A"),
  Min_length = c(min(width(seq_T)), min(width(seq_post))),
  Max_length = c(max(width(seq_T)), max(width(seq_post))),
  Matching_region_bp = 199
)

sequence_length_audit

# Confirm that every ASV in the relevant phyloseq objects has a FASTA record.
stopifnot(
  all(taxa_names(ps_far) %in% names(seq_T_199)),
  all(taxa_names(ps_JUV) %in% names(seq_post_199)),
  all(taxa_names(ps_ADU) %in% names(seq_post_199))
)

# ---- 10.2 Build the ASV-to-sequence mapping -----------------------------------

map_T <- tibble(
  ASV_ID = taxa_names(ps_far),
  Sequence_199 = as.character(
    seq_T_199[match(taxa_names(ps_far), names(seq_T_199))]
  ),
  Group = "T"
)

map_J <- tibble(
  ASV_ID = taxa_names(ps_JUV),
  Sequence_199 = as.character(
    seq_post_199[match(taxa_names(ps_JUV), names(seq_post_199))]
  ),
  Group = "J"
)

map_A <- tibble(
  ASV_ID = taxa_names(ps_ADU),
  Sequence_199 = as.character(
    seq_post_199[match(taxa_names(ps_ADU), names(seq_post_199))]
  ),
  Group = "A"
)

asv_mapping_clean <- bind_rows(
  map_T,
  map_J,
  map_A
)

# Assign one identifier to each unique sequence across the three groups.
sequence_catalog_clean <- asv_mapping_clean %>%
  distinct(Sequence_199) %>%
  mutate(
    Sequence_ID = sprintf("SEQ_%04d", row_number())
  )

asv_mapping_clean <- asv_mapping_clean %>%
  left_join(sequence_catalog_clean, by = "Sequence_199")

# Mapping integrity checks.
mapping_audit <- asv_mapping_clean %>%
  group_by(Group) %>%
  summarise(
    ASVs = n(),
    Unique_ASVs = n_distinct(ASV_ID),
    Missing_sequences = sum(is.na(Sequence_199)),
    Missing_sequence_IDs = sum(is.na(Sequence_ID)),
    Duplicated_ASV_IDs = sum(duplicated(ASV_ID)),
    .groups = "drop"
  )

mapping_audit

stopifnot(
  !anyNA(asv_mapping_clean$Sequence_ID),
  nrow(asv_mapping_clean) == nrow(distinct(asv_mapping_clean, Group, ASV_ID))
)

# ---- 10.3 Determine sequence presence across groups ---------------------------

# Each ASV in the mapping was detected in at least one sample of its group.
# Presence is therefore summarized at the sequence-equivalent level.
presence_groups <- asv_mapping_clean %>%
  distinct(Sequence_ID, Group) %>%
  mutate(Present = TRUE) %>%
  pivot_wider(
    names_from = Group,
    values_from = Present,
    values_fill = FALSE
  ) %>%
  mutate(
    across(
      c(T, J, A),
      ~ replace_na(.x, FALSE)
    ),
    Pattern = paste0(
      if_else(T, "T", ""),
      if_else(J, "J", ""),
      if_else(A, "A", "")
    )
  ) %>%
  select(Sequence_ID, T, J, A, Pattern)

sequence_pattern_summary <- presence_groups %>%
  count(Pattern, name = "N_sequences") %>%
  mutate(
    Percentage = round(100 * N_sequences / sum(N_sequences), 2)
  ) %>%
  arrange(desc(N_sequences))

sequence_pattern_summary

stopifnot(
  nrow(presence_groups) == n_distinct(asv_mapping_clean$Sequence_ID),
  !anyDuplicated(presence_groups$Sequence_ID),
  !anyNA(presence_groups$Pattern)
)

# ---- 10.4 Harmonize hub annotations from Cytoscape ----------------------------

# FAR_mods contains the 875-node T network annotation table.
# YL_mods and PAR_mods contain the J and A network annotation tables.
network_hubs <- bind_rows(
  FAR_mods %>%
    transmute(
      Group = "T",
      ASV_ID = as.character(name),
      Hub_degree = as.logical(as.character(Hub_degree)),
      Hub_between = as.logical(as.character(Hub_between))
    ),
  YL_mods %>%
    transmute(
      Group = "J",
      ASV_ID = as.character(name),
      Hub_degree = as.logical(as.character(Hub_degree)),
      Hub_between = as.logical(as.character(Hub_between))
    ),
  PAR_mods %>%
    transmute(
      Group = "A",
      ASV_ID = as.character(name),
      Hub_degree = as.logical(as.character(Hub_degree)),
      Hub_between = as.logical(as.character(Hub_between))
    )
) %>%
  mutate(
    Hub_status = case_when(
      Hub_degree & Hub_between ~ "Dual_hub",
      Hub_degree ~ "Degree_hub",
      Hub_between ~ "Betweenness_hub",
      TRUE ~ "Non_hub"
    )
  )

network_hubs_audit <- network_hubs %>%
  group_by(Group) %>%
  summarise(
    Nodes = n(),
    Unique_ASVs = n_distinct(ASV_ID),
    Missing_ASV_IDs = sum(is.na(ASV_ID) | ASV_ID == ""),
    Missing_degree_annotations = sum(is.na(Hub_degree)),
    Missing_betweenness_annotations = sum(is.na(Hub_between)),
    Missing_hub_status = sum(is.na(Hub_status)),
    .groups = "drop"
  )

network_hubs_audit

stopifnot(
  !anyDuplicated(network_hubs[c("Group", "ASV_ID")]),
  !anyNA(network_hubs$Hub_status)
)

# ---- 10.5 Assign topological states to sequence-equivalent units ---------------

sequence_hub_status <- asv_mapping_clean %>%
  left_join(
    network_hubs,
    by = c("Group", "ASV_ID")
  )

# All mapped ASVs should have a corresponding Cytoscape annotation.
sequence_annotation_audit <- sequence_hub_status %>%
  group_by(Group) %>%
  summarise(
    ASVs = n(),
    Missing_network_annotations = sum(is.na(Hub_status)),
    .groups = "drop"
  )

sequence_annotation_audit

stopifnot(!anyNA(sequence_hub_status$Hub_status))

# A sequence may map to multiple ASVs within a group. If those ASVs have
# different hub states, the sequence-equivalent unit is classified as Discordant.
sequence_hub_summary <- sequence_hub_status %>%
  group_by(Sequence_ID, Group) %>%
  summarise(
    Hub_status = {
      states <- unique(Hub_status[!is.na(Hub_status)])
      if (length(states) == 0) {
        "Not_in_network"
      } else if (length(states) > 1) {
        "Discordant"
      } else {
        states[[1]]
      }
    },
    ASVs_per_sequence = n_distinct(ASV_ID),
    .groups = "drop"
  )

# Build the complete descriptive matrix, retaining sequences detected in
# only one group and explicitly marking groups in which a sequence was absent.
sequence_topology_complete <- sequence_hub_summary %>%
  select(Sequence_ID, Group, Hub_status) %>%
  pivot_wider(
    names_from = Group,
    values_from = Hub_status
  ) %>%
  left_join(
    presence_groups %>%
      select(Sequence_ID, Pattern),
    by = "Sequence_ID"
  ) %>%
  mutate(
    across(
      c(T, J, A),
      ~ replace_na(.x, "Not_detected")
    )
  ) %>%
  select(Sequence_ID, Pattern, T, J, A)

# Matrix integrity and sequence-presence summaries.
dim(sequence_topology_complete)

sequence_topology_complete %>%
  count(Pattern, name = "N_sequences")

sequence_topology_complete %>%
  summarise(
    Missing_T = sum(is.na(T)),
    Missing_J = sum(is.na(J)),
    Missing_A = sum(is.na(A)),
    Missing_Pattern = sum(is.na(Pattern))
  )

stopifnot(
  nrow(sequence_topology_complete) == nrow(sequence_catalog_clean),
  !anyDuplicated(sequence_topology_complete$Sequence_ID),
  !anyNA(sequence_topology_complete[c("T", "J", "A", "Pattern")])
)

# Descriptive distribution of topological states within each group.
topology_state_summary <- sequence_topology_complete %>%
  select(T, J, A) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Group",
    values_to = "Hub_status"
  ) %>%
  count(Group, Hub_status, name = "N_sequences") %>%
  arrange(Group, Hub_status)

topology_state_summary

# Identify sequence-equivalent units with conflicting ASV-level hub states.
discordant_sequences <- sequence_topology_complete %>%
  filter(
    T == "Discordant" |
      J == "Discordant" |
      A == "Discordant"
  ) %>%
  select(Sequence_ID, Pattern, T, J, A)

discordant_sequence_details <- sequence_hub_status %>%
  semi_join(
    discordant_sequences,
    by = "Sequence_ID"
  ) %>%
  select(Sequence_ID, Group, ASV_ID, Hub_status) %>%
  arrange(Sequence_ID, Group, Hub_status, ASV_ID)

discordant_sequences
discordant_sequence_details

# Quantify sequence-equivalent units represented by multiple ASVs within
# each group, distinguishing concordant from discordant hub assignments.
sequence_multiplicity <- sequence_hub_status %>%
  group_by(Group, Sequence_ID) %>%
  summarise(
    n_ASVs = n_distinct(ASV_ID),
    n_states = n_distinct(Hub_status),
    Hub_states = paste(
      sort(unique(Hub_status)),
      collapse = " + "
    ),
    .groups = "drop"
  )

sequence_multiplicity_audit <- sequence_multiplicity %>%
  filter(n_ASVs > 1) %>%
  group_by(Group) %>%
  summarise(
    Multi_ASV_sequences = n(),
    Concordant = sum(n_states == 1),
    Discordant = sum(n_states > 1),
    .groups = "drop"
  )

sequence_multiplicity_audit

# ---- 10.6 Summarize shared sequence-equivalent units for Table 2 ---------------

# Only sequence-equivalent units detected in at least two groups are included.
# The table reports observed configurations; no statistical comparisons
# between independently inferred group networks are performed.
shared_topology_summary <- sequence_topology_complete %>%
  filter(Pattern %in% c("TJ", "JA", "TJA")) %>%
  count(Pattern, T, J, A, name = "n") %>%
  arrange(
    factor(Pattern, levels = c("JA", "TJ", "TJA")),
    desc(n),
    T, J, A
  )

shared_topology_summary

shared_sequence_summary <- sequence_topology_complete %>%
  filter(Pattern %in% c("TJ", "JA", "TJA")) %>%
  count(Pattern, name = "N_shared_sequences") %>%
  arrange(factor(Pattern, levels = c("JA", "TJ", "TJA")))

shared_sequence_summary

stopifnot(
  sum(shared_topology_summary$n) ==
    sum(sequence_topology_complete$Pattern %in% c("TJ", "JA", "TJA"))
)

# ---- 10.7 Export sequence-equivalence results ----------------------------------

write.table(
  asv_mapping_clean,
  file.path(RESULTS_DIR, "sequence_equivalence_ASV_mapping.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  sequence_topology_complete,
  file.path(RESULTS_DIR, "sequence_topology_complete.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  shared_topology_summary,
  file.path(RESULTS_DIR, "Table_2_topological_configurations.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  discordant_sequence_details,
  file.path(RESULTS_DIR, "discordant_sequence_ASVs.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  sequence_multiplicity_audit,
  file.path(RESULTS_DIR, "sequence_multiplicity_audit.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---- 11. Sensitivity analysis: homogeneous sample size (n = 4) --------------
# Five independent re-inferences are generated. Each iteration randomly selects
# four of the five T replicates and four of the ten J replicates. The A dataset
# contains four replicates and is retained without subsampling.
#
# This block evaluates the sensitivity of network inference to sample size.
# It does not perform inferential comparisons between independently inferred
# group networks.

infer_network_n4 <- function(ps_obj, stage = c("tadpole", "juvenile", "adult"), seed = NULL) {
  stage <- match.arg(stage)
  if (!is.null(seed)) set.seed(seed)

  if (stage == "tadpole") {
    ps_obj <- subset_samples(ps_obj, environment == "skin")
    ps_obj <- subset_taxa(ps_obj, taxa_names(ps_obj) %in% c(core_skin, skin_only_asvs))
    keep <- sample(sample_names(ps_obj), size = 4, replace = FALSE)
    ps_obj <- prune_samples(keep, ps_obj)
    ps_obj <- filter_taxa(ps_obj, function(x) sum(x > 0) >= 2, TRUE)
    se <- spiec.easi(
      ps_obj,
      method = "mb",
      lambda.min.ratio = 1e-2,
      nlambda = 20,
      pulsar.params = list(thresh = 0.05)
    )
  } else if (stage == "juvenile") {
    keep <- sample(sample_names(ps_obj), size = 4, replace = FALSE)
    ps_obj <- prune_samples(keep, ps_obj)
    se <- spiec.easi(
      ps_obj,
      method = "mb",
      lambda.min.ratio = 1e-2,
      nlambda = 20,
      pulsar.params = list(rep.num = 50)
    )
  } else {
    # The A dataset already contains n = 4.
    se <- spiec.easi(
      ps_obj,
      method = "mb",
      lambda.min.ratio = 1e-2,
      nlambda = 20,
      pulsar.params = list(rep.num = 50)
    )
    keep <- sample_names(ps_obj)
  }

  adj <- as.matrix(getRefit(se))
  taxa <- taxa_names(ps_obj)
  rownames(adj) <- taxa
  colnames(adj) <- taxa

  g <- graph_from_adjacency_matrix(
    adj,
    mode = "undirected",
    diag = FALSE
  )

  topology <- tibble(
    FeatureID = V(g)$name,
    Degree = degree(g),
    BetweennessCentrality = betweenness(g, normalized = TRUE)
  )

  list(
    stage = stage,
    samples = keep,
    phyloseq = ps_obj,
    spiec = se,
    adjacency = adj,
    graph = g,
    topology = topology,
    metrics = data.frame(
      Nodes = vcount(g),
      Edges = ecount(g),
      Connected_components = components(g)$no,
      Average_neighbors = mean(degree(g)),
      Network_density = edge_density(g),
      Clustering_coefficient = transitivity(g, type = "average"),
      Characteristic_path_length = suppressWarnings(
        mean_distance(g, directed = FALSE, unconnected = TRUE)
      )
    )
  )
}

set.seed(123)

sensitivity_n4 <- vector("list", 5)

for (i in seq_len(5)) {
  message("Sensitivity iteration ", i, "/5")

  sensitivity_n4[[i]] <- list(
    iteration = i,
    tadpole = infer_network_n4(ps_far, "tadpole", seed = 100 + i),
    juvenile = infer_network_n4(ps_JUV, "juvenile", seed = 200 + i),
    adult = infer_network_n4(ps_ADU, "adult", seed = 300 + i)
  )
}

sensitivity_metrics_n4 <- purrr::map_dfr(
  sensitivity_n4,
  function(x) {
    bind_rows(
      T = x$tadpole$metrics,
      J = x$juvenile$metrics,
      A = x$adult$metrics,
      .id = "Stage"
    ) %>%
      mutate(Iteration = x$iteration, .before = 1)
  }
)

sensitivity_samples_n4 <- purrr::map_dfr(
  sensitivity_n4,
  function(x) {
    bind_rows(
      T = tibble(Sample = x$tadpole$samples),
      J = tibble(Sample = x$juvenile$samples),
      A = tibble(Sample = x$adult$samples),
      .id = "Stage"
    ) %>%
      mutate(Iteration = x$iteration, .before = 1)
  }
)

write.table(
  sensitivity_metrics_n4,
  file.path(RESULTS_DIR, "sensitivity_metrics_n4.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  sensitivity_samples_n4,
  file.path(RESULTS_DIR, "sensitivity_samples_n4.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# NOTE: Table 1 network metrics and GLay module assignments were obtained in Cytoscape.
# The corresponding exported node/module tables are treated as analysis inputs.

# ---- 12. Save analysis workspace ---------------------------------------------
# The figures script loads this file so that statistical calculations are not
# silently repeated during figure generation.
save.image(
  file = file.path(RESULTS_DIR, "analysis_workspace.RData")
)

message("Analysis pipeline completed. Workspace saved to results/analysis_workspace.RData")
