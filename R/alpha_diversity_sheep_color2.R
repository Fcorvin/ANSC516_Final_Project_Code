##############################################################
# title: "Alpha Diversity & Stats - Sheep CODD Project"
# author: "Francis Corvin & Gemini AI"
##############################################################

## Load Packages
library(tidyverse)
library(qiime2R)
library(ggpubr)
library(viridis)

# Create output directory
if(!dir.exists("output")) dir.create("output")

## 1. DATA LOADING & MERGE
meta <- read_q2metadata("sheep_project/trial_1/unique_metadata.tsv")

# Load Diversity Vectors
evenness <- read_qza("sheep_project/trial_1/swab_core_metrics/evenness_vector.qza")$data %>% rownames_to_column("SampleID")
observed <- read_qza("sheep_project/trial_1/swab_core_metrics/observed_features_vector.qza")$data %>% rownames_to_column("SampleID")
shannon <- read_qza("sheep_project/trial_1/swab_core_metrics/shannon_vector.qza")$data %>% rownames_to_column("SampleID")
faith_pd <- read_qza("sheep_project/trial_1/swab_core_metrics/faith_pd_vector.qza")$data %>% rownames_to_column("SampleID")

# Merge metrics and standardize Evenness name
alpha_diversity <- faith_pd %>% 
  inner_join(evenness, by="SampleID") %>% 
  inner_join(observed, by="SampleID") %>% 
  inner_join(shannon, by="SampleID")

if("pielou_evenness" %in% colnames(alpha_diversity)) {
  alpha_diversity <- alpha_diversity %>% rename(pielou_e = pielou_evenness)
}

# Final Merge with metadata
meta <- meta %>% inner_join(alpha_diversity, by="SampleID")
meta$disease_state.ord <- factor(meta$disease_state, 
                                 levels = c("A_HEALTHY", "TRANSITION", "C_ID", "D_FOOTROT", "E_CODD", "TREATED"))

##############################################################
# 2. RAREFACTION CURVE (Requirement 2c)
##############################################################
# Use the exported CSV from Bell
rarefaction_raw <- read.csv("sheep_project/trial_1/swab_core_metrics/rarefaction_data.csv", check.names = FALSE)

# Standardize column names to lowercase
colnames(rarefaction_raw) <- tolower(colnames(rarefaction_raw))

rarefaction_plot <- rarefaction_raw %>%
  # FIX: Pivot ONLY numeric columns, excluding the first two (depth and iter)
  pivot_longer(
    cols = where(is.numeric) & -c(1, 2), 
    names_to = "SampleID", 
    values_to = "richness"
  ) %>%
  rename(depth = 1, iter = 2) %>% 
  filter(SampleID %in% meta$SampleID) %>%
  left_join(meta, by = "SampleID") %>%
  group_by(depth, disease_state.ord) %>%
  summarise(mean_richness = mean(richness, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = depth, y = mean_richness, color = disease_state.ord)) +
  geom_line(linewidth = 1.2) +
  theme_q2r() +
  scale_color_viridis_d() +
  labs(x = "Sequencing Depth", 
       y = "Observed Features", 
       title = "Rarefaction Curves (Justification for 10k Depth)")

ggsave("output/figure_rarefaction.png", rarefaction_plot, height = 5, width = 7)

##############################################################
# 3. BOXPLOTS WITH STATISTICAL TESTING (Requirement 3a)
##############################################################
# Define the 4 metrics for the loop
metrics <- c("shannon_entropy", "observed_features", "pielou_e", "faith_pd")
titles <- c("Shannon Entropy (Diversity)", "Observed ASVs (Richness)", 
            "Pielou's Evenness (Evenness)", "Faith's PD (Phylogeny)")

for(i in 1:4){
  p <- ggplot(meta, aes_string(x="disease_state.ord", y=metrics[i], fill="disease_state.ord")) +
    geom_boxplot(alpha=0.7, outlier.shape=NA) +
    geom_jitter(width=0.1, alpha=0.2) +
    # Requirement 3a: Kruskal-Wallis Test
    stat_compare_means(method = "kruskal.test", label.y.npc = "top") + 
    theme_q2r() +
    theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none") +
    scale_fill_viridis_d(option = "plasma") +
    labs(y=metrics[i], x="", title=titles[i])
  
  ggsave(paste0("output/figure_alpha_", metrics[i], ".png"), p, height=4, width=6)
}

message("Success! Rarefaction and 4 Alpha Diversity plots saved to 'output' folder.")

