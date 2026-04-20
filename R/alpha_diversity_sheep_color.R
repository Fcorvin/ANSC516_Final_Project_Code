##############################################################
# title: "Alpha diversity in R - Sheep Footrot Project (Enhanced Visuals)"
# author: "ANSC595 / Gemini Adaptive AI"
##############################################################

# Set your working directory if needed
# setwd("~/ANSC_516_Rstudio")

## Load Packages
library(tidyverse)
library(qiime2R)
library(ggpubr)

if(!dir.exists("output"))
  dir.create("output")

## Load metadata
meta <- read_q2metadata("sheep_project/trial_1/unique_metadata.tsv")

## Load alpha diversity vectors
evenness <- read_qza("sheep_project/trial_1/swab_core_metrics/evenness_vector.qza")
evenness <- evenness$data %>% rownames_to_column("SampleID")

observed_features <- read_qza("sheep_project/trial_1/swab_core_metrics/observed_features_vector.qza")
observed_features <- observed_features$data %>% rownames_to_column("SampleID")

shannon <- read_qza("sheep_project/trial_1/swab_core_metrics/shannon_vector.qza")
shannon <- shannon$data %>% rownames_to_column("SampleID")

faith_pd <- read_qza("sheep_project/trial_1/swab_core_metrics/faith_pd_vector.qza")
faith_pd <- faith_pd$data %>% rownames_to_column("SampleID")

## Merge all metrics and metadata
alpha_diversity <- merge(x=faith_pd, y=evenness, by = "SampleID")
alpha_diversity <- merge(alpha_diversity, observed_features, by = "SampleID")
alpha_diversity <- merge(alpha_diversity, shannon, by = "SampleID")

# Rename evenness for consistency with ANOVA code
if("pielou_evenness" %in% colnames(alpha_diversity)){
  alpha_diversity <- alpha_diversity %>% rename(pielou_e = pielou_evenness)
}

meta <- merge(meta, alpha_diversity, by = "SampleID")
row.names(meta) <- meta$SampleID

# Set factor order for disease progression
meta$disease_state.ord <- factor(meta$disease_state, 
                                 levels = c("A_HEALTHY", "TRANSITION", "C_ID", 
                                            "D_FOOTROT", "E_CODD", "TREATED"))

##############################################################
# 1. PIELOU'S EVENNESS (Normally Distributed)
##############################################################

# Stats
aov.evenness.disease <- aov(pielou_e ~ disease_state.ord, data=meta)
summary(aov.evenness.disease)
TukeyHSD(aov.evenness.disease)

# Plot
evenness_boxplot <- ggplot(meta, aes(x = disease_state.ord, y = pielou_e, fill = disease_state.ord)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + 
  geom_jitter(width = 0.15, alpha = 0.4) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(y="Pielou's Evenness", x="Disease State", title="Community Evenness by Disease State")

ggsave("output/evenness_color.png", evenness_boxplot, height=4, width=6)

##############################################################
# 2. FAITH'S PD (Non-Normal)
##############################################################

# Stats
kruskal.test(faith_pd ~ disease_state.ord, data=meta)
pairwise.wilcox.test(meta$faith_pd, meta$disease_state.ord, p.adjust.method="BH")

# Plot
faith_pd_boxplot <- ggplot(meta, aes(x = disease_state.ord, y = faith_pd, fill = disease_state.ord)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  scale_fill_brewer(palette = "Set3") + # Using a different palette for variety
  labs(y="Faith Phylogenetic Diversity", x="Disease State", title="Phylogenetic Diversity by Disease State")

ggsave("output/faith_pd_color.png", faith_pd_boxplot, height=4, width=6)

##############################################################
# 3. SHANNON ENTROPY (Non-Normal)
##############################################################

# Stats
kruskal.test(shannon_entropy ~ disease_state.ord, data=meta)
pairwise.wilcox.test(meta$shannon_entropy, meta$disease_state.ord, p.adjust.method="BH")

# Plot
shannon_boxplot <- ggplot(meta, aes(x = disease_state.ord, y = shannon_entropy, fill = disease_state.ord)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  scale_fill_viridis_d(option = "mako") + # Modern dark-to-light palette
  labs(y="Shannon Entropy", x="Disease State", title="Shannon Diversity by Disease State")

ggsave("output/shannon_color.png", shannon_boxplot, height=4, width=6)

##############################################################
# 4. OBSERVED FEATURES (Non-Normal)
##############################################################

# Stats
kruskal.test(observed_features ~ disease_state.ord, data=meta)
pairwise.wilcox.test(meta$observed_features, meta$disease_state.ord, p.adjust.method="BH")

# Plot
observed_boxplot <- ggplot(meta, aes(x = disease_state.ord, y = observed_features, fill = disease_state.ord)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(y="Observed ASVs", x="Disease State", title="Richness (Observed ASVs) by Disease State")

ggsave("output/observed_features_color.png", observed_boxplot, height=4, width=6)

##############################################################
# 5. CONTINUOUS TREND: EVENNESS OVER TIME
##############################################################

glm.evenness.day <- glm(pielou_e ~ day, data=meta)

evenness_time_plot <- ggplot(meta, aes(x = day, y = pielou_e)) +
  geom_point(aes(color = disease_state.ord), size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  theme_q2r() +
  labs(y="Pielou's Evenness", x="Day", color="Disease State", 
       title="Trend of Evenness over Time")

ggsave("output/evenness_time_trend.png", evenness_time_plot, height=4, width=6)