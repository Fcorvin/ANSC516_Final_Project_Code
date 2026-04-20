##############################################################
# title: "Professional Violin Plots - Sheep CODD Project"
# author: "Francis Corvin & Gemini AI"
##############################################################


library(qiime2R)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(viridis)

# Create output directory
if(!dir.exists("output/violin")) dir.create("output/violin", recursive = TRUE)

# Professional Disease Progression Palette
disease_colors <- c("A_HEALTHY"  = "#2E8B57", 
                    "TRANSITION" = "#8A2BE2", 
                    "C_ID"       = "#4169E1", 
                    "D_FOOTROT"  = "#FF8C00", 
                    "E_CODD"     = "#CD5C5C", 
                    "TREATED"    = "#2F4F4F")

##############################################################
# 1. LOAD & MERGE DATA (Fixed Column Name Logic)
##############################################################
meta <- read_q2metadata("sheep_project/trial_1/unique_metadata.tsv")
meta$disease_state.ord <- factor(meta$disease_state, 
                                 levels = c("A_HEALTHY", "TRANSITION", "C_ID", "D_FOOTROT", "E_CODD", "TREATED"))

# Helper function to load and rename columns without using tidyverse rename()
load_alpha <- function(path, new_name) {
  df <- read_qza(path)$data %>% rownames_to_column("SampleID")
  colnames(df)[2] <- new_name  # Force the 2nd column to our target name
  return(df)
}

# Load all vectors safely
shannon  <- load_alpha("sheep_project/trial_1/swab_core_metrics/shannon_vector.qza", "shannon_entropy")
faith_pd <- load_alpha("sheep_project/trial_1/swab_core_metrics/faith_pd_vector.qza", "faith_pd")
evenness <- load_alpha("sheep_project/trial_1/swab_core_metrics/evenness_vector.qza", "pielou_e")
observed <- load_alpha("sheep_project/trial_1/swab_core_metrics/observed_features_vector.qza", "observed_features")

# Merge into master plotting dataframe
alpha_data <- meta %>%
  inner_join(shannon, by="SampleID") %>%
  inner_join(faith_pd, by="SampleID") %>%
  inner_join(evenness, by="SampleID") %>%
  inner_join(observed, by="SampleID")

##############################################################
# 2. ALPHA DIVERSITY VIOLIN PLOTS
##############################################################
metrics <- c("shannon_entropy", "observed_features", "pielou_e", "faith_pd")
ylabels <- c("Shannon Entropy", "Observed ASVs", "Pielou's Evenness", "Faith's PD")

for(i in 1:4){
  p <- ggplot(alpha_data, aes(x=disease_state.ord, y=.data[[metrics[i]]], fill=disease_state.ord)) +
    geom_violin(trim=FALSE, alpha=0.6, color="black") +
    geom_boxplot(width=0.1, fill="white", outlier.shape=NA, color="black") +
    geom_jitter(width=0.1, alpha=0.2, size=1) +
    stat_compare_means(method = "kruskal.test", label.y.npc = "top", size=3) +
    theme_q2r() +
    scale_fill_manual(values=disease_colors) +
    labs(y=ylabels[i], x="", title=paste("Microbial Shifts:", ylabels[i])) +
    theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none")
  
  ggsave(paste0("output/violin/violin_alpha_", metrics[i], ".png"), p, height=5, width=7)
}

##############################################################
# 3. TAXA DISTRIBUTION VIOLIN PLOTS
##############################################################
ps <- qza_to_phyloseq(
  features = "sheep_project/trial_1/swab_core_metrics/rarefied_table.qza",
  taxonomy = "sheep_project/trial_1/taxonomy.qza",
  metadata = "sheep_project/trial_1/unique_metadata.tsv"
)

plot_taxa_violin <- function(ps_obj, rank, top_n=6) {
  ps_rel <- transform_sample_counts(ps_obj, function(x) x / sum(x))
  ps_glom <- tax_glom(ps_rel, taxrank = rank)
  
  melted <- psmelt(ps_glom)
  melted$disease_state.ord <- factor(melted$disease_state, 
                                     levels = c("A_HEALTHY", "TRANSITION", "C_ID", "D_FOOTROT", "E_CODD", "TREATED"))
  
  top_list <- melted %>%
    group_by(!!sym(rank)) %>%
    summarise(Mean = mean(Abundance)) %>%
    top_n(top_n, Mean) %>%
    pull(!!sym(rank))
  
  melted_sub <- melted %>% filter(!!sym(rank) %in% top_list)
  
  ggplot(melted_sub, aes(x=disease_state.ord, y=Abundance, fill=disease_state.ord)) +
    geom_violin(trim=FALSE, alpha=0.7) +
    geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
    facet_wrap(as.formula(paste0("~", rank)), scales="free_y", ncol=2) +
    scale_fill_manual(values=disease_colors) +
    theme_q2r() +
    labs(y="Relative Abundance", x="", title=paste("Top", rank, "Shifts Across Disease States")) +
    theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="none",
          strip.text = element_text(face="italic", size=8))
}

p_phylum_violin <- plot_taxa_violin(ps, "Phylum", top_n=4)
ggsave("output/violin/violin_taxa_phylum.pdf", p_phylum_violin, height=8, width=10)

p_genus_violin <- plot_taxa_violin(ps, "Genus", top_n=6)
ggsave("output/violin/violin_taxa_genus.pdf", p_genus_violin, height=10, width=12)

message("Success! All violin plots are now in 'output/violin'.")

