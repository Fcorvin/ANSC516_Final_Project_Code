##############################################################
# title: "Beta Diversity Analysis - Sheep CODD Project"
# author: "Francis Corvin & Gemini AI"
##############################################################

## Load Packages
library(tidyverse)
library(qiime2R)
library(vegan)
library(patchwork) # For combining plots

# Create output directory
if(!dir.exists("output")) dir.create("output")

## 1. DATA LOADING
metadata <- read_q2metadata("sheep_project/trial_1/unique_metadata.tsv")
metadata$disease_state.ord <- factor(metadata$disease_state, 
                                     levels = c("A_HEALTHY", "TRANSITION", "C_ID", "D_FOOTROT", "E_CODD", "TREATED"))

# Define professional color palette
disease_colors <- c("A_HEALTHY"  = "#2E8B57", # SeaGreen
                    "TRANSITION" = "#8A2BE2", # BlueViolet
                    "C_ID"       = "#4169E1", # RoyalBlue
                    "D_FOOTROT"  = "#FF8C00", # DarkOrange
                    "E_CODD"     = "#CD5C5C", # IndianRed
                    "TREATED"    = "#2F4F4F") # DarkSlateGray

# Define the metrics and their file prefixes
metrics <- c("jaccard", "bray_curtis", "unweighted_unifrac", "weighted_unifrac")
titles <- c("Jaccard (P/A)", "Bray-Curtis (Abundance)", 
            "Unweighted UniFrac (Phylo P/A)", "Weighted UniFrac (Phylo Abund)")

# List to store plots for combining later
plot_list <- list()

##############################################################
# 2. ANALYSIS LOOP (Ordination & Statistics)
##############################################################

for(i in 1:4){
  m <- metrics[i]
  
  # A. Load PCoA Results
  pcoa <- read_qza(paste0("sheep_project/trial_1/swab_core_metrics/", m, "_pcoa_results.qza"))
  
  # B. Build Plotting Data
  pcoa_df <- pcoa$data$Vectors %>%
    select(SampleID, PC1, PC2) %>%
    inner_join(metadata, by = "SampleID")
  
  # C. Run PERMANOVA (Requirement 3b)
  dist_mat <- read_qza(paste0("sheep_project/trial_1/swab_core_metrics/", m, "_distance_matrix.qza"))
  dm <- as.matrix(dist_mat$data)
  
  # Ensure metadata matches distance matrix order
  meta_sub <- metadata[match(rownames(dm), metadata$SampleID),]
  permanova <- adonis2(as.dist(dm) ~ disease_state.ord, data = meta_sub)
  
  # Save Stats
  write.csv(as.data.frame(permanova), paste0("output/stats_PERMANOVA_", m, ".csv"))
  p_val <- permanova$`Pr(>F)`[1]
  
  # D. Generate Professional Plot
  p <- ggplot(pcoa_df, aes(x=PC1, y=PC2, color=disease_state.ord)) +
    geom_point(aes(shape=foot), size=2.5, alpha=0.8) +
    stat_ellipse(linewidth=0.5, alpha=0.5) +
    theme_q2r() +
    scale_color_manual(values=disease_colors) +
    labs(title = titles[i],
         subtitle = paste0("PERMANOVA p = ", p_val),
         x = paste0("PC1 (", round(100*pcoa$data$ProportionExplained[1], 1), "%)"),
         y = paste0("PC2 (", round(100*pcoa$data$ProportionExplained[2], 1), "%)"),
         color = "Disease State", shape = "Foot") +
    theme(legend.position = "right", 
          plot.title = element_text(face="bold"))
  
  # Save Individual Plot
  ggsave(paste0("output/figure_beta_", m, ".png"), p, height=5, width=7)
  
  # Add to list (removing legend for the combined version)
  plot_list[[i]] <- p + theme(legend.position = "none")
}

##############################################################
# 3. COMBINED MASTER FIGURE
##############################################################
# Uses 'patchwork' to create one high-quality comparison image
combined_plot <- (plot_list[[1]] | plot_list[[2]]) / (plot_list[[3]] | plot_list[[4]]) +
  plot_annotation(title = "Beta Diversity Comparison", 
                  tag_levels = 'A')

ggsave("output/figure_beta_master_panel.pdf", combined_plot, height=10, width=12)

message("Success! 4 PERMANOVA reports and 5 images saved to 'output' folder.")

