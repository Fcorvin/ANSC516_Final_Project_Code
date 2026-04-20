##############################################################
# title: "Microbiome Synthesis - Beta Div, Composition & Heatmap"
# author: "Francis Corvin & Gemini AI"
##############################################################

# 1. PCoA PLOT (Beta Diversity - Bray-Curtis)
# This shows how the overall community shifts between disease states
dist_bray <- distance(ps, method="bray")
pcoa_res <- ordinate(ps, method="PCoA", distance=dist_bray)

p_pcoa <- plot_ordination(ps, pcoa_res, color="disease_state") +
  geom_point(size=4, alpha=0.8) +
  stat_ellipse(aes(group=disease_state), linetype = 2) +
  scale_color_manual(values=disease_colors) +
  theme_bw() +
  labs(title = "PCoA: Global Community Shift (Bray-Curtis)",
       subtitle = "Ellipses represent 95% confidence intervals") +
  theme(legend.position = "right", plot.title = element_text(face="bold"))

# 2. STACKED BAR PLOTS (Phylum Level Composition)
# This shows the "Ratio" of the main bacterial groups
ps_phylum <- tax_glom(ps, taxrank = "Phylum")
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

p_bar <- plot_bar(ps_phylum_rel, x="disease_state", fill="Phylum") +
  geom_bar(stat="identity", position="stack") +
  scale_fill_viridis_d(option = "mako") + 
  facet_grid(~disease_state, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(y = "Relative Abundance", x = "", title = "Community Composition by Phylum") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# 3. TOP 20 GENUS HEATMAP
# This creates a "signature" barcode of the most abundant players
top20_genera <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:20])
ps_top20 <- prune_taxa(top20_genera, ps)
ps_top20_rel <- transform_sample_counts(ps_top20, function(x) x / sum(x))

# Melt for ggplot heatmap
melted_heatmap <- psmelt(ps_top20_rel)
melted_heatmap$disease_state.ord <- factor(melted_heatmap$disease_state, 
                                           levels = c("A_HEALTHY", "TRANSITION", "C_ID", "D_FOOTROT", "E_CODD", "TREATED"))

p_heatmap <- ggplot(melted_heatmap, aes(x = Sample, y = Genus, fill = Abundance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue") +
  facet_grid(~disease_state.ord, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(title = "Heatmap: Top 20 Genera Abundance", x = "Samples grouped by Disease State") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "italic", size = 8))

# 4. SAVE ALL TO PDF
if(!dir.exists("output/synthesis")) dir.create("output/synthesis", recursive = TRUE)

ggsave("output/synthesis/PCoA_Beta_Diversity.pdf", p_pcoa, width=8, height=6)
ggsave("output/synthesis/Taxa_Barplots_Phylum.pdf", p_bar, width=10, height=6)
ggsave("output/synthesis/Taxa_Heatmap_Top20.pdf", p_heatmap, width=12, height=8)

message("Success! PCoA, Barplots, and Heatmap saved to 'output/synthesis'.")
