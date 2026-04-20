##############################################################
# title: "Microbial heatmap- multi color"
# author: "Francis Corvin & Gemini AI"
##############################################################

# 1. Reset and Load
graphics.off()
library(tidyverse)
library(phyloseq)
library(igraph)
library(ggraph)

# 2. Filter strictly for the "Top Players"
# We'll take the top 30 taxa by total abundance to keep the plot clean
top_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:30])
ps_top <- prune_taxa(top_taxa, ps)

# 3. Simplify Names
tax_table(ps_top)[, "Genus"][is.na(tax_table(ps_top)[, "Genus"])] <- "Unclassified"
taxa_names(ps_top) <- make.unique(as.character(tax_table(ps_top)[, "Genus"]))

# 4. Correlation & Graph
cor_mat <- cor(t(otu_table(ps_top)), method="spearman")
cor_mat[abs(cor_mat) < 0.6] <- 0 # Slightly lower threshold since we have fewer taxa
diag(cor_mat) <- 0

graph <- graph_from_adjacency_matrix(cor_mat, mode="undirected", weighted=TRUE)
E(graph)$direction <- ifelse(E(graph)$weight > 0, "Positive", "Negative")
graph <- delete_vertices(graph, V(graph)[degree(graph) == 0])

# 5. Simple Circle Plot (This avoids the 'physics' crash)
# A circular layout is much more stable than the 'fr' layout
network_simple <- ggraph(graph, layout = 'linear', circular = TRUE) +
  geom_edge_arc(aes(color = direction, alpha = abs(weight), width = abs(weight))) +
  geom_node_point(aes(color = name), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold") +
  scale_edge_color_manual(values = c("Positive" = "skyblue3", "Negative" = "firebrick")) +
  scale_color_viridis_d(option = "turbo") +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "Keystone Microbial Interactions",
       subtitle = "Top 30 Genus-level taxa by abundance")

# 6. Save using a standard GGSAVE (No PDF window needed)
ggsave("output/network/Keystone_Network_Circular.pdf", network_simple, width = 12, height = 12)

