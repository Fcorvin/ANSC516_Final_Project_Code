##############################################################
# title: "Microbial Co-occurrence Network"
# author: "Francis Corvin & Gemini AI"
##############################################################

# 1. Clear any active or broken plots
graphics.off() 
library(tidyverse)
library(phyloseq)
library(igraph)
library(ggraph)

# 2. Re-prepare the graph data (to ensure clean attributes)
ps_filt <- filter_taxa(ps, function(x) sum(x > 10) > (0.1*length(x)), TRUE)
tax_table(ps_filt)[, "Genus"][is.na(tax_table(ps_filt)[, "Genus"])] <- "Unclassified"
taxa_names(ps_filt) <- make.unique(as.character(tax_table(ps_filt)[, "Genus"]))

cor_mat <- cor(t(otu_table(ps_filt)), method="spearman")
cor_mat[abs(cor_mat) < 0.7] <- 0 
diag(cor_mat) <- 0

graph <- graph_from_adjacency_matrix(cor_mat, mode="undirected", weighted=TRUE)

# 3. Explicitly set edge attributes
# We use the 'direction' column for the edge colors
E(graph)$abs_weight <- abs(E(graph)$weight)
E(graph)$direction <- ifelse(E(graph)$weight > 0, "Positive", "Negative")
graph <- delete_vertices(graph, V(graph)[degree(graph) == 0])

# 4. The Plot
# NOTE: We use 'edge_colour = direction' inside aes() 
# and 'scale_edge_colour_manual' to match it.
network_plot_colored <- ggraph(graph, layout = 'fr', weights = abs_weight) +
  geom_edge_link(aes(edge_alpha = abs_weight, 
                     edge_width = abs_weight, 
                     edge_colour = direction)) + 
  geom_node_point(aes(color = name), size = 5) + 
  geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "bold") +
  
  # Match the aesthetic 'edge_colour'
  scale_edge_colour_manual(values = c("Positive" = "skyblue3", "Negative" = "firebrick"), 
                           name = "Correlation") +
  scale_color_viridis_d(option = "turbo") + 
  
  theme_void() +
  labs(title = "Genus-Level Co-occurrence Network",
       subtitle = "Nodes colored by Genus; Edges = Spearman |rho| > 0.7",
       color = "Genus") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = ggplot2::margin(20, 20, 20, 20))

# 5. Save the plot
# Using ggsave avoids the viewport errors associated with manual pdf() calls
if(!dir.exists("output/network/")) dir.create("output/network/", recursive=TRUE)

ggsave("output/network/Microbial_Network_By_Genus.pdf", 
       plot = network_plot_colored, 
       width = 14, height = 10, 
       device = "pdf",
       useDingbats = FALSE)

message("Success! The file should now be in output/network/")