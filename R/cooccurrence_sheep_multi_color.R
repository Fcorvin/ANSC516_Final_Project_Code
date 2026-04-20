##############################################################
# title: "Microbial Co-occurrence Network- multi color"
# author: "Francis Corvin & Gemini AI"
##############################################################

# 1. Load libraries
library(tidyverse)
library(phyloseq)
library(igraph)
library(ggraph)

# 2. Re-prepare Data
ps_filt <- filter_taxa(ps, function(x) sum(x > 10) > (0.1*length(x)), TRUE)
tax_table(ps_filt)[, "Genus"][is.na(tax_table(ps_filt)[, "Genus"])] <- "Unclassified"
taxa_names(ps_filt) <- make.unique(as.character(tax_table(ps_filt)[, "Genus"]))

# 3. Correlation
cor_mat <- cor(t(otu_table(ps_filt)), method="spearman")
cor_mat[abs(cor_mat) < 0.7] <- 0 
diag(cor_mat) <- 0

# 4. Build Graph 
graph <- graph_from_adjacency_matrix(cor_mat, mode="undirected", weighted=TRUE)

# Fix edge attributes as factors for clean mapping
E(graph)$abs_weight <- abs(E(graph)$weight)
E(graph)$direction <- factor(ifelse(E(graph)$weight > 0, "Positive", "Negative"), 
                             levels = c("Positive", "Negative"))

graph <- delete_vertices(graph, V(graph)[degree(graph) == 0])

# 5. Create Plot Object
network_plot_colored <- ggraph(graph, layout = 'fr', weights = abs_weight) +
  geom_edge_link(aes(edge_alpha = abs_weight, 
                     edge_width = abs_weight, 
                     color = direction)) + 
  geom_node_point(aes(color = name), size = 5) + 
  geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "bold") +
  
  scale_edge_color_manual(values = c("Positive" = "skyblue3", "Negative" = "firebrick"), 
                          name = "Correlation") +
  scale_color_viridis_d(option = "turbo") + 
  
  theme_void() +
  labs(title = "Genus-Level Co-occurrence Network",
       subtitle = "Nodes colored by Genus; Edges = Spearman |rho| > 0.7",
       color = "Genus") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        plot.margin = ggplot2::margin(20, 20, 20, 20))

# 6. SAVE DIRECTLY TO PDF
# We use the pdf() function here because it is immune to RStudio window size issues
if(!dir.exists("output/network/")) dir.create("output/network/", recursive=TRUE)

pdf("output/network/Microbial_Network_By_Genus.pdf", width = 14, height = 10)
print(network_plot_colored)
dev.off()

message("Success! Check output/network/ for your Microbial_Network_By_Genus.pdf")

