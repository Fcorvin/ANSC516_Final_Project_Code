##############################################################
# title: "Random Forest - Predictive Taxa"
# author: "Francis Corvin & Gemini AI"
##############################################################

library(tidyverse)
library(qiime2R)
library(randomForest)

# 1. Define Paths
metadata_path <- "sheep_project/trial_1/unique_metadata.tsv"
table_path    <- "sheep_project/trial_1/table_T1.qza"
taxonomy_path <- "sheep_project/trial_1/taxonomy.qza"

# 2. Load Data
meta <- read_q2metadata(metadata_path)
table <- read_qza(table_path)$data
taxonomy <- parse_taxonomy(read_qza(taxonomy_path)$data)

# 3. Align Samples
common_samples <- intersect(meta$SampleID, colnames(table))
meta_filtered <- meta %>% filter(SampleID %in% common_samples)
table_filtered <- table[, meta_filtered$SampleID]

# 4. Prepare Variables
X_data <- as.data.frame(t(table_filtered))
Y_groups <- factor(meta_filtered$disease_state, 
                   levels = c("A_HEALTHY", "TRANSITION", "C_ID", "D_FOOTROT", "E_CODD", "TREATED"))

# 5. Run Random Forest
set.seed(42)
rf_model <- randomForest(x = X_data, y = Y_groups, importance = TRUE, ntree = 500)

# --- NEW STEP: SAVE THE STATS ---
if(!dir.exists("output/random_forest")) dir.create("output/random_forest", recursive = TRUE)

# Save the accuracy/confusion matrix for your Results section
write.csv(rf_model$confusion, "output/random_forest/rf_confusion_matrix.csv")
# --------------------------------

# 6. Extract Top 25 Predictive ASVs
importance_df <- as.data.frame(importance(rf_model)) %>%
  rownames_to_column("ASV") %>%
  left_join(taxonomy %>% rownames_to_column("ASV"), by="ASV") %>%
  mutate(Label = paste0(Genus, " (", substr(ASV, 1, 4), "...)")) %>%
  top_n(25, MeanDecreaseGini)

# 7. Professional Visualization
rf_plot <- ggplot(importance_df, aes(x=reorder(Label, MeanDecreaseGini), y=MeanDecreaseGini, color=Phylum)) +
  geom_segment(aes(x=reorder(Label, MeanDecreaseGini), xend=reorder(Label, MeanDecreaseGini), y=0, yend=MeanDecreaseGini), color="grey") +
  geom_point(size=4) +
  coord_flip() + 
  theme_bw() +
  scale_color_viridis_d(option="mako") +
  labs(title="Top Microbial Predictors of Disease State",
       subtitle="Gini Importance: Contribution to model classification accuracy",
       x="", y="Mean Decrease Gini")

ggsave("output/random_forest/figure_RF_top25.pdf", rf_plot, height=8, width=10)

message("Success! You now have the plot AND the accuracy matrix in output/random_forest/")
