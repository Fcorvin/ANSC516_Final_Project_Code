##############################################################
# title: "Complete ANCOM-BC2 Analysis and Visualization"
# author: "Francis Corvin & Gemini AI"
##############################################################

# 1. Load All Required Libraries
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(ANCOMBC)
library(ggplot2)

# 2. Define File Paths and Create Output Directory
# Update these paths if your folder structure is different
input_table <- "sheep_project/trial_1/table_T1.qza"
input_taxonomy <- "sheep_project/trial_1/taxonomy.qza"
input_metadata <- "sheep_project/trial_1/unique_metadata.tsv"
output_dir <- "output/differential_abundance/"

if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

# 3. Load Data into Phyloseq Object
ps <- qza_to_phyloseq(
  features = input_table,
  taxonomy = input_taxonomy,
  metadata = input_metadata
)

# 4. Filter Samples and Setup Reference Group
# We focus on the disease progression and set Healthy as the baseline
ps_sub <- subset_samples(ps, disease_state %in% c("A_HEALTHY", "C_ID", "D_FOOTROT", "E_CODD"))
sample_data(ps_sub)$disease_state <- factor(
  sample_data(ps_sub)$disease_state, 
  levels = c("A_HEALTHY", "C_ID", "D_FOOTROT", "E_CODD")
)

# 5. Run ANCOM-BC2 Statistical Model
print("Running ANCOM-BC2 analysis... this may take a few minutes.")
out = ancombc2(
  data = ps_sub, 
  fix_formula = "disease_state", 
  p_adj_method = "holm", 
  group = "disease_state", 
  struc_zero = TRUE, 
  neg_lb = TRUE
)

# 6. Extract Results into a Data Frame
res_df <- as.data.frame(out$res)

# 7. Process Taxonomy for Labels
tax_raw <- read_qza(input_taxonomy)$data
taxonomy_clean <- parse_taxonomy(tax_raw)
tax_labels <- taxonomy_clean %>% 
  tibble::rownames_to_column("ASV") %>%
  dplyr::mutate(Label = paste0(Genus, " (", substr(ASV, 1, 4), ")"))

# 8. Prepare Data for the Multi-Disease Bar Chart
plot_data <- res_df %>%
  dplyr::rename(ASV = taxon) %>%
  tidyr::pivot_longer(
    cols = c("lfc_disease_stateC_ID", "lfc_disease_stateD_FOOTROT", "lfc_disease_stateE_CODD"), 
    names_to = "Stage", 
    values_to = "Beta"
  ) %>%
  dplyr::mutate(Stage = dplyr::case_when(
    grepl("C_ID", Stage) ~ "Infectious Dermatitis (ID)",
    grepl("D_FOOTROT", Stage) ~ "Footrot",
    grepl("E_CODD", Stage) ~ "CODD"
  )) %>%
  dplyr::left_join(tax_labels, by="ASV")

# 9. Filter for Significant Biomarkers (q < 0.05 and Effect Size > 1.5)
sig_asvs <- res_df %>%
  dplyr::filter(q_disease_stateC_ID < 0.05 | q_disease_stateD_FOOTROT < 0.05 | q_disease_stateE_CODD < 0.05) %>%
  dplyr::pull(taxon)

plot_final_sig <- plot_data %>%
  dplyr::filter(ASV %in% sig_asvs) %>%
  dplyr::group_by(ASV) %>%
  dplyr::filter(any(abs(Beta) > 1.5)) %>% 
  dplyr::ungroup()

# --- GRAPHIC 1: MULTI-DISEASE FACETED BAR CHART ---
biomarker_plot <- ggplot(plot_final_sig, aes(x=reorder(Label, Beta), y=Beta, fill=Phylum)) +
  geom_bar(stat="identity") +
  coord_flip() +
  facet_wrap(~Stage) +
  theme_bw() +
  scale_fill_viridis_d(option="turbo") +
  labs(title="Microbial Biomarkers of Disease Progression",
       subtitle="Significant Taxa (q < 0.05) compared to Healthy Sheep",
       x="", y="Log-Fold Change (Effect Size)") +
  theme(axis.text.y = element_text(size=6, face="italic"),
        legend.position = "bottom")

ggsave(paste0(output_dir, "figure_MultiDisease_Biomarkers.pdf"), biomarker_plot, width=14, height=10)

# --- GRAPHIC 2: VOLCANO PLOT (CODD VS HEALTHY) ---
volcano_data <- res_df %>%
  dplyr::rename(ASV = taxon) %>%
  dplyr::left_join(tax_labels, by="ASV") %>%
  dplyr::mutate(Significant = q_disease_stateE_CODD < 0.05 & abs(lfc_disease_stateE_CODD) > 1.5)

volcano_plot <- ggplot(volcano_data, aes(x=lfc_disease_stateE_CODD, y=-log10(q_disease_stateE_CODD), color=Significant)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("grey70", "red")) +
  theme_bw() +
  labs(title="Volcano Plot: CODD vs Healthy",
       subtitle="Red points indicate q < 0.05 and Fold Change > 1.5",
       x="Log Fold Change", y="-Log10 Q-value")

ggsave(paste0(output_dir, "figure_Volcano_CODD.pdf"), volcano_plot, width=8, height=7)

# 10. Save CSV of raw results for your records
write.csv(res_df, paste0(output_dir, "ANCOMBC_results_final.csv"), row.names=FALSE)

message("Success! Script complete. Check the output/differential_abundance/ folder.")

# --- ADD THIS TO THE END OF YOUR SCRIPT TO GENERATE THE SIMPLIFIED PLOT ---

# 1. Prepare Data with Stricter Filters
# We increase Beta threshold to 2.5 and limit to the Top 15 taxa per stage
plot_data_simplified <- res_df %>%
  dplyr::rename(ASV = taxon) %>%
  tidyr::pivot_longer(
    cols = c("lfc_disease_stateC_ID", "lfc_disease_stateD_FOOTROT", "lfc_disease_stateE_CODD"), 
    names_to = "Stage", 
    values_to = "Beta"
  ) %>%
  dplyr::mutate(Stage = dplyr::case_when(
    grepl("C_ID", Stage) ~ "ID",
    grepl("D_FOOTROT", Stage) ~ "Footrot",
    grepl("E_CODD", Stage) ~ "CODD"
  )) %>%
  dplyr::left_join(tax_labels, by="ASV") %>%
  # Filter for high impact: Log-Fold Change > 2.5
  dplyr::filter(abs(Beta) > 2.5) %>% 
  # Take only the Top 15 per group to ensure the Y-axis is readable
  dplyr::group_by(Stage) %>%
  dplyr::slice_max(order_by = abs(Beta), n = 15) %>% 
  dplyr::ungroup()

# 2. Create the Simplified Visualization
simplified_biomarker_plot <- ggplot(plot_data_simplified, aes(x=reorder(Label, Beta), y=Beta, fill=Phylum)) +
  geom_bar(stat="identity", width = 0.7) +
  coord_flip() +
  # 'free_y' allows each panel to show only its specific top taxa
  facet_wrap(~Stage, scales = "free_y") + 
  theme_bw() +
  scale_fill_viridis_d(option="turbo") +
  labs(title="Key Microbial Signatures of Foot Disease",
       subtitle="Top 15 Significant Biomarkers per stage (|LFC| > 2.5)",
       x="", y="Log-Fold Change (Effect Size)") +
  theme(axis.text.y = element_text(size=9, face="bold.italic"), # Larger, clearer font
        axis.text.x = element_text(size=10),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "bottom")

# 3. Save the new PDF
ggplot2::ggsave("output/differential_abundance/figure_Biomarkers_Simplified.pdf", 
                simplified_biomarker_plot, width=12, height=9)

message("Success! The simplified plot is saved as figure_Biomarkers_Simplified.pdf")
