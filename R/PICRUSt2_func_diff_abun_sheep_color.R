##############################################################
# title: "Functional Differential Abundance (PICRUSt2)"
# author: "Francis Corvin & Gemini AI"
##############################################################


library(tidyverse)
library(ANCOMBC)
library(phyloseq)
library(ggplot2)

# 1. Set Directory
setwd("C:/Users/Francis/Documents/ANSC_516_Rstudio")

# 2. Load Pathway Data (.gz)
path_file <- "sheep_project/trial_1/swab_core_metrics/picrust2_output/pathways_out/path_abun_unstrat.tsv.gz"
path_df <- read_delim(path_file, delim="\t", show_col_types = FALSE) %>%
  column_to_rownames("pathway")

# 3. Load Metadata (Using your newly cleaned SampleID column)
metadata <- read_delim("sheep_project/trial_1/unique_metadata.tsv", delim="\t", show_col_types = FALSE) %>%
  filter(SampleID %in% colnames(path_df)) %>%
  filter(disease_state %in% c("A_HEALTHY", "C_ID", "D_FOOTROT", "E_CODD")) %>%
  column_to_rownames("SampleID")

# 4. Build Phyloseq and Run ANCOM-BC2
path_df_sub <- path_df[, rownames(metadata)]
ps_path <- phyloseq(otu_table(as.matrix(path_df_sub), taxa_are_rows = TRUE), sample_data(metadata))

print("Running analysis...")
out_path <- ancombc2(data = ps_path, fix_formula = "disease_state", p_adj_method = "holm", 
                     group = "disease_state", struc_zero = TRUE, neg_lb = TRUE)

# 5. Extract and Plot Top 15 Pathways per Stage
res_path <- as.data.frame(out_path$res) %>%
  rename(Pathway = taxon) %>%
  pivot_longer(cols = contains("lfc_disease_state"), names_to = "Stage", values_to = "LFC") %>%
  mutate(Stage = case_when(grepl("C_ID", Stage) ~ "ID", grepl("D_FOOTROT", Stage) ~ "Footrot", TRUE ~ "CODD"),
         Pathway_Clean = gsub("_", " ", Pathway)) %>%
  filter(abs(LFC) > 1.5) %>% 
  group_by(Stage) %>%
  slice_max(order_by = abs(LFC), n = 15) %>%
  ungroup()

functional_plot <- ggplot(res_path, aes(x=reorder(Pathway_Clean, LFC), y=LFC, fill=LFC > 0)) +
  geom_bar(stat="identity") + coord_flip() + facet_wrap(~Stage, scales="free_y") +
  theme_bw() + scale_fill_manual(values=c("skyblue3", "orangered"), labels=c("Lower", "Higher")) +
  labs(title="Metabolic Shifts", x="", y="Log-Fold Change") +
  theme(axis.text.y = element_text(size=7), legend.position = "bottom")

ggsave("output/functional_abundance/figure_Functional_Pathways.pdf", functional_plot, width=14, height=9)

# 5. Extract and Plot Top 15 Pathways per Stage
# We use the explicit 'dplyr::' prefix to make sure R finds the column names
res_path <- as.data.frame(out_path$res) %>%
  dplyr::rename(Pathway = taxon) %>%
  tidyr::pivot_longer(cols = contains("lfc_disease_state"), 
                      names_to = "Stage", 
                      values_to = "LFC") %>%
  dplyr::mutate(Stage = dplyr::case_when(
    grepl("C_ID", Stage) ~ "ID", 
    grepl("D_FOOTROT", Stage) ~ "Footrot", 
    TRUE ~ "CODD"
  ),
  Pathway_Clean = gsub("_", " ", Pathway)) %>%
  dplyr::filter(abs(LFC) > 1.5) %>% 
  dplyr::group_by(Stage) %>%
  dplyr::slice_max(order_by = abs(LFC), n = 15) %>%
  dplyr::ungroup()

# Generate the plot
functional_plot <- ggplot(res_path, aes(x=reorder(Pathway_Clean, LFC), y=LFC, fill=LFC > 0)) +
  geom_bar(stat="identity") + 
  coord_flip() + 
  facet_wrap(~Stage, scales="free_y") +
  theme_bw() + 
  scale_fill_manual(values=c("skyblue3", "orangered"), labels=c("Lower", "Higher")) +
  labs(title="Predicted Metabolic Shifts", 
       subtitle="Top 15 Functional Biomarkers per stage | ANCOM-BC2",
       x="", y="Log-Fold Change (vs Healthy)") +
  theme(axis.text.y = element_text(size=7, face="bold"), 
        legend.position = "bottom")

# Save the PDF
if(!dir.exists("output/functional_abundance/")) dir.create("output/functional_abundance/", recursive=TRUE)
ggsave("output/functional_abundance/figure_Functional_Pathways.pdf", functional_plot, width=14, height=9)

message("Success! The functional plot is now in your output folder.")

# 1. Create a quick lookup for common MetaCyc prefixes to make them slightly more readable
# or search for a standard mapping. 
# For now, we will clean the IDs and ensure the plot is legible.

res_path_mapped <- as.data.frame(out_path$res) %>%
  dplyr::rename(Pathway = taxon) %>%
  tidyr::pivot_longer(cols = contains("lfc_disease_state"), 
                      names_to = "Stage", 
                      values_to = "LFC") %>%
  dplyr::mutate(
    Stage = dplyr::case_when(
      grepl("C_ID", Stage) ~ "ID", 
      grepl("D_FOOTROT", Stage) ~ "Footrot", 
      TRUE ~ "CODD"
    ),
    # This removes common prefixes like 'PWY-' or 'METACYC-' to make them look cleaner
    # while keeping the ID so you can still search it.
    Pathway_Label = gsub("PWY-|METACYC-", "", Pathway),
    Pathway_Label = gsub("_", " ", Pathway_Label)
  ) %>%
  dplyr::filter(abs(LFC) > 1.5) %>% 
  dplyr::group_by(Stage) %>%
  dplyr::slice_max(order_by = abs(LFC), n = 15) %>%
  dplyr::ungroup()

# 2. Update the Plot with the new labels
functional_plot_readable <- ggplot(res_path_mapped, aes(x=reorder(Pathway_Label, LFC), y=LFC, fill=LFC > 0)) +
  geom_bar(stat="identity", width = 0.7) + 
  coord_flip() + 
  facet_wrap(~Stage, scales="free_y") +
  theme_bw() + 
  scale_fill_manual(values=c("skyblue3", "orangered"), labels=c("Lower in Disease", "Higher in Disease")) +
  labs(title="Key Metabolic Functions by Disease Stage", 
       subtitle="Top 15 Pathways (|LFC| > 1.5). Names derived from MetaCyc IDs.",
       x="Metabolic Pathway", y="Log-Fold Change (vs Healthy)") +
  theme(axis.text.y = element_text(size=8, face="bold"), 
        strip.text = element_text(size=12, face="bold"),
        legend.position = "bottom")

# 3. Save the new version
ggsave("output/functional_abundance/figure_Functional_Readable.pdf", 
       functional_plot_readable, width=14, height=9)

message("Success! A more readable PDF is ready.")
