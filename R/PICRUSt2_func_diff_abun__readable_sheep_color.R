##############################################################
# title: "Complete Functional Analysis with Descriptive Names"
# author: "Francis Corvin & Gemini AI"
##############################################################

# 1. Load All Required Libraries
library(tidyverse)
library(ANCOMBC)
library(phyloseq)
library(ggplot2)

# 2. SET WORKING DIRECTORY
setwd("C:/Users/Francis/Documents/ANSC_516_Rstudio")

# 3. Define File Paths
path_file    <- "sheep_project/trial_1/swab_core_metrics/picrust2_output/pathways_out/path_abun_unstrat.tsv.gz"
meta_file    <- "sheep_project/trial_1/unique_metadata.tsv"
mapping_file <- "sheep_project/trial_1/swab_core_metrics/picrust2_output/metacyc_pathways_info.txt.gz"
output_dir   <- "output/functional_abundance/"

if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

# 4. Load Pathway Abundance Data
# R handles .gz automatically; "pathway" is the ID column
path_df <- read_delim(path_file, delim="\t", show_col_types = FALSE) %>%
  column_to_rownames("pathway")

# 5. Load Metadata (Using your manual fix: SampleID)
metadata <- read_delim(meta_file, delim="\t", show_col_types = FALSE) %>%
  filter(SampleID %in% colnames(path_df)) %>%
  filter(disease_state %in% c("A_HEALTHY", "C_ID", "D_FOOTROT", "E_CODD")) %>%
  column_to_rownames("SampleID")

# 6. Load the Mapping File (The one you just downloaded)
# This file typically has two columns: the ID and the English description
path_metadata <- read_delim(mapping_file, delim="\t", col_names = c("taxon", "Description"), show_col_types = FALSE)

# 7. Build Phyloseq and Run ANCOM-BC2
path_df_sub <- path_df[, rownames(metadata)]
ps_path <- phyloseq(otu_table(as.matrix(path_df_sub), taxa_are_rows = TRUE), 
                    sample_data(metadata))

print("Running ANCOM-BC2 math... please wait.")
out_path <- ancombc2(data = ps_path, fix_formula = "disease_state", p_adj_method = "holm", 
                     group = "disease_state", struc_zero = TRUE, neg_lb = TRUE)

# 8. Merge Results with Human-Readable Names
res_path <- as.data.frame(out_path$res) %>%
  # Join with the mapping file to get the 'Description' column
  dplyr::left_join(path_metadata, by = "taxon") %>%
  # Fallback: if a name isn't in the map, use the ID so the code doesn't break
  dplyr::mutate(Description = ifelse(is.na(Description), taxon, Description)) %>%
  tidyr::pivot_longer(cols = contains("lfc_disease_state"), 
                      names_to = "Stage", 
                      values_to = "LFC") %>%
  dplyr::mutate(Stage = dplyr::case_when(
    grepl("C_ID", Stage) ~ "ID", 
    grepl("D_FOOTROT", Stage) ~ "Footrot", 
    TRUE ~ "CODD"
  )) %>%
  # Filter for high impact and top 15 for readability
  dplyr::filter(abs(LFC) > 1.5) %>% 
  dplyr::group_by(Stage) %>%
  dplyr::slice_max(order_by = abs(LFC), n = 15) %>%
  dplyr::ungroup()

# 9. Create the Final Descriptive Plot
functional_plot <- ggplot(res_path, aes(x=reorder(Description, LFC), y=LFC, fill=LFC > 0)) +
  geom_bar(stat="identity", width = 0.7) + 
  coord_flip() + 
  facet_wrap(~Stage, scales="free_y") +
  theme_bw() + 
  scale_fill_manual(values=c("skyblue3", "orangered"), labels=c("Lower in Disease", "Higher in Disease")) +
  labs(title="Key Metabolic Functions by Disease Stage", 
       subtitle="Top 15 Pathways per stage (|LFC| > 1.5) with Full Descriptions",
       x="", y="Log-Fold Change (vs Healthy)",
       fill="Direction") +
  theme(axis.text.y = element_text(size=7, face="bold"), 
        strip.text = element_text(size=11, face="bold"),
        legend.position = "bottom")

# 10. Save Output
ggsave(paste0(output_dir, "figure_Functional_Descriptive.pdf"), functional_plot, width=16, height=10)
write.csv(res_path, paste0(output_dir, "ANCOMBC_pathway_descriptive_results.csv"), row.names=FALSE)

message("Success! Your final descriptive functional plot is ready.")

# 1. Ensure the output directory exists
output_dir <- "output/functional_abundance/individual_plots/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

# 2. Create a list of the stages we want to plot
stages <- c("ID", "Footrot", "CODD")

# 3. Use a loop to create and save each plot separately
for (s in stages) {
  
  # Filter data for just this stage
  stage_data <- res_path %>%
    dplyr::filter(Stage == s)
  
  # Create the plot
  p <- ggplot(stage_data, aes(x=reorder(Description, LFC), y=LFC, fill=LFC > 0)) +
    geom_bar(stat="identity", width = 0.7) + 
    coord_flip() + 
    theme_bw() + 
    scale_fill_manual(values=c("skyblue3", "orangered"), labels=c("Lower", "Higher")) +
    labs(title=paste("Metabolic Shifts in", s), 
         subtitle="Top 15 Functional Biomarkers | ANCOM-BC2 (q < 0.05)",
         x="", y="Log-Fold Change (vs Healthy)",
         fill="Direction") +
    theme(axis.text.y = element_text(size=10, face="bold"), # Larger text now that there is room
          plot.title = element_text(size=16, face="bold"),
          legend.position = "bottom")
  
  # Save the PDF with a specific name
  file_name <- paste0(output_dir, "Functional_Pathways_", s, ".pdf")
  ggsave(file_name, p, width=12, height=8)
  
  message(paste("Saved:", file_name))
}

message("All individual PDFs are ready in output/functional_abundance/individual_plots/")
