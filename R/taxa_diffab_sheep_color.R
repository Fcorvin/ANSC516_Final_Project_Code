##############################################################
# title: "Multi-Level Taxonomy Barplots - Sheep CODD Project"
# author: "ANSC595"
##############################################################
library(tidyverse)
library(qiime2R)

# 1. Load and Parse Data
meta <- read_q2metadata("sheep_project/trial_1/unique_metadata.tsv")
table <- read_qza("sheep_project/trial_1/table_T1.qza")$data
tax_raw <- read_qza("sheep_project/trial_1/taxonomy.qza")$data

# This fix is critical for reaching Family and Genus levels
taxonomy <- parse_taxonomy(tax_raw) 

# Set disease state order
meta$disease_state <- factor(meta$disease_state, 
                             levels = c("A_HEALTHY", "TRANSITION", "C_ID", "D_FOOTROT", "E_CODD", "TREATED"))

if(!dir.exists("output")) dir.create("output")

# 2. Loop through each taxonomic level
levels_to_plot <- c("Phylum", "Family", "Genus", "ASV")

for (lvl in levels_to_plot) {
  
  # Step A: Summarize data at the specific level
  if (lvl == "ASV") {
    # For ASVs, we use the raw table data
    st <- table %>% as.data.frame() %>% rownames_to_column("ASV")
  } else {
    # For others, we use the summarized taxonomy levels
    st <- summarize_taxa(table, taxonomy)[[lvl]] %>% rownames_to_column(lvl)
  }
  
  # Step B: Format for plotting and calculate means
  plot_data <- st %>%
    pivot_longer(-1, names_to="SampleID", values_to="Abundance") %>%
    inner_join(meta, by="SampleID") %>%
    group_by(disease_state, !!sym(lvl)) %>%
    summarize(MeanAbundance = mean(Abundance), .groups="drop")
  
  # Step C: Group rare taxa into "Other" (keeps the plot clean/professional)
  # We keep the top 12 most abundant taxa at each level
  top_taxa <- plot_data %>%
    group_by(!!sym(lvl)) %>%
    summarize(Total = sum(MeanAbundance)) %>%
    top_n(12, Total) %>%
    pull(1)
  
  plot_data <- plot_data %>%
    mutate(DisplayTaxa = ifelse(!!sym(lvl) %in% top_taxa, as.character(!!sym(lvl)), "Other"))
  
  # Step D: Create the plot
  p <- ggplot(plot_data, aes(x=disease_state, y=MeanAbundance, fill=DisplayTaxa)) +
    geom_bar(stat="identity", position="fill") +
    scale_fill_viridis_d(option="turbo") + 
    theme_q2r() +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          legend.text = element_text(size=7)) +
    labs(y="Relative Abundance", x="", 
         title=paste(lvl, "Level Composition"),
         fill=lvl)
  
  # Step E: Save
  ggsave(paste0("output/figure_taxa_", lvl, ".pdf"), p, height=7, width=10)
}

message("Success! Phylum, Family, Genus, and ASV plots are now in the output folder.")

