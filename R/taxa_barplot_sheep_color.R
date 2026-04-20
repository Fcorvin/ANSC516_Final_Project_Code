##############################################################
# title: "Comprehensive Taxonomy Barplots"
# author: "Francis Corvin & Gemini AI"
##############################################################
library(tidyverse)
library(qiime2R)

# Load data
meta <- read_q2metadata("sheep_project/trial_1/unique_metadata.tsv")
taxonomy <- read_qza("sheep_project/trial_1/taxonomy.qza")$data
table <- read_qza("sheep_project/trial_1/table_T1.qza")$data

# Define levels to plot
levels <- c("Phylum", "Family", "Genus", "ASV")

for (lvl in levels) {
  # Summarize taxa at the current level
  if (lvl == "ASV") {
    tax_table <- table %>% as.data.frame() %>% rownames_to_column("ASV")
  } else {
    tax_table <- summarize_taxa(table, taxonomy)[[lvl]] %>% rownames_to_column(lvl)
  }
  
  # Prepare for plotting: Top 15 taxa, others as "Other"
  plot_data <- tax_table %>%
    pivot_longer(-1, names_to="SampleID", values_to="Abundance") %>%
    inner_join(meta, by="SampleID") %>%
    group_by(disease_state, !!sym(lvl)) %>%
    summarize(MeanAbund = mean(Abundance), .groups="drop")
  
  # Identify top 15 taxa globally to keep plot clean
  top_taxa <- plot_data %>%
    group_by(!!sym(lvl)) %>%
    summarize(Total = sum(MeanAbund)) %>%
    top_n(15, Total) %>%
    pull(1)
  
  plot_data <- plot_data %>%
    mutate(Taxa_Fixed = ifelse(!!sym(lvl) %in% top_taxa, as.character(!!sym(lvl)), "Other"))
  
  # Plot
  p <- ggplot(plot_data, aes(x=disease_state, y=MeanAbund, fill=Taxa_Fixed)) +
    geom_bar(stat="identity", position="fill") +
    scale_fill_viridis_d(option="turbo") +
    theme_q2r() +
    theme(axis.text.x = element_text(angle=45, hjust=1), legend.text = element_text(size=8)) +
    labs(y="Relative Abundance", x="", title=paste(lvl, "Level Composition"), fill=lvl)
  
  ggsave(paste0("output/figure_taxa_", lvl, ".pdf"), p, height=7, width=10)
}