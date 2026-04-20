##############################################################
# title: "Multi-Disease Biomarker Analysis (ID, Footrot, CODD)"
# author: "Francis Corvin & Gemini AI"
##############################################################
library(tidyverse)
library(qiime2R)

# 1. Load Data
# (Assuming taxonomy_clean is still in your environment from the previous steps)
# If not: taxonomy_clean <- parse_taxonomy(read_qza("sheep_project/trial_1/taxonomy.qza")$data)

tax_labels <- taxonomy_clean %>% 
  rownames_to_column("ASV") %>%
  mutate(Label = paste0(Genus, " (", substr(ASV, 1, 4), ")")) %>%
  select(ASV, Label, Phylum)

# 2. Extract ANCOM-BC Beta Coefficients for all major states
# We read the beta CSV we generated earlier
res_abc <- read.csv("output/differential_abundance/ANCOMBC_beta.csv") %>%
  rename(ASV = X)

# 3. Pivot the data to compare ID, Footrot, and CODD side-by-side
comparison_plot_data <- res_abc %>%
  select(ASV, disease_stateC_ID, disease_stateD_FOOTROT, disease_stateE_CODD) %>%
  pivot_longer(cols = starts_with("disease_state"), 
               names_to = "Condition", 
               values_to = "Beta") %>%
  mutate(Condition = case_when(
    Condition == "disease_stateC_ID" ~ "Infectious Dermatitis (ID)",
    Condition == "disease_stateD_FOOTROT" ~ "Footrot",
    Condition == "disease_stateE_CODD" ~ "CODD"
  )) %>%
  left_join(tax_labels, by="ASV")

# 4. Filter for Top Biomarkers
# We identify ASVs that have a strong effect (Beta > 2) in at least one condition
top_biomarkers <- comparison_plot_data %>%
  group_by(ASV) %>%
  filter(any(abs(Beta) > 2.5)) %>% # Adjust this number to show more/fewer taxa
  ungroup()

# 5. Create the Faceted Plot
biomarker_plot <- ggplot(top_biomarkers, aes(x=reorder(Label, Beta), y=Beta, fill=Phylum)) +
  geom_bar(stat="identity") +
  coord_flip() +
  facet_wrap(~Condition, scales = "free_y") + # Shows each disease in its own panel
  theme_bw() +
  scale_fill_viridis_d(option="turbo") +
  labs(title="Key Microbial Biomarkers Across Foot Disease Progression",
       subtitle="Beta coefficients relative to A_HEALTHY baseline",
       x="", y="Log-linear Differential Abundance (Effect Size)") +
  theme(axis.text.y = element_text(face="italic", size=7),
        strip.background = element_rect(fill="grey90"),
        legend.position = "bottom")

# 6. Save
ggsave("output/differential_abundance/figure_MultiDisease_Biomarkers.pdf", biomarker_plot, width=14, height=10)

message("Success! Multi-disease biomarker plot saved.")