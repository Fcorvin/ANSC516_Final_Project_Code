##############################################################
# title: "ANCOM-BC Volcano Plot: Split Labels for Gained/Lost"
# author: "Francis Corvin & Gemini AI"
##############################################################

library(tidyverse)
library(ggrepel)

# 1. IDENTIFY THE SOURCE DATA
# Based on your ls(), your ANCOM-BC result is likely named 'out'
# We will use 'ps' to get the Taxonomy names
ancom_source <- out 

# 2. EXTRACT & MERGE TAXONOMY
# We take the statistics (LFC, Q-values) and join them with the Genus names
# from your phyloseq object 'ps'
res <- ancom_source$res
tax_df <- as.data.frame(tax_table(ps)) %>% rownames_to_column("TaxaID")

ancombc_labeled <- data.frame(
  TaxaID = row.names(res$lfc),
  lfc = res$lfc[,2],        # Log Fold Change (Column 2 is usually the comparison)
  q_value = res$q_val[,2],  # Adjusted p-values
  significant = res$diff_abn[,2]
) %>%
  left_join(tax_df, by = "TaxaID") %>%
  # If Genus is NA, use the ID; otherwise, use Genus
  mutate(LabelName = ifelse(is.na(Genus), TaxaID, as.character(Genus)))

# 3. SPLIT LABELING LOGIC
# Most significant 'Gained' (Top Right)
top_right <- ancombc_labeled %>%
  filter(lfc > 0, significant == TRUE) %>%
  slice_min(q_value, n = 5)

# Most significant 'Lost' (Top Left)
top_left <- ancombc_labeled %>%
  filter(lfc < 0, significant == TRUE) %>%
  slice_min(q_value, n = 5)

to_label <- bind_rows(top_right, top_left)

# 4. CREATE THE VOLCANO PLOT
volcano_split <- ggplot(ancombc_labeled, aes(x = lfc, y = -log10(q_value))) +
  # Background points: Gray for non-sig, Red for sig
  geom_point(aes(color = significant), alpha = 0.4, size = 2) +
  
  # Add the labels for the Top 5 Gained and Top 5 Lost
  geom_label_repel(data = to_label,
                   aes(label = LabelName),
                   box.padding = 0.6, 
                   point.padding = 0.5,
                   max.overlaps = Inf,
                   segment.color = 'grey30',
                   size = 3.5,
                   fontface = "italic",
                   fill = "white") +
  
  # Visual formatting
  scale_color_manual(values = c("TRUE" = "#CD5C5C", "FALSE" = "grey80"),
                     name = "Significant Shift",
                     labels = c("TRUE" = "q < 0.05", "FALSE" = "n.s.")) +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "royalblue") +
  
  labs(title = "Differential Abundance: Taxa Gained vs. Lost",
       subtitle = "Labels identify top 5 most significant taxa in each direction",
       x = "Log Fold Change (LFC)",
       y = expression(-log[10]~(q-value))) +
  
  theme(legend.position = "bottom",
        plot.title = element_text(face="bold", size = 14),
        axis.title = element_text(face="bold"))

# 5. SAVE
if(!dir.exists("output/ancombc")) dir.create("output/ancombc", recursive = TRUE)
ggsave("output/ancombc/Volcano_Split_Labels.pdf", volcano_split, width=10, height=8)

# Show the result
print(volcano_split)