##############################################################
# title: "PICRUSt2 Functional Analysis & DESeq2"
# author: "Francis Corvin & Gemini AI"
##############################################################
library(tidyverse)
library(qiime2R)
library(DESeq2)

# 1. Define the correct path based on your folder structure
pathway_path <- "sheep_project/trial_1/swab_core_metrics/picrust2_output/pathways_out/path_abun_unstrat.tsv.gz"
metadata_path <- "sheep_project/trial_1/unique_metadata.tsv"

# 2. Load Pathway Abundance (R handles .gz automatically)
if(!file.exists(pathway_path)) stop("Pathway file not found! Check the path.")
pathway_table <- read.delim(pathway_path, row.names=1, check.names=F)

# 3. Load Metadata
meta <- read_q2metadata(metadata_path)

# 4. MATCH SAMPLES
# Ensure we only use samples present in both the pathway table and metadata
common_samples <- intersect(colnames(pathway_table), meta$SampleID)
pathway_filtered <- pathway_table[, common_samples]
meta_filtered <- meta %>% filter(SampleID %in% common_samples)

# Reorder metadata to match the columns of the pathway table exactly
meta_filtered <- meta_filtered[match(colnames(pathway_filtered), meta_filtered$SampleID),]

# 5. Prepare DESeq2 Object
# We round the PICRUSt2 counts because DESeq2 requires integers
ds_counts <- round(as.matrix(pathway_filtered))

dds <- DESeqDataSetFromMatrix(
  countData = ds_counts,
  colData = meta_filtered,
  design = ~ disease_state
)

# 6. Run Differential Abundance
dds <- DESeq(dds)

# 7. Extract Results (Comparing CODD vs Healthy)
res <- results(dds, contrast=c("disease_state", "E_CODD", "A_HEALTHY"))
res_df <- as.data.frame(res) %>%
  rownames_to_column("Pathway") %>%
  filter(padj < 0.05) %>%
  arrange(log2FoldChange)

# 8. Professional Visualization of Top Pathways
if(!dir.exists("output/picrust")) dir.create("output/picrust", recursive = TRUE)

picrust_plot <- res_df %>%
  top_n(30, abs(log2FoldChange)) %>%
  ggplot(aes(x=reorder(Pathway, log2FoldChange), y=log2FoldChange, fill=log2FoldChange > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values=c("#4169E1", "#CD5C5C"), labels=c("Down-regulated", "Up-regulated")) +
  theme_bw() +
  labs(title="Metabolic Pathways: CODD vs Healthy",
       subtitle="Significant pathways (p-adj < 0.05) via DESeq2",
       x="MetaCyc Pathway ID", y="Log2 Fold Change", fill="Direction")

ggsave("output/picrust/figure_PICRUSt2_DESeq2.pdf", picrust_plot, height=10, width=12)

write.csv(res_df, "output/picrust/DESeq2_full_results.csv")

message("Success! PICRUSt2 results saved to output/picrust/")
