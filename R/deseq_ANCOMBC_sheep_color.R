##############################################################
# title: "Differential Abundance: DESeq2 and ANCOM-BC"
# author: "Francis Corvin & Gemini AI"
##############################################################


library(tidyverse)
library(qiime2R)
library(phyloseq)
library(DESeq2)
library(ANCOMBC)

# 1. Define Paths
meta_p <- "sheep_project/trial_1/unique_metadata.tsv"
tab_p  <- "sheep_project/trial_1/table_T1.qza"
tax_p  <- "sheep_project/trial_1/taxonomy.qza"

# 2. Load and Parse Taxonomy (CRITICAL FIX)
tax_raw <- read_qza(tax_p)$data
taxonomy_clean <- parse_taxonomy(tax_raw)

# 3. Create Phyloseq Object
ps <- qza_to_phyloseq(features=tab_p, metadata=meta_p)
tax_table(ps) <- as.matrix(taxonomy_clean) # Overwrite with clean taxonomy

# 4. Set Reference Group (Healthy)
sample_data(ps)$disease_state <- factor(
  sample_data(ps)$disease_state, 
  levels = c("A_HEALTHY", "TRANSITION", "C_ID", "D_FOOTROT", "E_CODD", "TREATED")
)

# 5. Method 1: DESeq2 (With the Zero-Fix)
print("Running DESeq2...")
ds <- phyloseq_to_deseq2(ps, ~ disease_state)

# FIX: Use 'poscounts' to handle the zeros in microbiome data
ds <- estimateSizeFactors(ds, type="poscounts") 
ds <- DESeq(ds, test="Wald", fitType="parametric")

# 6. Method 2: ANCOM-BC
print("Running ANCOM-BC...")
out = ancombc(
  data = ps, 
  formula = "disease_state", 
  p_adj_method = "holm", 
  group = "disease_state", 
  struc_zero = TRUE,
  neg_lb = TRUE
)

# 7. Save Results
if(!dir.exists("output/differential_abundance")) dir.create("output/differential_abundance", recursive=TRUE)

write.csv(as.data.frame(results(ds)), "output/differential_abundance/DESeq2_results.csv")
write.csv(as.data.frame(out$res$beta), "output/differential_abundance/ANCOMBC_beta.csv")

message("Success! Both models finished and results are saved.")
