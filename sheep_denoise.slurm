#!/bin/bash

#SBATCH --nodes=1                      ## number of requested nodes
#SBATCH -n 12                          ## number of requested CPUs
#SBATCH -p cpu                         ## partition
#SBATCH --time=20:00:00                ## 20 hours for 30M reads
#SBATCH --job-name sheep_dada2         ## name of job
#SBATCH -A microbiome                  ## account name

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

# Step 1: Move to your project directory
cd /scratch/bell/fcorvin/sheep_project/

# Step 2: Load the environment (Ensure this matches your version)
module load bioinformatics
module load qiime2

# Step 3: Create folders to keep the results organized
mkdir -p trial_1 trial_2

# --- TRIAL 1: Primers trimmed, No truncation ---
echo "Running Trial 1..."
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 24 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-n-threads 12 \
  --o-table trial_1/table_T1.qza \
  --o-representative-sequences trial_1/rep-seqs_T1.qza \
  --o-denoising-stats trial_1/stats_T1.qza

# --- TRIAL 2: Primers trimmed, Slight Reverse Truncation (200bp) ---
# (This helps merging if the very end of the reverse reads have errors)
echo "Running Trial 2..."
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 24 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 200 \
  --p-n-threads 12 \
  --o-table trial_2/table_T2.qza \
  --o-representative-sequences trial_2/rep-seqs_T2.qza \
  --o-denoising-stats trial_2/stats_T2.qza

# Step 4: Generate visualizations for tomorrow's analysis
echo "Generating visualizations..."
qiime metadata tabulate --m-input-file trial_1/stats_T1.qza --o-visualization trial_1/stats_T1.qzv
qiime metadata tabulate --m-input-file trial_2/stats_T2.qza --o-visualization trial_2/stats_T2.qzv

echo "End time"
date +"%d %B %Y %H:%M:%S"