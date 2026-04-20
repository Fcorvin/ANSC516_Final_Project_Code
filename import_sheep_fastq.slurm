#!/bin/bash

#SBATCH --nodes=1
#SBATCH -n 20
#SBATCH --time=08:00:00
#SBATCH --job-name=sra_download_fcorvin
#SBATCH -A microbiome
#SBATCH -p cpu

echo "Start time"
date +"%d %B %Y %H:%M:%S"

# 1. Load the SRA Toolkit module
module load biocontainers
module load sra-tools
module list

# 2. Go to your scratch directory
cd /scratch/bell/fcorvin/sheep_project/data_download/

# 3. Create folders
# mkdir -p raw
# mkdir -p seq

# 4. Use prefetch to get the SRA files
# The --resume yes flag allows it to skip the 45 you already downloaded!

prefetch -O raw --resume yes --option-file SRR_Acc_List.txt

cp raw/SRR*/*.sra seq/
cd seq
ls
for i in *.sra; do fasterq-dump ${i} --split-3 -0 .; done
ls

echo "End time"
date +"%d %B %Y %H:%M:%S"
