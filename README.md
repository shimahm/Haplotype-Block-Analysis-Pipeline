# Haplotype-Block-Analysis-Pipeline

This repository contains scripts for extracting genomic regions around GWAS top SNPs, converting them for haplotype block analysis, and visualizing linkage disequilibrium (LD) structure. Everything is automated for reproducibility and ease of use.

---

## 1. Extract Regional VCFs Around Top SNPs

**Purpose:**  
Given a list of top GWAS SNPs (chromosome & position), extract ±100 kb regions from a compressed VCF file using `bcftools`. This script also handles potential Windows-style file formatting issues in the SNP list.

### a. Preprocess Input Files

```bash
# Index VCF if not already indexed
tabix -p vcf gwas_ready_filtered.vcf.gz

# Check for Windows line endings in SNP list
cat -A top_snps.txt

# If you see "^M" at line ends, convert to Unix format:
dos2unix top_snps.txt
```

### b. Extract Regions Script

```bash
#!/bin/bash

# --- Extraction of ±100 kb regions around each SNP ---

# Input files
VCF=gwas_ready_filtered.vcf.gz
SNP_LIST=top_snps.txt
WINDOW=100000  # 100 kb up/downstream

# Output directory
mkdir -p extracted_regions

# Read SNP list and extract ±100 kb regions
while IFS=$'\t' read -r chrom pos; do
    region_start=$((pos - WINDOW))
    region_end=$((pos + WINDOW))

    # Ensure start is not less than 1
    if [ "$region_start" -lt 1 ]; then
        region_start=1
    fi

    region_name="${chrom}_${pos}"
    echo "Extracting region around ${chrom}:${pos}..."

    # Extract region from VCF and index it
    bcftools view -r "${chrom}:${region_start}-${region_end}" "$VCF" -Oz -o extracted_regions/${region_name}.vcf.gz
    tabix -p vcf extracted_regions/${region_name}.vcf.gz

done < "$SNP_LIST"

echo "✅ All regions extracted to 'extracted_regions/'"
```

---

## 2. Convert to PLINK Format & Run Haplotype Block Analysis

**Purpose:**  
Convert each regional VCF to PLINK binary format and perform haplotype block analysis using PLINK.

```bash
#!/bin/bash

VCF_DIR="extracted_regions"
OUT_DIR="plink_results"
mkdir -p "$OUT_DIR"

for vcf in "$VCF_DIR"/*.vcf.gz; do
    region_name=$(basename "$vcf" .vcf.gz)

    echo "🔄 Converting $region_name to PLINK format..."
    plink --vcf "$vcf" \
          --make-bed \
          --out "$OUT_DIR/$region_name" \
          --allow-extra-chr \
          --double-id  # Required for some datasets

    # Only proceed if PLINK conversion succeeded
    if [[ -f "$OUT_DIR/$region_name.bed" ]]; then
        echo "🧬 Running haplotype block analysis for $region_name..."
        plink --bfile "$OUT_DIR/$region_name" \
              --blocks \
              --out "$OUT_DIR/${region_name}_blocks" \
              --allow-extra-chr
    else
        echo "❌ Failed to generate PLINK files for $region_name. Skipping block analysis."
    fi
done

echo "✅ All haplotype blocks processed (or skipped if conversion failed)."
```

---

## 3. Calculate LD Matrices for Each Region

**Purpose:**  
For each converted region, calculate pairwise LD (r²) using PLINK.

```bash
for region_file in plink_results/*.bed; do
    region_name=$(basename "$region_file" .bed)
    echo "Calculating LD for $region_name"
    plink --bfile plink_results/$region_name \
          --allow-extra-chr \
          --r2 --ld-window 99999 --ld-window-kb 200 --ld-window-r2 0 \
          --out plink_results/${region_name}_ld
done
```

---

## 4. Visualize LD Structure in R

**Purpose:**  
Plot LD heatmaps using [LDheatmap](https://cran.r-project.org/web/packages/LDheatmap/index.html) from PLINK output.

```r
# Set working directory
setwd("Y:/Bigdata/computing/Shima/b_carinata/haplotype_analysis")

# Load required package
if (!require("LDheatmap")) install.packages("LDheatmap"); library(LDheatmap)

# Path to your PLINK LD output file
ld_file <- "plink_results/CM081019.1_13765389_ld.ld"
ld_data <- read.table(ld_file, header = TRUE)

# Fix missing SNP IDs: create IDs from chromosome and position if missing
ld_data$SNP_A <- ifelse(ld_data$SNP_A == ".", paste0(ld_data$CHR_A, "_", ld_data$BP_A), ld_data$SNP_A)
ld_data$SNP_B <- ifelse(ld_data$SNP_B == ".", paste0(ld_data$CHR_B, "_", ld_data$BP_B), ld_data$SNP_B)

# Combine unique SNPs
snps <- unique(c(ld_data$SNP_A, ld_data$SNP_B))

# Create LD matrix
ld_matrix <- matrix(NA, nrow = length(snps), ncol = length(snps))
rownames(ld_matrix) <- snps
colnames(ld_matrix) <- snps

# Fill LD matrix
for (i in 1:nrow(ld_data)) {
  snp_a <- ld_data$SNP_A[i]
  snp_b <- ld_data$SNP_B[i]
  r2 <- ld_data$R2[i]
  
  if (!is.na(snp_a) && !is.na(snp_b) && snp_a %in% snps && snp_b %in% snps) {
    ld_matrix[snp_a, snp_b] <- r2
    ld_matrix[snp_b, snp_a] <- r2
  }
}

# Replace NA with 0 (no LD)
ld_matrix[is.na(ld_matrix)] <- 0

# Build real genetic map based on base pair positions
bp_info <- unique(rbind(
  data.frame(SNP = ld_data$SNP_A, BP = ld_data$BP_A),
  data.frame(SNP = ld_data$SNP_B, BP = ld_data$BP_B)
))
bp_info <- bp_info[match(rownames(ld_matrix), bp_info$SNP), ]  # match order

# Create genetic map
genetic_map <- bp_info$BP
names(genetic_map) <- bp_info$SNP

# Plot LD heatmap
LDheatmap(ld_matrix,
          genetic_map,
          color = heat.colors(20),
          name = "LD",
          title = "LD heatmap (r²) using real coordinates",
          add.key = TRUE)
```

---

## Notes

- **File Names:** All scripts assume input/output file names as shown above. Rename as needed.
- **Dependencies:**  
  - [`bcftools`](http://samtools.github.io/bcftools/bcftools.html)  
  - [`plink`](https://www.cog-genomics.org/plink/)  
  - [`LDheatmap` R package](https://cran.r-project.org/web/packages/LDheatmap/index.html)  
- **Reference Genome:** If allele checks are needed, provide reference genome path accordingly.
- **SNP List Format:** Tab-delimited with columns: `chromosome<TAB>position`

