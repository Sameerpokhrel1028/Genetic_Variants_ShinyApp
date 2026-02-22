# Interactive Visualization of Variants from VCF files

## Quick Start

This workflow allows you to explore **genome-wide variant density and divergence patterns** from a multi-sample VCF file.

### 1. Clone the repository

```bash
git clone https://github.com/Sameerpokhrel1028/Genetic_Variants_ShinyApp.git
cd Genetic_Variants_ShinyApp
```

### 2. Generate a window-based variant matrix (TSV)

```bash
bash scripts/parse_vcf_to_matrix.sh -v my.vcf.gz -o my_10kb.tsv
```

### 3. Open the hosted Shiny application

(https://sameerpokhrel.shinyapps.io/rshinyappvariantsobservations/)

Upload `my_10kb.tsv` and begin exploring genome-wide patterns.

---

## Workflow Overview

```
Pangenome VCF
    ↓
Variant counting in fixed genomic windows (e.g. 10 kb or 100 kb)
    ↓
Sample-wise variant density matrix (TSV)
    ↓
Interactive Shiny app for visualization and exploration
```

This workflow enables structured, scalable exploration of genome-wide variant density and divergence patterns across multiple samples or assemblies.

---

# Step 1: Generate the Window-Based Variant Matrix

The parsing script:

```
scripts/parse_vcf_to_matrix.sh
```

Converts a multi-sample VCF into a fixed-window variant count matrix.

For each sample, the script counts **non-reference genotypes per window**.

A genotype is counted if it contains any ALT allele:

```
0/1, 1/1, 0|1, 1|0, 1/2, 2/2, etc.
```

Genotypes `0/0` and `./.` are ignored.

### Output format

```
Chrom   Start   End   Sample1   Sample2   Sample3 ...
```

This TSV file is the direct input for the Shiny application.

---

## Requirements

The parsing step requires:

- bcftools  
- bedtools  

On HPC systems, you may need:

```bash
module load BCFtools
module load BEDTools
```

For best performance, index compressed VCFs:

```bash
tabix -p vcf my.vcf.gz
```

---

## Usage

### Basic command

```bash
bash scripts/parse_vcf_to_matrix.sh -v input.vcf.gz -o output.tsv
```

### Required arguments

- `-v` Input multi-sample VCF (bgzipped `.vcf.gz` recommended)  
- `-o` Output TSV file  

### Optional arguments

- `-w` Window size in base pairs (default: 10000)  
- `-t` Threads for bcftools (default: 4)  
- `-h` Show help  

---

## Example Commands

Create a 10 kb matrix:

```bash
bash scripts/parse_vcf_to_matrix.sh -v my.vcf.gz -o my_10kb.tsv
```

Create a 100 kb matrix using 16 threads:

```bash
bash scripts/parse_vcf_to_matrix.sh -v my.vcf.gz -w 100000 -t 16 -o my_100kb.tsv
```

Window size controls resolution:

- Smaller windows provide fine-scale detail  
- Larger windows smooth variation and reduce runtime  

---

# Step 2: Visualize in the Shiny App (Recommended)

Most users should use the hosted Shiny application.

Open:

https://sameerpokhrel.shinyapps.io/rshinyappvariantsobservations/

Upload the generated TSV file.

The Shiny app provides:

- Genome-wide and chromosome-level variant density heatmaps  
- Dynamic chromosome and sample selection  
- Hierarchical clustering based on window profiles  
- PCA visualization with percent variance explained  

No local R installation is required when using the hosted version.

---

# Running the Shiny App Locally (Optional)

Advanced users may run the app locally.

Install required R packages:

```r
install.packages(c("shiny", "dplyr", "readr", "ggplot2"))
install.packages("ComplexHeatmap")
install.packages("circlize")
```

Launch the app:

```r
library(shiny)
runApp("app.R")
```

Running locally allows modification of visualization settings or integration into larger workflows.

---

# Repository Structure

```
Genetic_Variants_ShinyApp/
├── app.R
├── functions.R
├── README.md
├── scripts/
│   └── parse_vcf_to_matrix.sh
└── 10kbvariants.tsv
```
