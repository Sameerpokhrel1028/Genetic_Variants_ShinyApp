# Genetic Variants Shiny App  
## Interactive Visualization of Genomic Variants and Genetic Relationships

This repository provides a complete workflow to:

- Generate window-based variant density matrices from VCF files  
- Compute IBS-based genetic distances (1 − IBS) across multiple genomic scales  
- Explore genome-wide, chromosome-level, and regional variation  
- Visualize results interactively using an R Shiny application  

The workflow is designed for large-scale genomics datasets and is species-independent.

---

# Quick Start

## 0. Load Required Software

Required tools:

- bcftools  
- bedtools  
- plink (for IBS)  
- tabix (recommended)  

On HPC systems:

    module load BCFtools
    module load BEDTools
    module load PLINK

Index your VCF (recommended):

    tabix -p vcf my.vcf.gz

---

## 1. Clone Repository

    git clone https://github.com/Sameerpokhrel1028/Genetic_Variants_ShinyApp.git
    cd Genetic_Variants_ShinyApp

---

#  Step 1: Variant Density Matrix

Script:

    scripts/parse_vcf_to_matrix.sh

Basic usage:

    bash scripts/parse_vcf_to_matrix.sh \
      -v my.vcf.gz \
      -o my_10kb.tsv

Optional parameters:

    -w 10000    # window size (default: 10 kb)
    -t 16       # threads (optional)

Output format:

    Chrom   Start   End   Sample1   Sample2   ...

What is counted:

    0/1, 1/1, 0|1, 1|0, 1/2, etc.

Ignored:

    0/0
    ./.

---

#  Step 2: IBS Distance (Genome + Chromosome + Region)

Script:

    scripts/compute_ibs_chromosomes.sh

Basic usage:

    bash scripts/compute_ibs_chromosomes.sh \
      --vcf my.vcf.gz \
      --out ibs_scope.tsv

Optional parameters:

    --region-mb 5     # region bin size (default: 5 Mb)
    --max-miss 0.20   # missingness filter
    --memory 16000    # PLINK memory

---

Output format:

    Scope   Sample1   Sample2   Distance

Scope types:

- Genome  
- Chromosome (e.g., TRv2Chr.01)  
- Region bins (e.g., TRv2Chr.01:1-5000000)  

Distance definition:

    Distance = 1 − IBS

---

#  Step 3: Shiny Application

Hosted App:

https://sameerpokhrel.shinyapps.io/rshinyappvariantsobservations/

Upload:

- Variant matrix (my_10kb.tsv)  
- IBS file (ibs_scope.tsv)  

---

# 🔬 Workflow Overview

    Multi-sample VCF
            ↓
    ----------------------------------
    | Variant Density Pipeline       |
    ----------------------------------
    Non-ref genotype extraction
            ↓
    Window-based counting (bedtools)
            ↓
    Variant density matrix (TSV)
            ↓
    Heatmap visualization

    ----------------------------------
    | IBS Distance Pipeline          |
    ----------------------------------
    Biallelic SNP filtering
            ↓
    Missingness filtering
            ↓
    PLINK IBS calculation
            ↓
    1 - IBS conversion
            ↓
    Genome + Chromosome + Region bins
            ↓
    Clustering / PCA / Heatmap

---

#  Variant Density Module

Script:

    scripts/parse_vcf_to_matrix.sh

Features:

- Window-based counting (default 10 kb)  
- Per-sample variant matrix  
- Efficient for large VCFs  
- Species-independent  

---

#  IBS Clustering Module

Script:

    scripts/compute_ibs_chromosomes.sh

Features:

- Genome-wide IBS  
- Chromosome-level IBS  
- Region-based IBS (5 Mb bins default)  
- Fast PLINK-based computation  

---

#  R Shiny Application Features

## Variant Density Heatmap

- Variant density visualization (variants/kb)  
- Log-transformed scale  
- Multi-sample comparison  
- Chromosome selection  
- Publication-quality output  

---

## Genome Clustering / IBS

Neighbor-Joining Tree:
- IBS-based clustering  
- Thick branches and bold labels  

Distance Heatmap:
- Pairwise similarity visualization  
- Optional clustering  

PCA (MDS):
- Dimensionality reduction  
- Population structure visualization  

---

#  Hierarchical Scope Selection

| Level | Description |
|------|------------|
| All | Whole genome |
| Chromosome | Single chromosome |
| Region bins | Sub-chromosomal windows |

Region workflow:

    Select "Region bins"
            ↓
    Select chromosome
            ↓
    Select region (0–5 Mb, 5–10 Mb, etc.)

---

#  Design Philosophy

- Heavy computation done outside Shiny  
- Shiny used only for visualization  

This ensures:

- Scalability  
- Speed  
- Flexibility  

---

#  Generalizability


Chromosome naming is handled dynamically. Tested with a few VCF files in peanut. Theoretically, it should work with other species.

---

#  Requirements

Variant Density:

    bcftools
    bedtools
    tabix (recommended)

IBS Distance:

    bcftools
    plink

---

# Summary

This workflow provides:

- Efficient VCF → matrix conversion  
- Multi-scale IBS distance computation  
- Interactive visualization  
- Publication-ready outputs  
