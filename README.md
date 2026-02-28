# Interactive Visualization of Variants from VCF Files

This workflow enables genome-wide exploration of:

- Window-based variant density
- IBS-derived genetic distance (1-IBS)
- Genome-wide and chromosome-level divergence patterns
- Interactive clustering and PCA visualization

It converts a multi-sample VCF into structured matrices for scalable analysis and visualization.


## Quick Start

### 0. Load Required Modules / Software

The workflow requires:

- bcftools
- bedtools
- plink (for IBS distance)
- python3
- tabix (recommended for indexed VCFs)

On HPC systems:

```bash
module load BCFtools
module load BEDTools
module load PLINK
```

Make sure your VCF is indexed:

```bash
tabix -p vcf my.vcf.gz
```

---

### 1. Clone the repository

```bash
git clone https://github.com/Sameerpokhrel1028/Genetic_Variants_ShinyApp.git
cd Genetic_Variants_ShinyApp
```

### 2. Generate a Window-Based Variant Matrix (Variant Density)

Counts non-reference genotypes per genomic window (default: 10 kb).

```bash
bash scripts/parse_vcf_to_matrix.sh -v my.vcf.gz -o my_10kb.tsv
```

Optional: change resolution and threads

```bash
bash scripts/parse_vcf_to_matrix.sh -v my.vcf.gz -w 100000 -t 16 -o my_100kb.tsv
```

### 3. Compute IBS-Derived Genetic Distance (1-IBS)

Computes genome-wide and per-chromosome genetic distance using PLINK.

```bash
bash scripts/compute_ibs_chromosomes.sh \
  --vcf my.vcf.gz \
  --out ibs_scope.tsv
```

Optional parameters:

- `--max-miss 0.20`  Maximum allowed missingness (default: 0.20)
- `--memory 16000`   PLINK memory in MB

Output format:

```
Scope   Sample1   Sample2   Distance
```

Where:

- `Scope` = Genome or chromosome ID (chr01, chr02, etc.)
- `Distance` = 1 - IBS


### 4. Open the Hosted Shiny Application

https://sameerpokhrel.shinyapps.io/rshinyappvariantsobservations/

Upload:

- The window-based variant matrix (e.g., `my_10kb.tsv`)
- The IBS distance file (`ibs_scope.tsv`)


## Workflow Overview

```
Multi-sample VCF
        ↓
(1) Window-based non-reference genotype counts
        ↓
Variant Density Matrix (TSV)
        ↓
Heatmaps and clustering in Shiny

AND

Multi-sample VCF
        ↓
Biallelic SNP filtering + missingness filter
        ↓
PLINK IBS computation
        ↓
1-IBS distance (Genome + Chromosome)
        ↓
Hierarchical clustering and PCA
```


## Variant Density Module

Script:

```
scripts/parse_vcf_to_matrix.sh
```

For each sample and each genomic window, counts genotypes containing any ALT allele:

```
0/1, 1/1, 0|1, 1|0, 1/2, 2/2, etc.
```

Ignored genotypes:

```
0/0
./.
```

Output format:

```
Chrom   Start   End   Sample1   Sample2   ...
```


## Identity Distance (1-IBS) Module

Script:

```
scripts/compute_ibs_chromosomes.sh
```

Processing steps:

1. Filter to biallelic SNPs
2. Remove SNPs exceeding missingness threshold
3. Compute IBS using PLINK
4. Convert similarity to genetic distance (1 - IBS)

Output (long format):

```
Scope   Sample1   Sample2   Distance
```

Includes:

- Genome-wide distance
- Per-chromosome distance


## Requirements

### Variant Density

- bcftools
- bedtools

On HPC systems:

```bash
module load BCFtools
module load BEDTools
```

Recommended:

```bash
tabix -p vcf my.vcf.gz
```


### IBS Distance

- bcftools
- plink (1.9 recommended)
- python3
- gunzip

On HPC systems:

```bash
module load BCFtools
module load PLINK
```

