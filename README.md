# Peanut-PanMAGIC Variant Density Platform

This repository provides an end-to-end workflow for exploring genome-wide
variant density and divergence patterns across peanut genomes derived from
a pangenome analysis.

The platform consists of three linked components:
1) Variant extraction from a pangenome VCF
2) Window-based variant summarization and filtering
3) Interactive visualization using a Shiny application

------------------------------------------------------------
Workflow Overview
------------------------------------------------------------

Pangenome VCF
  ↓
Variant counting (10 kb / 100 kb windows)
  ↓
Filtering (>40 variants per kb)
  ↓
Sample-wise TSV matrix
  ↓
Shiny app (visualization & exploration)

Each step is modular and can be run independently.

------------------------------------------------------------
1. Variant Counting from Pangenome VCF
------------------------------------------------------------

Input:
- Multi-sample pangenome VCF
  e.g. peanutpan.vcf.phased.gz

Script:
- scripts/parse_vcf_to_10kb_matrix.sh

What this script does:
- Extracts non-reference variants for each genome
- Bins variants into fixed genomic windows (default: 10 kb)
- Generates a wide matrix of variant counts per window per sample

Output:
- samplewise_variant_matrix_10kb.tsv

This step is intended to be run on an HPC system or workstation
with bcftools and bedtools installed.

------------------------------------------------------------
2. Filtering and Normalization
------------------------------------------------------------

Hypervariable windows (>40 variants per kb), which often reflect
repetitive regions or mapping artifacts, are filtered downstream.

Filtering and normalization (variants per kb) are performed within
the Shiny application to keep preprocessing lightweight and flexible.

------------------------------------------------------------
3. Interactive Visualization (Shiny App)
------------------------------------------------------------

App files:
- shiny_app/app.R
- shiny_app/functions.R

Input file:
- TSV matrix generated from the variant parsing step
  e.g. samplewise_variant_matrix_10kb.tsv

What users can do in the app:
- Select chromosome (chr01–chr20)
- Select one or multiple genomes
- Visualize genome-wide variant density
- Explore genome divergence using dendrograms and PCA

------------------------------------------------------------
Live Application
------------------------------------------------------------

The interactive Shiny application is publicly available at:

https://sameerpokhrel.shinyapps.io/shinyapp1kbvariants/

No local installation is required to explore the results.

------------------------------------------------------------
About This App
------------------------------------------------------------

This Shiny application visualizes variant density and clusters genotypes
for a peanut 16-way MAGIC population developed in the Ozias-Akins lab
at the University of Georgia.

Tabs in the app:

- Variants Heatmap
  Shows variant-density heatmaps for chromosomes 1–20 across
  selected peanut genomes.

- Genome Clusters
  Displays genome divergence using hierarchical clustering
  (dendrogram) and principal component analysis (PCA).

- About
  Describes the purpose of the app, input data, and outputs.

Input:
- Variant counts per 10 kb windows for 18 peanut genomes derived
  from a pangenome VCF.

Output:
- Variant-density heatmaps
- PCA plots
- Dendrograms for identifying divergence patterns and potential
  introgressed genomic regions.

------------------------------------------------------------
Notes
------------------------------------------------------------

- This repository contains scripts and documentation only
- Raw VCFs and large intermediate files are not included
- Window size (10 kb, 100 kb, etc.) can be modified in the parsing script
- Variant parsing and visualization are intentionally separated
  for performance, clarity, and reuse
