# Alzheimer's Transcriptomic and Network Analysis

### Transcriptomic and network analysis of hippocampal gene expression in Alzheimer's disease using GSE1297

## Overview

This project analyzes hippocampal gene expression in Alzheimer's disease using the publicly available GSE1297 microarray dataset. The analysis compares 9 control samples with 7 severe Alzheimer's disease samples using the Affymetrix Human Genome U133A platform.

The workflow combines differential gene expression analysis, functional enrichment analysis, and protein-protein interaction (PPI) network analysis to explore molecular patterns associated with severe Alzheimer's disease.

## Objectives

- Identify genes showing differential expression between control and severe Alzheimer's disease hippocampal samples.
- Investigate biological processes and pathways represented among the top-ranked genes.
- Construct a protein-protein interaction network and identify highly connected candidate hub genes.
- Explore publication enrichment associated with the analyzed gene set.

## Dataset

Dataset: GSE1297

Platform: Affymetrix Human Genome U133A (HG-U133A)

Samples analyzed: 16

- Control: 9
- Severe Alzheimer's disease: 7

A total of 22,283 probes were included in the differential expression analysis.

## Analysis Workflow

GSE1297
↓
Data preprocessing
↓
Log2 transformation
↓
Differential expression analysis
↓
Top 150 probes by nominal P-value
↓
Probe annotation
↓
143 unique genes
↓
GO enrichment
↓
KEGG enrichment
↓
STRING PPI network
↓
Hub gene analysis
↓
Publication enrichment

## Tools and Resources

- R 4.6.1 - Statistical analysis
- GEOquery - GEO dataset retrieval
- limma - Differential expression analysis
- hgu133a.db - Probe annotation
- Enrichr - GO and KEGG enrichment
- STRING v12.0 - PPI network and functional enrichment

## Key Results

### Differential Expression

A total of 3,482 probes had nominal P < 0.05. However, no probe reached FDR-adjusted P < 0.05, with the smallest adjusted P-value being 0.117.

Therefore, the top 150 probes were retained as an exploratory ranked gene set, yielding 143 unique genes after annotation and removal of duplicate gene entries.

### GO Enrichment

The 143-gene set was tested against 1,273 GO Biological Process terms.

The lowest adjusted P-value was 0.0947, meaning that no GO term reached statistical significance after multiple-testing correction.

The highest-ranked terms included negative regulation of cell growth, miRNA-related processes, and sex-hormone-related processes.

### KEGG Enrichment

The gene set was tested against 175 KEGG pathways.

The strongest nominal enrichment was observed for transcriptional misregulation in cancer, but no pathway remained significant after correction. The lowest adjusted P-value was 0.42.

### PPI Network

STRING mapped 137 of the 143 genes to nodes, producing a network containing 81 edges.

The most connected genes were:

NRXN1 - degree 8

CSF2 - degree 7

GATA3 - degree 7

GFAP - degree 6

CD68 - degree 6

The PPI enrichment P-value was 0.060, indicating that the network did not contain significantly more interactions than expected by chance.

## Interpretation

The analysis did not identify statistically significant individual genes or GO/KEGG enrichment after multiple-testing correction. Therefore, the downstream findings should be considered exploratory rather than definitive evidence of Alzheimer's disease-associated molecular changes.

STRING functional enrichment identified several neuronal and synaptic categories, but these results used STRING's default whole-genome background and were therefore interpreted cautiously.

## Report

The complete analysis and results are available in the PDF report.

[Download the full analysis report](./Transcriptomic%20and%20Network%20Analysis%20of%20Alzheimer's%20Disease%20Hippocampal%20Gene%20Expression.pdf)

## Limitations

- Small sample size (16 samples).
- No individual probe reached FDR-adjusted significance.
- The 143-gene set was selected based on nominal P-value ranking and is therefore exploratory.
- GO and KEGG analyses did not produce statistically significant results after correction.
- STRING enrichment used a whole-genome background rather than an array-specific background.

## Project Status

Completed - exploratory transcriptomic and network analysis.
