---
type: analysis-report
title: RNA-seq Differential Expression Analysis
date: 2025-12-28
author: Research Team
version: 1.0
status: Draft
assessment: pass-fail
---

# RNA-seq Differential Expression Analysis

## Objective

Identify differentially expressed genes between treatment and control groups in the RNA-seq dataset to understand the molecular mechanisms of drug response.

### Key Questions

1. Which genes show significant differential expression (FDR < 0.05)?
2. What biological pathways are enriched in the differentially expressed genes?
3. Are the expression changes consistent with the expected mechanism of action?

## Methods

### Data Sources

Raw sequencing data from 24 samples (12 treatment, 12 control) were obtained from the core facility. Quality control was performed using FastQC and MultiQC.

### Analysis Pipeline

```bash
# Alignment with STAR
STAR --genomeDir /ref/hg38 --readFilesIn sample_R1.fq sample_R2.fq

# Quantification with featureCounts
featureCounts -a genes.gtf -o counts.txt aligned.bam

# Differential expression with DESeq2
Rscript run_deseq2.R counts.txt metadata.csv
```

## Results

### Quality Control

| Check | Criterion | Result | Notes |
|-------|-----------|--------|-------|
| Sample quality | RIN > 7 | ✓ | Mean RIN: 8.4 |
| Sequence depth | >20M reads | ✓ | Avg: 32M reads |
| Mapping rate | >80% | ✓ | Mean: 92.1% |
| Duplicate rate | <30% | ✗ | Mean: 34.2% |
| rRNA contamination | <5% | ✓ | Mean: 2.1% |

### Differential Expression

Analysis revealed 540 differentially expressed genes (342 upregulated, 198 downregulated) at FDR < 0.05 threshold.

## Summary

### Key Findings

#### Finding 1: Significant upregulation in oxidative phosphorylation genes

Analysis revealed 342 significantly upregulated genes (FDR < 0.05, log2FC > 1) in the treatment group. The top enriched pathway was oxidative phosphorylation (KEGG, p = 3.2e-12).

#### Finding 2: Downregulation of inflammatory markers

A total of 198 genes showed significant downregulation, with enrichment in inflammatory response pathways (GO:0006954, p = 1.8e-9).

#### Finding 3: Validation targets identified

Three genes (BRCA1, TP53, EGFR) were prioritized for RT-qPCR validation based on fold change and biological relevance.

## Discussion

> **Important:** High duplicate rates in some samples may affect quantification accuracy. Consider additional library preparation QC for future experiments.

The observed upregulation pattern is consistent with cellular stress response pathways, suggesting the treatment induces adaptive metabolic changes.

## Recommendations

- [ ] Complete RT-qPCR validation
- [ ] Prepare figures for manuscript
- [x] Upload raw data to GEO
- [ ] Draft methods section
