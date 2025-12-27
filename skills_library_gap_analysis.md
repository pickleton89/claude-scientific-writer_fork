# Bioinformatics Skills Library: Gap Analysis and Roadmap

## Executive Summary cur 

Your current skills library is an exceptional **Science Communication Engine**, excelling at the "tail end" of the scientific process—transforming results into papers, posters, slides, and diagrams. However, for a bioinformatic scientist, the critical blind spot is the **Execution Phase**: the bridge between raw data and the final figure. There is a significant gap between `hypothesis-generation` (the idea) and `plotting-libraries` (the visualization).

This analysis identifies the missing skill areas needed to transform your library from a **Writer** into a full **Research Assistant** that supports the complete bioinformatics workflow.

---

## Current Strengths

Your library demonstrates strong coverage across several domains:

| Domain | Coverage | Supporting Skills |
|--------|----------|-------------------|
| Scientific Writing | ✅ Strong | scientific-writing, literature-review, hypothesis-generation |
| Presentations | ✅ Strong | scientific-slides, latex-posters, pptx-posters |
| Literature & Research | ✅ Strong | research-lookup, citation-management, peer-review, scholar-evaluation |
| Visualization | ✅ Moderate | plotting-libraries, scientific-schematics, visual-design |
| Document Processing | ✅ Strong | markitdown, venue-templates, document-skills/* |
| Critical Evaluation | ✅ Strong | scientific-critical-thinking (methodology critique, evidence quality) |

---

## Identified Gaps

### Critical Gaps (No Current Coverage)

#### 1. Sequence & Omics Analysis

No skills exist for the core work of bioinformatics:

- Sequence analysis (DNA/RNA/protein)
- Genomics pipelines (alignment, variant calling, annotation)
- Transcriptomics (RNA-seq, differential expression, scRNA-seq)
- Proteomics and metabolomics workflows
- Common file formats: FASTA, FASTQ, BAM/SAM, VCF, GFF/GTF, BED

#### 2. Bioinformatics Tools & Pipelines

No reference exists for standard tools:

- **Aligners:** BWA, STAR, Bowtie2, minimap2
- **Variant callers:** GATK, bcftools, DeepVariant
- **Assemblers:** SPAdes, Flye, MEGAHIT
- **Analysis suites:** samtools, bedtools, VEP
- **Workflow managers:** Snakemake, Nextflow, WDL

#### 3. Bio-Data Wrangling

You have `plotting-libraries` to visualize data, but bioinformatic data (FASTQ, BAM, VCF, PDB) requires heavy cleaning and specific parsing before visualization. This gap represents the "glue code" that consumes 80% of a bioinformatician's time.

**Required capabilities:**
- Format parsing with Biopython (Python) and Bioconductor (R)
- Tidy data transformations using Pandas/Polars or dplyr
- Sanity checking for NaNs, sample mismatches, and gene name standardization (Ensembl vs. Symbol)

### Major Gaps (Partial or Weak Coverage)

#### 4. Statistical Methods for High-Dimensional Data

`scientific-critical-thinking` evaluates logic, but you lack a skill to generate actual statistical testing code. A chart shows a trend; statistics prove it.

**Missing capabilities:**
- Test selection decision trees (t-test vs. Wilcoxon, ANOVA vs. Kruskal-Wallis)
- Normalization methods (TMM, DESeq2, limma-voom)
- Multiple testing correction (FDR, Benjamini-Hochberg, Bonferroni, permutation)
- Dimensionality reduction (PCA, t-SNE, UMAP)
- Clustering methods (hierarchical, k-means, Louvain)
- Batch effect correction (ComBat, Harmony)
- Survival analysis
- Power analysis and sample size calculation

#### 5. R & Bioconductor Ecosystem

`plotting-libraries` is Python-only, but bioinformatics relies heavily on R:

- Bioconductor packages (DESeq2, edgeR, clusterProfiler, Seurat)
- ggplot2 for publication figures
- R Markdown/Quarto for reproducible reports
- tidyverse for data wrangling

#### 6. Biological Databases & Resources

No skill exists for navigating essential databases:

- **Sequence:** NCBI (GenBank, SRA, GEO), Ensembl, UCSC
- **Protein:** UniProt, PDB, AlphaFold DB
- **Pathways:** KEGG, Reactome, GO, MSigDB
- **Disease:** OMIM, ClinVar, GWAS Catalog
- API access patterns for programmatic retrieval

#### 7. Reproducible Research & Data Management

In bioinformatics, if the code isn't reproducible, the science isn't valid. No skill is dedicated to environment or pipeline management.

**Missing capabilities:**
- Environment configuration: `environment.yml` (Conda) or Dockerfile definitions
- Project structure: Cookiecutter data science directory layouts
- Workflow managers: Nextflow or Snakemake scaffolds
- Data versioning (DVC, git-annex)
- FAIR principles implementation
- Repository submission (GEO, SRA, Zenodo, Figshare)
- Metadata standards (MIAME, MINSEQE)

### Moderate Gaps

#### 8. Domain-Specific Visualization

Missing bioinformatics-specific plot types:

- Genome browsers / track plots (IGV, pyGenomeTracks)
- Heatmaps with clustering dendrograms
- Volcano plots, MA plots, dispersion plots
- Phylogenetic trees
- Network/interaction graphs
- Circos/chord diagrams for genomic features

#### 9. Code Documentation

`scientific-writing` focuses on prose for humans (journals). Bioinformaticians also write for machines and other developers.

**Missing capabilities:**
- Docstring generation: Python (Google/NumPy style) and R (roxygen2)
- README generation with installation, usage, and badges
- Beginner-friendly inline commenting
- Notebook best practices (Jupyter, R Markdown, Quarto)

### Weak Areas (Existing Skills with Gaps)

| Current Skill | Gap |
|---------------|-----|
| scientific-writing | Needs bioinformatics depth: HGVS nomenclature, gene naming conventions, computational methods sections |
| reporting_guidelines.md | Clinical focus (CONSORT/STROBE). Missing: MIAME (arrays), MINSEQE (sequencing), MIQE (qPCR), MIxS (metagenomics), STORMS |
| hypothesis-generation | Generic experimental design. Missing: power analysis for omics, sample size for RNA-seq, batch design |
| scientific-critical-thinking | Missing omics-specific statistical pitfalls (normalization artifacts, confounders in expression data) |

---

## Recommended New Skills (Priority Order)

| Priority | Skill Name | Description | Integration |
|----------|------------|-------------|-------------|
| 1 | **bio-data-wrangling** | Parse, clean, and reshape biological data into tidy formats; handle FASTQ/BAM/VCF/PDB | Feeds directly into plotting-libraries |
| 2 | **statistical-analysis** | Test selection, code generation (scipy/statsmodels/DESeq2), FDR correction | Results feed into scientific-writing |
| 3 | **r-bioconductor** | R for bioinformatics, Bioconductor packages, ggplot2, tidyverse | Parallel to Python plotting |
| 4 | **bioinformatics-pipelines** | Standard workflows (RNA-seq, WGS, variant calling), tool selection | Core analysis capability |
| 5 | **biological-databases** | Navigate NCBI/Ensembl/UniProt, programmatic access, data retrieval | Supports literature-review |
| 6 | **reproducible-research** | Conda, containers, workflow managers, project structure | Generates Data/Code Availability statements |
| 7 | **bioinformatics-visualization** | Domain-specific plots (volcano, heatmaps, genome tracks, networks) | Extends plotting-libraries |
| 8 | **code-documentation** | Docstrings, READMEs, beginner-friendly comments | Works with markitdown |
| 9 | **sequence-analysis** | Alignment, assembly, annotation, file formats | Foundation for pipelines |

---

## Quick Wins (Extensions to Existing Skills)

These additions provide immediate value with minimal effort:

1. **Add to plotting-libraries:** Bioinformatics plot types, seaborn clustermap, volcano plots, survival curves
2. **Add to scientific-writing/field_terminology.md:** Gene nomenclature (HGNC), variant notation (HGVS), species conventions
3. **Add to reporting_guidelines.md:** MIAME, MINSEQE, MIQE, STORMS (metagenomics)
4. **Add to venue-templates:** Bioinformatics journals (Nucleic Acids Research, Bioinformatics, Genome Biology)

---

## Recommended Starting Point

Given your current learning trajectory in R and Python, start with **`bio-data-wrangling`**. This skill will provide the most immediate help in your daily workflow by helping you write the "glue code" that consumes the majority of a bioinformatician's time. It bridges the gap between your raw data and your visualization capabilities.

---

## The New Workflow Map

With these skills implemented, your library transforms from a Writer to a full Research Assistant:

```
┌─────────────────────────┐
│  hypothesis-generation  │
└───────────┬─────────────┘
            │
            ▼
      ┌─────────────┐
      │  Experiment │
      └──────┬──────┘
             │
             ▼
┌────────────────────────┐
│   bio-data-wrangling   │
└───────────┬────────────┘
            │
     ┌──────┴──────┐
     ▼             ▼
┌──────────┐  ┌─────────────────┐
│ plotting │  │   statistical   │
│ libraries│  │    analysis     │
└────┬─────┘  └───────┬─────────┘
     │                │
     └───────┬────────┘
             ▼
┌────────────────────────┐
│  reproducible-research │
└───────────┬────────────┘
            │
     ┌──────┴──────┐
     ▼             ▼
┌──────────────┐  ┌──────────────────┐
│  scientific  │  │ code-documentation│
│   writing    │  └──────────────────┘
└──────┬───────┘
       │
       ▼
┌──────────────────────┐
│ slides / posters     │
└──────────────────────┘
```

---

## Summary

Your library excels at scientific communication but lacks the computational infrastructure for the execution phase of bioinformatics research. The priority path forward is:

1. **Immediate:** Start with `bio-data-wrangling` to address daily workflow needs
2. **Short-term:** Add `statistical-analysis` and `r-bioconductor` to enable rigorous quantitative work
3. **Medium-term:** Build out `reproducible-research` and pipeline skills for production-ready analyses
4. **Ongoing:** Extend existing skills with bioinformatics-specific content (quick wins)

This roadmap transforms your skills library into a cohesive ecosystem that supports the complete journey from hypothesis to publication.
