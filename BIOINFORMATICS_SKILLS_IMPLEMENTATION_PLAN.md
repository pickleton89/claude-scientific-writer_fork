# Bioinformatics Skills Implementation Plan

> Detailed roadmap for extending the Claude Scientific Writer skills library
> Based on gap analysis and architectural decisions from 2025-12-27
> Target: Papers, grants, figures, and erudite analysis workflows

---

## Executive Summary

This plan extends the skills library to support bioinformatics research workflows while preserving its core identity as a **Science Communication Engine**. The implementation adds **3 new skills** and **4 extensions** to existing skills, focusing on documentation and communication rather than pipeline execution.

### Scope Decisions (Confirmed)

| Question | Decision | Implication |
|----------|----------|-------------|
| Maintain R/Bioconductor docs? | **Yes** | Add R references to plotting-libraries |
| Run or document pipelines? | **Document** | Focus on Methods sections, not execution |
| Primary output? | **Papers/grants/figures** | Prioritize statistical interpretation and visual communication |

### Implementation Overview

| Category | Items | Estimated Files | Priority |
|----------|-------|-----------------|----------|
| New Skills | 3 | ~18 new files | High |
| Skill Extensions | 4 | ~12 new files | Medium |
| Quick Wins | 3 | ~3 new files | Immediate |
| **Total** | **10 work items** | **~33 files** | — |

---

## Phase 0: Quick Wins (Immediate Value)

**Objective**: Deliver immediate value with minimal effort by extending existing skills.

**Timeline**: Can be completed in a single session.

---

### 0.1 Extend `research-lookup` with Bioinformatics Databases

**File**: `skills/research-lookup/references/bioinformatics_databases.md`

**Purpose**: Enable Claude to guide users through programmatic data retrieval from core bioinformatics resources.

**Content Outline**:

```markdown
# Bioinformatics Database Reference

## Sequence Databases
### NCBI Ecosystem
- GenBank: Nucleotide sequences
- RefSeq: Curated reference sequences
- SRA: Raw sequencing data
- GEO: Gene expression data
- API patterns: Entrez E-utilities, EDirect CLI

### Ensembl
- Genome browser navigation
- BioMart for bulk queries
- REST API patterns
- Species assemblies and versioning (GRCh38 vs hg38 notation)

### UCSC Genome Browser
- Track selection and visualization
- Table Browser for data extraction
- bigBed/bigWig format handling

## Protein Databases
### UniProt
- Swiss-Prot vs TrEMBL distinction
- Programmatic access patterns
- ID mapping between databases

### Protein Data Bank (PDB)
- Structure retrieval
- RCSB vs PDBe vs PDBj
- Integration with AlphaFold DB

### AlphaFold Database
- Predicted structure retrieval
- Confidence score interpretation (pLDDT)
- Bulk download patterns

## Pathway & Ontology Databases
### Gene Ontology (GO)
- BP/MF/CC categories
- Evidence codes and their reliability
- Enrichment analysis resources (DAVID, Enrichr, clusterProfiler)

### KEGG
- Pathway maps
- API access (requires subscription for some features)
- KEGG Orthology (KO) identifiers

### Reactome
- Human-focused pathway database
- Visualization tools
- Analysis service API

### MSigDB
- Hallmark gene sets
- C2 curated vs C5 GO sets
- GSEA integration

## Disease & Variant Databases
### ClinVar
- Variant interpretation
- Significance levels
- Submission sources and conflicts

### OMIM
- Mendelian disease catalog
- Gene-phenotype relationships

### GWAS Catalog
- Association data retrieval
- Effect size interpretation

## Common Patterns
### Identifier Mapping
- Cross-reference strategies (UniProt ID mapping, biomaRt)
- Gene symbol standardization (HGNC)
- Version-specific identifiers (Ensembl ENSG with versions)

### Bulk Data Retrieval
- FTP vs API trade-offs
- Rate limiting considerations
- Caching strategies for reproducibility
```

**Success Criteria**: User can ask "How do I retrieve all human kinase structures from PDB?" and receive actionable guidance.

---

### 0.2 Extend `scientific-writing/references/field_terminology.md`

**File**: Extend existing `skills/scientific-writing/references/field_terminology.md`

**Purpose**: Add bioinformatics-specific nomenclature standards for accurate scientific communication.

**Content to Add**:

```markdown
## Bioinformatics Nomenclature

### Gene Naming Conventions

#### HGNC (Human Gene Nomenclature Committee)
- Official human gene symbols: uppercase, italicized in text (e.g., *TP53*, *BRCA1*)
- Protein products: uppercase, not italicized (e.g., TP53, BRCA1)
- Never use aliases in formal writing without defining them

#### Model Organism Conventions
| Organism | Gene | Protein | Example |
|----------|------|---------|---------|
| Human | *GENE* (italic, uppercase) | GENE (roman, uppercase) | *TP53* → TP53 |
| Mouse | *Gene* (italic, sentence case) | Gene (roman, sentence case) | *Trp53* → Trp53 |
| Zebrafish | *gene* (italic, lowercase) | Gene (roman, sentence case) | *tp53* → Tp53 |
| Drosophila | *gene* (italic, lowercase) | Gene (roman, sentence case) | *p53* → P53 |
| C. elegans | *gene-#* (italic, lowercase-number) | GENE-# (roman, uppercase) | *cep-1* → CEP-1 |
| Yeast | *GENE#* (italic, uppercase) | Gene#p (roman, sentence case + p) | *TRP1* → Trp1p |

#### Common Errors to Avoid
- Using unofficial aliases (p53 vs *TP53*)
- Inconsistent capitalization
- Missing italics for genes
- Confusing gene and protein notation

### Variant Nomenclature (HGVS)

#### Sequence Types
- `c.` - coding DNA
- `g.` - genomic DNA
- `m.` - mitochondrial DNA
- `n.` - non-coding RNA
- `p.` - protein

#### Common Variant Notation
| Type | Format | Example |
|------|--------|---------|
| Substitution | c.76A>T | Adenine to Thymine at position 76 |
| Deletion | c.76delA | Deletion of A at position 76 |
| Insertion | c.76_77insT | Insertion of T between positions 76 and 77 |
| Duplication | c.76dupA | Duplication of A at position 76 |
| Frameshift | p.Arg97Profs*23 | Frameshift at Arg97, stop at position 23 |

#### Protein Changes
- Three-letter amino acid codes preferred in publications
- Use `Ter` or `*` for stop codons, never `X`
- Example: p.Arg248Gln (not p.R248Q in formal text)

### Genomic Coordinates

#### Genome Assemblies
- Always specify assembly version: GRCh38 (preferred) or hg38 (UCSC)
- For mouse: GRCm39 or mm39
- Include in Methods: "Coordinates are reported relative to GRCh38"

#### Coordinate Systems
- **1-based, fully-closed** (used by: Ensembl, HGVS, VCF)
- **0-based, half-open** (used by: UCSC, BED, BAM)
- Always state which system when reporting coordinates

### Species Conventions
- First mention: Full name italicized (*Homo sapiens*)
- Subsequent: Abbreviated (*H. sapiens*) or common name
- Never: Non-italicized genus (*Homo sapiens* is wrong as *Homo sapiens*)

### Sequence Data Conventions
- RNA-seq read counts: Report as raw counts for DESeq2/edgeR, or TPM/FPKM for visualization
- Coverage: Report as mean ± SD or median with range
- Quality scores: Phred scale (Q30 = 99.9% accuracy)
```

**Success Criteria**: Manuscripts produced with this skill correctly use *TP53* (gene) vs TP53 (protein) notation.

---

### 0.3 Extend `venue-templates` with Bioinformatics Journals

**File**: `skills/venue-templates/references/bioinformatics_journals.md`

**Purpose**: Provide submission guidelines for major bioinformatics journals.

**Content Outline**:

```markdown
# Bioinformatics Journal Formatting Guide

## Nucleic Acids Research (NAR)

### Article Types
- **Database issue** (annual, October deadline): 4-8 pages
- **Web Server issue** (annual, February deadline): 4-6 pages
- **Methods Online**: New computational methods
- **Standard research articles**: No page limit

### Formatting Requirements
- Abstract: 250 words max, structured optional
- Data availability: MANDATORY section after Conclusion
- Code availability: GitHub/Zenodo DOI required
- Supplementary data: Hosted by journal, not external links

### Key Sections Required
1. Introduction
2. Materials and Methods (detailed, reproducible)
3. Results
4. Discussion
5. Data Availability
6. Supplementary Data
7. Funding
8. Conflict of Interest

### Figure Specifications
- Format: TIFF, EPS, or PDF (vector preferred)
- Resolution: 300 dpi minimum
- Width: Single column (86mm) or double column (178mm)
- Color: Free, but ensure grayscale readability

### Citation Style
- Numbered, in order of appearance
- Format: Author1, Author2, and Author3 (Year) Title. *Journal*, **volume**, pages.

---

## Bioinformatics (Oxford)

### Article Types
- **Original Paper**: Novel methods, 7 pages
- **Application Note**: Software/database, 2 pages
- **Review**: By invitation

### Key Requirements
- Abstract: 150-250 words
- Web link for software REQUIRED
- Supplementary materials for validation data

### Specific Expectations
- Comparison with existing methods mandatory
- Runtime and memory benchmarks expected
- Example datasets with expected outputs

---

## Genome Biology

### Article Types
- **Research**: Major advances, no length limit
- **Method**: New computational/experimental methods
- **Software**: Open-source tools

### Formatting
- Open access (APC ~$3,000)
- Structured abstract required
- Data availability: GEO/SRA accession mandatory

### Review Criteria
- Novelty of biological insight
- Technical rigor
- Broad interest to genomics community

---

## Genome Research

### Focus Areas
- Genome structure and function
- Comparative genomics
- New sequencing technologies

### Requirements
- Resource papers must deposit all data
- Methods must be reproducible
- Code must be openly available

---

## Nature Methods

### Scope
- Transformative methodological advances
- Not incremental improvements

### Requirements
- Validation on multiple datasets
- Comparison with gold standards
- Code/software with documentation
- Reporting summary (Life Sciences Reporting)

---

## PLOS Computational Biology

### Article Types
- **Research Article**: Novel computational methods/analyses
- **Software**: Open-source tools
- **Education**: Tutorials (Ten Simple Rules series)

### Requirements
- Open access
- FAIR data principles
- Open source code (GitHub + archived DOI)

---

## Common Bioinformatics Journal Requirements

### Data Deposition (Universal)
| Data Type | Repository |
|-----------|------------|
| Raw sequencing | SRA, ENA, DDBJ |
| Processed expression | GEO |
| Variants | ClinVar, EVA |
| Structures | PDB, EMDB |
| Proteomics | PRIDE |
| Metabolomics | MetaboLights |
| Code | GitHub + Zenodo DOI |

### Software Availability Checklist
- [ ] Source code on GitHub/GitLab
- [ ] Archived version with DOI (Zenodo)
- [ ] Installation instructions
- [ ] Example data with expected output
- [ ] License clearly stated
- [ ] Version number in manuscript

### Common Rejection Reasons
1. Missing data availability statement
2. No comparison with existing methods
3. Irreproducible results (missing parameters, versions)
4. Overstated claims without validation
5. No biological validation of computational predictions
```

**Success Criteria**: User can ask "What are NAR Database issue requirements?" and receive complete submission guidance.

---

## Phase 1: New Skill — `statistical-analysis`

**Objective**: Enable Claude to guide statistical test selection, generate analysis code, and interpret results for bioinformatics data.

**Priority**: HIGH — Fills critical gap between data and interpretation.

---

### 1.1 Skill Architecture

```
skills/statistical-analysis/
├── SKILL.md                           # Core skill definition
└── references/
    ├── test_decision_framework.md     # Which test to use when
    ├── python_scipy_patterns.md       # scipy.stats code patterns
    ├── r_stats_patterns.md            # R base stats + specialized packages
    ├── multiple_testing_correction.md # FDR, Bonferroni, permutation
    ├── normalization_methods.md       # TMM, DESeq2, limma-voom, batch correction
    ├── dimensionality_reduction.md    # PCA, t-SNE, UMAP theory + code
    ├── clustering_methods.md          # Hierarchical, k-means, Louvain
    ├── survival_analysis.md           # Kaplan-Meier, Cox regression
    └── power_sample_size.md           # Power analysis for omics studies
```

---

### 1.2 SKILL.md Content

```markdown
---
name: statistical-analysis
description: "Statistical methods for scientific data analysis. Test selection decision trees, code generation (Python scipy/statsmodels, R), multiple testing correction, and result interpretation. Use for hypothesis testing, experimental design validation, and quantitative analysis in papers and grants."
---

<objective>
Guide rigorous statistical analysis for scientific research, with emphasis on bioinformatics and high-dimensional data. This skill provides decision frameworks for test selection, code patterns for implementation, and interpretation guidance for manuscripts.
</objective>

<scope>
**In Scope:**
- Test selection based on data characteristics
- Code generation for Python (scipy, statsmodels) and R
- Multiple testing correction strategies
- Effect size calculation and interpretation
- Power analysis and sample size estimation
- Common statistical pitfalls and how to avoid them

**Out of Scope:**
- Machine learning model development (see separate resources)
- Bayesian inference (specialized topic)
- Time series analysis (specialized topic)
</scope>

<decision_framework>
## Quick Test Selection

### Comparing Two Groups

| Data Type | Distribution | Paired? | Test |
|-----------|--------------|---------|------|
| Continuous | Normal | No | Independent t-test |
| Continuous | Normal | Yes | Paired t-test |
| Continuous | Non-normal | No | Mann-Whitney U |
| Continuous | Non-normal | Yes | Wilcoxon signed-rank |
| Categorical | — | No | Chi-square / Fisher's exact |
| Proportions | — | No | Z-test for proportions |

### Comparing Multiple Groups

| Data Type | Distribution | Test | Post-hoc |
|-----------|--------------|------|----------|
| Continuous | Normal | One-way ANOVA | Tukey HSD |
| Continuous | Non-normal | Kruskal-Wallis | Dunn's test |
| Continuous | Normal, repeated | Repeated measures ANOVA | — |
| Continuous | Non-normal, repeated | Friedman test | Nemenyi |

### Relationships

| Variables | Test | Output |
|-----------|------|--------|
| Two continuous | Pearson correlation | r, p-value |
| Two continuous (non-normal) | Spearman correlation | ρ, p-value |
| Continuous outcome, predictors | Linear regression | β, R², p-values |
| Binary outcome, predictors | Logistic regression | OR, p-values |
| Time-to-event | Cox proportional hazards | HR, p-values |

</decision_framework>

<usage_patterns>
## When to Apply This Skill

1. **During analysis planning**: "What test should I use for comparing gene expression between treatment groups?"
2. **Code generation**: "Generate Python code for Mann-Whitney U test with effect size"
3. **Results interpretation**: "Is my p-value of 0.03 significant after testing 20,000 genes?"
4. **Manuscript writing**: "How do I report ANOVA results in a methods section?"
5. **Sample size planning**: "How many samples do I need to detect a 2-fold change in RNA-seq?"

</usage_patterns>

<common_pitfalls>
## Statistical Errors to Avoid

1. **Multiple testing without correction**: Testing 20,000 genes at α=0.05 yields ~1,000 false positives
2. **Pseudoreplication**: Treating technical replicates as biological replicates
3. **Normality assumption violations**: Using t-test on count data or skewed distributions
4. **P-value misinterpretation**: P < 0.05 does not mean 95% probability hypothesis is true
5. **Effect size neglect**: Statistically significant ≠ biologically meaningful
6. **Batch effects**: Confounding technical variation with biological signal
7. **Overfitting in high-dimensional data**: More features than samples

</common_pitfalls>

<cross_references>
- **scientific-critical-thinking**: Conceptual framework for evidence evaluation
- **plotting-libraries**: Visualizing statistical results
- **scientific-writing**: Reporting statistical methods and results
- **reproducible-research**: Documenting analysis parameters
</cross_references>
```

---

### 1.3 Reference File: `test_decision_framework.md`

**Purpose**: Detailed decision trees for test selection.

**Content Structure**:

```markdown
# Statistical Test Decision Framework

## Step 1: Identify Your Question Type

### Comparison Questions
"Is there a difference between groups?"
→ Go to Section 2

### Association Questions
"Is there a relationship between variables?"
→ Go to Section 3

### Prediction Questions
"Can I predict an outcome from predictors?"
→ Go to Section 4

---

## Section 2: Comparison Tests

### 2.1 How Many Groups?

```
                    ┌─────────────────┐
                    │ How many groups?│
                    └────────┬────────┘
                             │
              ┌──────────────┼──────────────┐
              ▼              ▼              ▼
         Two groups    Three+ groups   Before/After
              │              │              │
              ▼              ▼              ▼
         Section 2.2    Section 2.3    Section 2.4
```

### 2.2 Two Independent Groups

**Decision Tree:**
```
Is your outcome variable continuous?
├── Yes → Is the data normally distributed in each group?
│         ├── Yes → Are variances equal? (Levene's test)
│         │         ├── Yes → Independent samples t-test
│         │         └── No → Welch's t-test
│         └── No → Mann-Whitney U test
└── No → Is it categorical?
          ├── Binary → Fisher's exact test (n<20) or Chi-square (n≥20)
          └── Ordinal → Mann-Whitney U test
```

**Normality Assessment:**
1. Visual: Histogram, Q-Q plot
2. Statistical: Shapiro-Wilk test (n < 50), Anderson-Darling (n ≥ 50)
3. Rule of thumb: n > 30 often approximates normal via CLT

**Code Examples:**

Python (scipy):
```python
from scipy import stats

# Normality check
stat, p_normal = stats.shapiro(group1)

# If normal: t-test
t_stat, p_value = stats.ttest_ind(group1, group2)

# If non-normal: Mann-Whitney
u_stat, p_value = stats.mannwhitneyu(group1, group2, alternative='two-sided')

# Effect size (Cohen's d)
cohens_d = (group1.mean() - group2.mean()) / np.sqrt(
    ((len(group1)-1)*group1.std()**2 + (len(group2)-1)*group2.std()**2) /
    (len(group1) + len(group2) - 2)
)
```

R:
```r
# Normality check
shapiro.test(group1)

# If normal: t-test
t.test(group1, group2)

# If non-normal: Mann-Whitney (Wilcoxon rank-sum)
wilcox.test(group1, group2)

# Effect size
library(effectsize)
cohens_d(group1, group2)
```

[Continue with detailed trees for each scenario...]
```

---

### 1.4 Reference File: `multiple_testing_correction.md`

**Purpose**: Critical for omics data analysis.

**Key Sections**:
- Why correction is necessary (family-wise error rate)
- Bonferroni: When and why it's too conservative
- Benjamini-Hochberg FDR: Standard for genomics
- q-value interpretation
- Permutation-based approaches
- Code patterns for both Python and R

---

### 1.5 Reference File: `normalization_methods.md`

**Purpose**: Pre-requisite for valid statistical testing of sequencing data.

**Key Sections**:
- Why raw counts can't be compared directly
- Library size normalization (CPM, RPKM/FPKM, TPM)
- TMM normalization (edgeR)
- DESeq2 size factors
- Limma-voom transformation
- When to use which method
- Batch effect correction (ComBat, Harmony, limma removeBatchEffect)

---

### 1.6 Remaining Reference Files

| File | Purpose | Key Content |
|------|---------|-------------|
| `python_scipy_patterns.md` | Code reference | Complete scipy.stats API patterns with examples |
| `r_stats_patterns.md` | Code reference | Base R stats + common packages (car, lme4) |
| `dimensionality_reduction.md` | Theory + code | PCA interpretation, t-SNE perplexity, UMAP parameters |
| `clustering_methods.md` | Method selection | When to use hierarchical vs k-means vs graph-based |
| `survival_analysis.md` | Time-to-event | Kaplan-Meier plots, log-rank test, Cox regression |
| `power_sample_size.md` | Study design | Power calculations for RNA-seq, proteomics |

---

## Phase 2: New Skill — `code-documentation`

**Objective**: Enable Claude to generate high-quality documentation for analysis code, following community standards.

**Priority**: HIGH — Bridges code and communication.

---

### 2.1 Skill Architecture

```
skills/code-documentation/
├── SKILL.md                         # Core skill definition
└── references/
    ├── python_docstrings.md         # Google, NumPy, Sphinx styles
    ├── r_documentation.md           # roxygen2, pkgdown
    ├── readme_templates.md          # Analysis project READMEs
    ├── notebook_best_practices.md   # Jupyter, R Markdown, Quarto
    ├── inline_commenting.md         # When and how to comment
    └── api_documentation.md         # For reusable code/packages
```

---

### 2.2 SKILL.md Content

```markdown
---
name: code-documentation
description: "Documentation for scientific code. Docstrings (Python Google/NumPy, R roxygen2), README generation, notebook best practices, and inline commenting. Use when documenting analysis scripts, creating project READMEs, or writing reproducible notebooks."
---

<objective>
Generate clear, consistent documentation for scientific code that serves both human readers and automated tools. This skill covers function-level documentation (docstrings), project-level documentation (READMEs), and narrative documentation (notebooks).
</objective>

<documentation_hierarchy>
## Levels of Documentation

### 1. Inline Comments
- **When**: Complex logic that isn't self-evident
- **Style**: Brief, explains "why" not "what"
- **Anti-pattern**: Commenting obvious code (`i += 1  # increment i`)

### 2. Docstrings
- **When**: Every public function, class, module
- **Style**: Structured format (Google, NumPy, or roxygen2)
- **Content**: Purpose, parameters, returns, examples

### 3. README
- **When**: Every project/repository
- **Content**: Purpose, installation, usage, examples, citation

### 4. Notebooks
- **When**: Analysis workflows, tutorials, exploratory work
- **Style**: Narrative flow with code cells

### 5. API Documentation
- **When**: Reusable packages/libraries
- **Tools**: Sphinx (Python), pkgdown (R)
</documentation_hierarchy>

<style_selection>
## Which Docstring Style?

| Context | Recommended Style | Why |
|---------|-------------------|-----|
| Scientific Python | NumPy style | Standard in numpy, scipy, pandas ecosystem |
| General Python | Google style | Clean, readable, less verbose |
| Sphinx documentation | NumPy or Google | Both have Sphinx extensions |
| R packages | roxygen2 | Required for CRAN submission |
| Quick scripts | Minimal | One-line description sufficient |

</style_selection>

<cross_references>
- **reproducible-research**: Project structure and environment documentation
- **scientific-writing**: Methods section describing analysis
- **markitdown**: Converting documentation formats
</cross_references>
```

---

### 2.3 Reference File: `python_docstrings.md`

**Content Structure**:

```markdown
# Python Docstring Reference

## NumPy Style (Recommended for Scientific Code)

### Function Template

```python
def calculate_fold_change(treatment: np.ndarray,
                          control: np.ndarray,
                          log2: bool = True,
                          pseudocount: float = 1.0) -> np.ndarray:
    """
    Calculate fold change between treatment and control conditions.

    Computes the ratio of treatment to control values, optionally
    log2-transformed. Adds pseudocount to avoid division by zero
    and log of zero.

    Parameters
    ----------
    treatment : np.ndarray
        Expression values for treatment condition. Shape (n_genes,)
        or (n_genes, n_samples).
    control : np.ndarray
        Expression values for control condition. Must match treatment shape.
    log2 : bool, optional
        If True, return log2 fold change. Default is True.
    pseudocount : float, optional
        Value added before division/log. Default is 1.0.

    Returns
    -------
    np.ndarray
        Fold change values. Same shape as input arrays.

    Raises
    ------
    ValueError
        If treatment and control shapes don't match.

    See Also
    --------
    scipy.stats.ttest_ind : Statistical test for fold change significance.

    Notes
    -----
    Log2 fold change of 1 means 2-fold increase.
    Log2 fold change of -1 means 2-fold decrease.

    Examples
    --------
    >>> treatment = np.array([10, 20, 5])
    >>> control = np.array([5, 10, 10])
    >>> calculate_fold_change(treatment, control)
    array([ 1., 1., -1.])
    """
    if treatment.shape != control.shape:
        raise ValueError(f"Shape mismatch: {treatment.shape} vs {control.shape}")

    fc = (treatment + pseudocount) / (control + pseudocount)

    if log2:
        return np.log2(fc)
    return fc
```

### Class Template

```python
class DifferentialExpression:
    """
    Differential expression analysis container.

    Stores results from differential expression analysis and provides
    methods for filtering, visualization, and export.

    Parameters
    ----------
    gene_ids : array-like
        Gene identifiers (Ensembl IDs or symbols).
    log2fc : array-like
        Log2 fold change values.
    pvalues : array-like
        Raw p-values from statistical test.
    padj : array-like, optional
        Adjusted p-values (FDR). Computed if not provided.

    Attributes
    ----------
    results : pd.DataFrame
        Combined results with all columns.
    n_genes : int
        Total number of genes tested.
    n_significant : int
        Number of genes passing default thresholds.

    Examples
    --------
    >>> de = DifferentialExpression(
    ...     gene_ids=['BRCA1', 'TP53', 'MYC'],
    ...     log2fc=[2.5, -1.2, 0.1],
    ...     pvalues=[0.001, 0.02, 0.8]
    ... )
    >>> de.get_significant(padj_threshold=0.05, lfc_threshold=1.0)
    """
```

## Google Style

[Include Google style examples...]

## Quick Reference Card

| Section | Required? | Content |
|---------|-----------|---------|
| Short summary | Yes | One line, imperative mood |
| Extended summary | If needed | Paragraph explaining details |
| Parameters | If any | Name, type, description |
| Returns | If any | Type, description |
| Raises | If any | Exception types and when |
| Examples | Recommended | Doctest-compatible |
| See Also | Optional | Related functions |
| Notes | Optional | Implementation details |
| References | Optional | Citations |
```

---

### 2.4 Reference File: `readme_templates.md`

**Purpose**: Templates for analysis project READMEs.

**Key Templates**:
1. Bioinformatics analysis project README
2. Software tool README
3. Data repository README

Each template includes:
- Badge suggestions (CI status, license, DOI)
- Standard sections (Installation, Usage, Data, Citation)
- Example content

---

### 2.5 Reference File: `notebook_best_practices.md`

**Key Sections**:
- Narrative structure (introduction, methods, results, conclusions)
- Code cell organization (imports first, one concept per cell)
- Markdown formatting for scientific notebooks
- Output management (hide code for reports, show for methods)
- Reproducibility (seeds, version printing, environment capture)
- Jupyter vs R Markdown vs Quarto comparison

---

## Phase 3: New Skill — `reproducible-research`

**Objective**: Enable Claude to guide FAIR data practices, environment management, and reproducibility documentation.

**Priority**: MEDIUM — Critical for modern bioinformatics but less immediate than statistical analysis.

---

### 3.1 Skill Architecture

```
skills/reproducible-research/
├── SKILL.md                          # Core skill definition
└── references/
    ├── conda_environments.md         # Conda/Mamba patterns
    ├── docker_containers.md          # Containerization for reproducibility
    ├── fair_principles.md            # Findable, Accessible, Interoperable, Reusable
    ├── data_repositories.md          # GEO, SRA, Zenodo, Figshare
    ├── metadata_standards.md         # MIAME, MINSEQE, MIQE
    ├── project_structure.md          # Directory organization
    ├── version_control_data.md       # DVC, git-annex, git-lfs
    └── workflow_documentation.md     # Snakemake/Nextflow documentation (not execution)
```

---

### 3.2 SKILL.md Content

```markdown
---
name: reproducible-research
description: "Reproducibility and FAIR data practices for scientific research. Environment management (Conda, Docker), data repositories (GEO, SRA, Zenodo), metadata standards (MIAME, MINSEQE), and project organization. Use for Methods sections, Data Availability statements, and analysis documentation."
---

<objective>
Guide the documentation and infrastructure needed for reproducible scientific computing. This skill focuses on what to document and where to deposit data—not on executing analyses, but on making them reproducible by others.
</objective>

<scope>
**In Scope:**
- Environment specification (Conda environment.yml, requirements.txt, Dockerfiles)
- Data deposition guidance (which repository, what metadata)
- Project structure templates
- FAIR principles implementation
- Data Availability statement writing
- Code Availability statement writing

**Out of Scope:**
- Executing Snakemake/Nextflow pipelines (operational)
- Setting up CI/CD systems (DevOps)
- Cloud infrastructure (AWS, GCP)
</scope>

<fair_quick_reference>
## FAIR Principles Summary

### Findable
- Assign persistent identifiers (DOI)
- Rich metadata describing the data
- Registered in searchable resource

### Accessible
- Retrievable by identifier using standard protocol
- Metadata accessible even if data restricted
- Authentication/authorization where needed

### Interoperable
- Use formal, shared language (ontologies)
- Use FAIR vocabularies
- Include qualified references to other data

### Reusable
- Clear usage license
- Detailed provenance
- Meet domain-relevant standards
</fair_quick_reference>

<data_deposition_guide>
## Where to Deposit Data

| Data Type | Primary Repository | Alternative |
|-----------|-------------------|-------------|
| Raw sequencing (FASTQ) | SRA (NCBI) | ENA (EBI), DDBJ |
| Processed expression | GEO | ArrayExpress |
| Variants | ClinVar, EVA | dbSNP, dbVar |
| Proteomics | PRIDE | MassIVE |
| Metabolomics | MetaboLights | Metabolomics Workbench |
| Structures | PDB | EMDB (cryo-EM) |
| General/Other | Zenodo | Figshare, Dryad |
| Code | GitHub + Zenodo | GitLab, Bitbucket + archive |

</data_deposition_guide>

<cross_references>
- **scientific-writing**: Data Availability statement text
- **code-documentation**: README and environment docs
- **venue-templates**: Journal-specific deposition requirements
</cross_references>
```

---

### 3.3 Key Reference Files

#### `conda_environments.md`
- `environment.yml` structure
- Pinning versions for reproducibility
- Cross-platform considerations
- Mamba for faster solving
- Exporting and recreating environments

#### `data_repositories.md`
- Step-by-step GEO submission guide
- SRA submission workflow
- Zenodo for code archiving
- DOI minting and citation

#### `metadata_standards.md`
- MIAME for microarray data
- MINSEQE for sequencing data
- MIQE for qPCR data
- STORMS for metagenomics
- Checklist templates for each

---

## Phase 4: Extend `plotting-libraries`

**Objective**: Add R ecosystem and bioinformatics-specific visualizations.

**Priority**: MEDIUM — Complements existing Python coverage.

---

### 4.1 New Reference Files

```
skills/plotting-libraries/
└── references/
    ├── matplotlib.md          # EXISTS
    ├── seaborn.md             # EXISTS
    ├── ggplot2.md             # NEW: R visualization
    ├── bioconductor_viz.md    # NEW: ComplexHeatmap, EnhancedVolcano
    ├── volcano_plots.md       # NEW: Differential expression visualization
    ├── heatmaps.md            # NEW: Expression heatmaps with clustering
    ├── survival_curves.md     # NEW: Kaplan-Meier visualization
    └── genome_tracks.md       # NEW: pyGenomeTracks, Gviz
```

---

### 4.2 Reference File: `ggplot2.md`

**Content Structure**:

```markdown
# ggplot2 Reference for Scientific Figures

## Grammar of Graphics Fundamentals

### Basic Structure
```r
ggplot(data, aes(x = var1, y = var2)) +
  geom_point() +
  theme_minimal()
```

### Core Components
| Component | Purpose | Example |
|-----------|---------|---------|
| `data` | The dataset | `ggplot(df, ...)` |
| `aes()` | Aesthetic mappings | `aes(x = time, y = value, color = group)` |
| `geom_*` | Geometric objects | `geom_point()`, `geom_line()`, `geom_boxplot()` |
| `stat_*` | Statistical transformations | `stat_smooth()`, `stat_summary()` |
| `scale_*` | Axis and legend scales | `scale_color_manual()`, `scale_y_log10()` |
| `facet_*` | Multi-panel layouts | `facet_wrap(~group)`, `facet_grid(row~col)` |
| `theme_*` | Visual appearance | `theme_minimal()`, `theme_classic()` |

## Publication-Ready Figures

### Theme for Publications
```r
theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
  theme(
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", size = 0.2)
  )
}
```

### Saving for Publication
```r
ggsave("figure1.pdf", width = 7, height = 5, units = "in", dpi = 300)
ggsave("figure1.tiff", width = 180, height = 120, units = "mm", dpi = 300)
```

## Common Bioinformatics Plots in ggplot2

### Expression Boxplot
```r
ggplot(expr_data, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Condition", y = "Expression (log2 TPM)") +
  theme_publication()
```

### PCA Plot
```r
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  labs(x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_publication()
```

[Continue with more examples...]
```

---

### 4.3 Reference File: `volcano_plots.md`

**Purpose**: Standard visualization for differential expression.

**Content Includes**:
- Both Python (matplotlib/seaborn) and R (ggplot2, EnhancedVolcano) implementations
- Significance thresholds (padj < 0.05, |log2FC| > 1)
- Gene labeling strategies
- Color schemes for up/down/non-significant
- Publication formatting

---

### 4.4 Reference File: `heatmaps.md`

**Purpose**: Expression heatmaps with hierarchical clustering.

**Content Includes**:
- Python: seaborn clustermap, matplotlib imshow
- R: ComplexHeatmap (preferred), pheatmap
- Clustering method selection
- Color palette selection (diverging for expression)
- Row/column annotation
- Scaling and centering decisions

---

## Phase 5: Update SKILL.md Files with Cross-References

**Objective**: Ensure skills reference each other appropriately.

---

### 5.1 Skills to Update

| Skill | Add Cross-Reference To |
|-------|------------------------|
| `scientific-writing` | `statistical-analysis` (for methods), `reproducible-research` (for data availability) |
| `hypothesis-generation` | `statistical-analysis` (for power analysis) |
| `scientific-critical-thinking` | `statistical-analysis` (for quantitative evaluation) |
| `plotting-libraries` | `statistical-analysis` (for result visualization) |
| `research-lookup` | Bioinformatics databases reference |
| `peer-review` | `statistical-analysis` (for evaluating stats in papers) |

---

## Implementation Timeline

### Week 1: Quick Wins (Phase 0) ✅ COMPLETED
- [x] Create `research-lookup/references/bioinformatics_databases.md`
- [x] Extend `scientific-writing/references/field_terminology.md`
- [x] Create `venue-templates/references/bioinformatics_journals.md`
- [x] Commit and update CHANGELOG

### Week 2: Statistical Analysis Skill (Phase 1) ✅ COMPLETED
- [x] Create `statistical-analysis/SKILL.md`
- [x] Create `references/test_decision_framework.md`
- [x] Create `references/multiple_testing_correction.md`
- [x] Create `references/normalization_methods.md`
- [x] Create remaining reference files (python_scipy_patterns, r_stats_patterns, dimensionality_reduction, clustering_methods, survival_analysis, power_sample_size)
- [x] Commit and update CHANGELOG

### Week 3: Code Documentation Skill (Phase 2) ✅ COMPLETED
- [x] Create `code-documentation/SKILL.md`
- [x] Create `references/python_docstrings.md`
- [x] Create `references/r_documentation.md`
- [x] Create `references/readme_templates.md`
- [x] Create `references/notebook_best_practices.md`
- [x] Create `references/inline_commenting.md`
- [x] Create `references/api_documentation.md`
- [x] Commit and update CHANGELOG

### Week 4: Reproducible Research Skill (Phase 3)
- [ ] Create `reproducible-research/SKILL.md`
- [ ] Create core reference files (FAIR, repositories, metadata)
- [ ] Create environment reference files (Conda, Docker)
- [ ] Commit and update CHANGELOG

### Week 5: Plotting Extensions (Phase 4)
- [ ] Create `plotting-libraries/references/ggplot2.md`
- [ ] Create `plotting-libraries/references/bioconductor_viz.md`
- [ ] Create specialized plot references (volcano, heatmap, survival)
- [ ] Update SKILL.md cross-references (Phase 5)
- [ ] Final commit and CHANGELOG update

---

## Success Criteria

### Functional Criteria
1. **Statistical Analysis**: Claude correctly recommends Mann-Whitney U for non-normal, unpaired continuous data
2. **Code Documentation**: Claude generates valid NumPy-style docstrings with examples
3. **Reproducible Research**: Claude produces complete Data Availability statements with repository names
4. **Plotting Extensions**: Claude provides both Python and R code for volcano plots

### Quality Criteria
1. All reference files follow consistent markdown structure
2. Code examples are tested and functional
3. Cross-references between skills are bidirectional
4. CHANGELOG documents all additions

### Integration Criteria
1. New skills appear in skill listings
2. Existing skills reference new skills appropriately
3. No duplicate content between skills
4. Clear scope boundaries prevent overlap

---

## Appendix: File Checklist

### New Files (Total: ~33)

**Phase 0: Quick Wins (3 files)** ✅ COMPLETED
- [x] `skills/research-lookup/references/bioinformatics_databases.md`
- [x] `skills/scientific-writing/references/field_terminology.md` (extend)
- [x] `skills/venue-templates/references/bioinformatics_journals.md`

**Phase 1: statistical-analysis (10 files)** ✅ COMPLETED
- [x] `skills/statistical-analysis/SKILL.md`
- [x] `skills/statistical-analysis/references/test_decision_framework.md`
- [x] `skills/statistical-analysis/references/multiple_testing_correction.md`
- [x] `skills/statistical-analysis/references/normalization_methods.md`
- [x] `skills/statistical-analysis/references/python_scipy_patterns.md`
- [x] `skills/statistical-analysis/references/r_stats_patterns.md`
- [x] `skills/statistical-analysis/references/dimensionality_reduction.md`
- [x] `skills/statistical-analysis/references/clustering_methods.md`
- [x] `skills/statistical-analysis/references/survival_analysis.md`
- [x] `skills/statistical-analysis/references/power_sample_size.md`

**Phase 2: code-documentation (7 files)** ✅ COMPLETED
- [x] `skills/code-documentation/SKILL.md`
- [x] `skills/code-documentation/references/python_docstrings.md`
- [x] `skills/code-documentation/references/r_documentation.md`
- [x] `skills/code-documentation/references/readme_templates.md`
- [x] `skills/code-documentation/references/notebook_best_practices.md`
- [x] `skills/code-documentation/references/inline_commenting.md`
- [x] `skills/code-documentation/references/api_documentation.md`

**Phase 3: reproducible-research (8 files)**
- [ ] `skills/reproducible-research/SKILL.md`
- [ ] `skills/reproducible-research/references/conda_environments.md`
- [ ] `skills/reproducible-research/references/docker_containers.md`
- [ ] `skills/reproducible-research/references/fair_principles.md`
- [ ] `skills/reproducible-research/references/data_repositories.md`
- [ ] `skills/reproducible-research/references/metadata_standards.md`
- [ ] `skills/reproducible-research/references/project_structure.md`
- [ ] `skills/reproducible-research/references/version_control_data.md`
- [ ] `skills/reproducible-research/references/workflow_documentation.md`

**Phase 4: plotting-libraries extensions (6 files)**
- [ ] `skills/plotting-libraries/references/ggplot2.md`
- [ ] `skills/plotting-libraries/references/bioconductor_viz.md`
- [ ] `skills/plotting-libraries/references/volcano_plots.md`
- [ ] `skills/plotting-libraries/references/heatmaps.md`
- [ ] `skills/plotting-libraries/references/survival_curves.md`
- [ ] `skills/plotting-libraries/references/genome_tracks.md`

---

*Document generated: 2025-12-27*
*Based on gap analysis: skills_library_gap_analysis.md*
