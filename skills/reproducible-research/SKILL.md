---
name: reproducible-research
version: 2.1.0
description: "Reproducibility and FAIR data practices with quantified completeness thresholds. Environment management (Conda, Docker), data repositories (GEO, SRA, Zenodo), metadata standards (MIAME, MINSEQE), and project organization. Use for Methods sections, Data Availability statements, and analysis documentation."
shared-thresholds: "../QUANTIFICATION_THRESHOLDS.md"
---

> **Quantified Thresholds:** This skill references shared thresholds from [`QUANTIFICATION_THRESHOLDS.md`](../QUANTIFICATION_THRESHOLDS.md) §5 (Documentation Completeness) and §8 (Iteration & Stopping Criteria).

<objective>
Guide the documentation and infrastructure needed for reproducible scientific computing with quantified completeness levels. Achieve ≥75% reproducibility checklist compliance for publication-ready projects. This skill focuses on what to document and where to deposit data—not on executing analyses, but on making them reproducible by others.
</objective>

<scope>
**In Scope:**
- Environment specification (Conda environment.yml, requirements.txt, Dockerfiles)
- Data deposition guidance (which repository, what metadata)
- Project structure templates
- FAIR principles implementation
- Data Availability statement writing
- Code Availability statement writing
- Metadata standards compliance (MIAME, MINSEQE, MIQE)
- Version control for data and large files

**Out of Scope:**
- Executing Snakemake/Nextflow pipelines (operational)
- Setting up CI/CD systems (DevOps)
- Cloud infrastructure (AWS, GCP)
- Pipeline development (focus is on documentation)
</scope>

<fair_quick_reference>
## FAIR Principles Summary

### FAIR Compliance Scoring

| Principle | Minimum Requirements | Score |
|-----------|---------------------|-------|
| **Findable** | DOI assigned, ≥5 metadata fields, indexed in repository | 25% |
| **Accessible** | Standard protocol (HTTP/FTP), auth documented if restricted | 25% |
| **Interoperable** | Standard format, ontology terms used, qualified references | 25% |
| **Reusable** | License specified, provenance documented, domain standards met | 25% |

**FAIR Compliance Thresholds:**
- **Non-compliant:** <50% (fails ≥2 principles)
- **Partially compliant:** 50-74% (1 principle incomplete)
- **Compliant:** 75-99% (all principles met at minimum)
- **Exemplary:** 100% (exceeds all requirements)

### Findable (25%)
- [ ] Persistent identifier assigned (DOI) — **Required**
- [ ] Rich metadata: title, authors, date, description, keywords (≥5 fields)
- [ ] Registered in searchable resource (repository, catalog)

### Accessible (25%)
- [ ] Retrievable by identifier using standard protocol (HTTP, FTP)
- [ ] Metadata accessible even if data restricted
- [ ] Authentication/authorization documented where needed

### Interoperable (25%)
- [ ] Use formal, shared language (ontologies: GO, SO, CHEBI)
- [ ] Standard file formats (FASTQ, BAM, CSV, not proprietary)
- [ ] Include qualified references to related datasets

### Reusable (25%)
- [ ] Clear usage license specified (CC-BY, MIT, GPL) — **Required**
- [ ] Detailed provenance (how data was generated)
- [ ] Meet domain-relevant standards (MIAME, MINSEQE, MIQE)

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

<environment_documentation>
## Environment Specification Hierarchy

### Quantified Reproducibility Levels

| Level | Name | Requirements | Reproducibility Score |
|-------|------|--------------|----------------------|
| 1 | Minimal | requirements.txt with `>=` versions | 25-49% |
| 2 | Standard | environment.yml with pinned versions | 50-74% |
| 3 | Complete | Docker/Singularity container | 75-89% |
| 4 | Exemplary | Container + CI tests + archived DOI | 90-100% |

**Scoring Thresholds** (see [`QUANTIFICATION_THRESHOLDS.md`](../QUANTIFICATION_THRESHOLDS.md) §5):
- **Not reproducible:** <50% checklist items
- **Partially reproducible:** 50-74%
- **Reproducible (publication-ready):** 75-89%
- **Highly reproducible:** 90-100%

### Level 1: Minimal (requirements.txt)
```
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
```
- Quick start for collaborators
- Insufficient for exact reproducibility (score: 25-49%)
- **Version pinning:** 0% exact pins

### Level 2: Standard (environment.yml)
```yaml
name: analysis-env
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.10.12
  - numpy=1.24.3
  - pandas=2.0.2
```
- Pinned versions (score: 50-74%)
- Cross-platform (mostly)
- **Version pinning:** 100% exact pins required

### Level 3: Complete (Docker/Singularity)
```dockerfile
FROM continuumio/miniconda3:23.5.2-0
COPY environment.yml .
RUN conda env create -f environment.yml
```
- Fully reproducible (score: 75-89%)
- Required for HPC/cloud
- **Validation:** Must pass `docker build` + `docker run` without errors

### Level 4: Exemplary (Archived + Tested)
- Container image archived with DOI
- CI pipeline validates reproducibility
- Produces identical output on 3 consecutive runs
- Score: 90-100%

</environment_documentation>

<project_structure>
## Standard Project Layout

```
project/
├── README.md              # Project overview, installation, usage
├── LICENSE                # Usage terms (MIT, GPL, CC-BY, etc.)
├── environment.yml        # Conda environment specification
├── Dockerfile             # Container definition (optional)
│
├── data/
│   ├── raw/               # Original, immutable data
│   ├── processed/         # Cleaned, transformed data
│   └── external/          # Data from other sources
│
├── notebooks/             # Jupyter/R Markdown notebooks
│   ├── 01_exploration.ipynb
│   ├── 02_analysis.ipynb
│   └── 03_figures.ipynb
│
├── src/                   # Source code
│   ├── __init__.py
│   ├── data_processing.py
│   ├── analysis.py
│   └── visualization.py
│
├── results/
│   ├── figures/           # Publication figures
│   └── tables/            # Summary tables
│
├── docs/                  # Additional documentation
│   ├── methods.md         # Detailed methods
│   └── data_dictionary.md # Variable definitions
│
└── tests/                 # Unit tests (optional)
```

</project_structure>

<usage_patterns>
## When to Apply This Skill

1. **Starting a new project**: "How should I organize my RNA-seq analysis project?"
2. **Preparing for publication**: "What do I need for a Data Availability statement?"
3. **Documenting environments**: "How do I create a reproducible Conda environment?"
4. **Data deposition**: "Where should I deposit my proteomics data?"
5. **Writing Methods sections**: "What computational details should I include?"
6. **Archiving code**: "How do I get a DOI for my GitHub repository?"

</usage_patterns>

<common_pitfalls>
## Reproducibility Errors to Avoid

| Error | Severity | Impact on Score | Fix |
|-------|----------|-----------------|-----|
| Unpinned dependencies | 🔴 Critical | -20% | Pin all versions with `==` |
| Missing random seeds | 🔴 Critical | -15% | Set seed at script start |
| Hardcoded paths | 🟡 Major | -10% | Use relative paths or config |
| Undocumented preprocessing | 🟡 Major | -10% | Script all transformations |
| Missing software versions | 🟡 Major | -10% | Include in environment file |
| Broken data links | 🔴 Critical | -15% | Use DOIs and archived copies |
| No environment file | 🔴 Critical | -25% | Create environment.yml |
| Uncommitted changes | 🟡 Major | -10% | Commit before running |
| Missing intermediate files | 🟢 Minor | -5% | Document regeneration |
| No license | 🟢 Minor | -5% | Add LICENSE file |

**Detailed Checklist (20 items):**

| Category | Item | Points |
|----------|------|--------|
| **Environment (25%)** | requirements.txt or environment.yml present | 5 |
| | All dependencies pinned to exact versions | 10 |
| | Python/R version specified | 5 |
| | System requirements documented | 5 |
| **Data (25%)** | Raw data archived with DOI | 10 |
| | Data processing scripts included | 5 |
| | Random seeds set for all stochastic operations | 5 |
| | Checksums for input files | 5 |
| **Code (25%)** | All analysis code in repository | 10 |
| | No hardcoded absolute paths | 5 |
| | README with usage instructions | 5 |
| | LICENSE file present | 5 |
| **Documentation (25%)** | Methods sufficient for replication | 10 |
| | Parameter values documented | 5 |
| | Expected outputs described | 5 |
| | Execution order documented | 5 |

**Minimum for publication:** ≥75 points (15/20 items)

</common_pitfalls>

<data_availability_templates>
## Data Availability Statement Templates

### Template 1: Public Repository
```
The sequencing data generated in this study have been deposited in
the NCBI Gene Expression Omnibus (GEO) under accession number GSExxxxx.
Processed data and analysis code are available at
https://github.com/username/project (DOI: 10.5281/zenodo.xxxxxxx).
```

### Template 2: Controlled Access
```
Raw sequencing data are available through the European Genome-phenome
Archive (EGA) under accession EGAS00001xxxxxx. Access requires approval
from the Data Access Committee due to patient privacy considerations.
Processed, de-identified summary statistics are freely available at
[repository URL].
```

### Template 3: Existing Data
```
This study uses publicly available data from [Source] (accession: xxx).
No new data were generated. Analysis code is available at
https://github.com/username/project.
```

</data_availability_templates>

<cross_references>
## Related Skills

**Skill Selection:** See `SKILL_ROUTER.md` for decision trees when multiple skills may apply.

- **scientific-writing**: Data Availability statement text, Methods section writing
- **code-documentation**: README templates, docstring standards
- **venue-templates**: Journal-specific deposition requirements
- **statistical-analysis**: Reproducible analysis code patterns
</cross_references>
