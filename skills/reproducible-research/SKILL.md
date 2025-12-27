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

<environment_documentation>
## Environment Specification Hierarchy

### Level 1: Minimal (requirements.txt)
```
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
```
- Quick start for collaborators
- Often insufficient for exact reproducibility

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
- Pinned versions
- Cross-platform (mostly)
- Standard for bioinformatics

### Level 3: Complete (Docker/Singularity)
```dockerfile
FROM continuumio/miniconda3:23.5.2-0
COPY environment.yml .
RUN conda env create -f environment.yml
```
- Fully reproducible
- Required for HPC/cloud
- Publication-grade reproducibility

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

1. **Unpinned dependencies**: `pip install numpy` vs `numpy==1.24.3`
2. **Missing random seeds**: Results vary between runs
3. **Hardcoded paths**: `/Users/myname/data/` breaks on other systems
4. **Undocumented preprocessing**: Filtering steps not recorded
5. **Missing software versions**: "We used Python" (which version?)
6. **Broken data links**: URLs that expire or change
7. **No environment file**: "It works on my machine"
8. **Uncommitted changes**: Analysis ran on modified code
9. **Missing intermediate files**: Can't reproduce from raw data
10. **No license**: Others can't legally reuse code

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
- **scientific-writing**: Data Availability statement text, Methods section writing
- **code-documentation**: README templates, docstring standards
- **venue-templates**: Journal-specific deposition requirements
- **statistical-analysis**: Reproducible analysis code patterns
</cross_references>
