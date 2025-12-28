# Quantification Thresholds Reference

> Shared numeric thresholds replacing qualitative terms across all skills
> Version: 1.0.0
> Last Updated: 2025-12-27

---

## Purpose

This document provides **quantified definitions** for qualitative terms used across skills. When a skill uses terms like "comprehensive," "sufficient," or "appropriate," refer to these tables for specific numeric thresholds.

**Usage:** Skills should reference this document rather than defining thresholds locally, ensuring consistency across the library.

---

## 1. Literature Coverage Thresholds

### By Document Type

| Document Type | Min Papers | Min Databases | Time Range | Citation Density |
|---------------|------------|---------------|------------|------------------|
| Original Research | 30 | 2 | 5-10 years | 2-4 per paragraph |
| Review Article | 100 | 4 | 10+ years | 4-6 per paragraph |
| Systematic Review | 50+ | 3 | Per PRISMA | All relevant |
| Meta-Analysis | 10+ studies | 3 | 10+ years | All included |
| Methods Paper | 20 | 2 | 5 years | Key comparisons |
| Commentary/Letter | 10-15 | 1 | 2-3 years | Selective |

### Database Coverage Requirements

| Field | Required Databases | Recommended Additional |
|-------|-------------------|----------------------|
| Biomedical | PubMed, Embase | Cochrane, Web of Science |
| Life Sciences | PubMed, Web of Science | Scopus, Google Scholar |
| Bioinformatics | PubMed, arXiv | bioRxiv, IEEE Xplore |
| Clinical | PubMed, Cochrane | CINAHL, PsycINFO |
| Chemistry | SciFinder, Web of Science | Reaxys, PubChem |

### Search Saturation Criteria

A literature search is **"comprehensive"** when:
- [ ] ≥3 databases searched with documented strategy
- [ ] Forward/backward citation search performed
- [ ] Gray literature checked (preprints, theses, conference proceedings)
- [ ] ≤5% new relevant papers found in final search iteration
- [ ] Date range explicitly specified and justified

---

## 2. Visual Quality Thresholds

### Resolution by Context

| Context | Min Resolution | Recommended | Max File Size |
|---------|---------------|-------------|---------------|
| Print (journal) | 300 DPI | 600 DPI | Per journal spec |
| Print (poster) | 150 DPI | 300 DPI | — |
| Screen/web | 72 DPI | 150 DPI | 2 MB |
| Presentation | 150 DPI | 300 DPI | 5 MB per slide |

### Accessibility Standards (WCAG 2.1)

| Element | Minimum | Target | Excellent |
|---------|---------|--------|-----------|
| Text contrast ratio | 4.5:1 | 7:1 | 10:1 |
| Large text contrast | 3:1 | 4.5:1 | 7:1 |
| Non-text contrast | 3:1 | 4.5:1 | — |
| Focus indicator | 3:1 | 4.5:1 | — |

### Font Size by Context

| Context | Minimum | Recommended | Maximum |
|---------|---------|-------------|---------|
| Figure labels | 8pt | 10-12pt | 14pt |
| Figure captions | 9pt | 10pt | 12pt |
| Poster body | 24pt | 28-32pt | 48pt |
| Poster title | 72pt | 96pt | 120pt |
| Slide body | 18pt | 24pt | 32pt |
| Slide title | 28pt | 36pt | 44pt |

### Color Usage

| Requirement | Threshold | Tool |
|-------------|-----------|------|
| Colorblind-safe palette | ≤8 distinct colors | Colorbrewer, viridis |
| Distinguishable colors | ΔE ≥ 20 (CIE2000) | Color contrast checker |
| Sequential palettes | Monotonic lightness | viridis, plasma, inferno |
| Diverging palettes | Symmetric around midpoint | RdBu, PiYG, coolwarm |

---

## 3. Replication & Sample Size Thresholds

### Biological Replication

| Experiment Type | Min Biological Reps | Min Technical Reps | Notes |
|-----------------|--------------------|--------------------|-------|
| Cell culture (in vitro) | 3 | 3 | Independent passages |
| Primary cells | 3-5 | 2-3 | Different donors |
| Mouse studies | 6-10 | — | Per power calculation |
| Rat studies | 6-8 | — | Per power calculation |
| Zebrafish | 10-20 | — | High variability |
| Drosophila | 30-50 | — | Population genetics |
| Human samples | Per power calc | 1-2 | IRB-approved |
| Clinical trial | Per power calc | 1 | Pre-registered |

### Statistical Power Requirements

| Study Type | Minimum Power | Effect Size Basis |
|------------|---------------|-------------------|
| Confirmatory | 80% | Pre-specified |
| Exploratory | 80% | Pilot data or literature |
| Replication | 95% | Original study |
| Non-inferiority | 80% | Clinically meaningful difference |

### Sample Size Red Flags

| Indicator | Concern Level | Action |
|-----------|---------------|--------|
| n < 3 per group | 🔴 Critical | Cannot perform statistics |
| n = 3-5 per group | 🟡 Major | Justify or increase |
| No power calculation | 🟡 Major | Request justification |
| Post-hoc power | 🔴 Critical | Invalid—reject |

---

## 4. Issue Severity Classification

### Manuscript Review Severity

| Severity | Definition | Examples | Recommended Action |
|----------|------------|----------|-------------------|
| **Critical** | Invalidates conclusions | Data fabrication, fundamental flaw in design, wrong statistical test invalidating results, ethical violation | Reject |
| **Major** | Significantly affects interpretation | Missing essential control, underpowered study, incomplete methods, major statistical error | Major revision required |
| **Minor** | Affects clarity/completeness | Unclear writing, formatting issues, minor missing details, typos | Minor revision |
| **Optional** | Suggestions for improvement | Style preferences, additional analyses that would strengthen | Author discretion |

### Quantified Severity Criteria

**Critical Issues (any one = reject consideration):**
- [ ] Data cannot be independently verified
- [ ] Statistical test fundamentally inappropriate (e.g., t-test on categorical data)
- [ ] Missing data exceeds 30% with no sensitivity analysis
- [ ] Conclusions contradicted by presented data
- [ ] Ethical approval missing for human/animal research

**Major Issues (≥2 = major revision):**
- [ ] Key control missing from ≥50% of experiments
- [ ] Sample size <80% of power calculation requirement
- [ ] Methods insufficient for replication (fails reproducibility checklist)
- [ ] Statistical assumptions violated without robust alternative
- [ ] Results section missing ≥20% of stated outcomes

**Minor Issues (address in revision):**
- [ ] Grammar/spelling errors (>5 per page)
- [ ] Figure quality below journal standards
- [ ] Reference formatting inconsistent
- [ ] Acronyms undefined on first use
- [ ] Minor gaps in Methods (recoverable from supplementary)

---

## 5. Documentation Completeness Thresholds

### Code Documentation

| Level | Coverage | Description |
|-------|----------|-------------|
| Minimal | — | One-line script description |
| Basic | 50% | Public functions documented |
| Standard | 80% | All public API documented |
| Comprehensive | 95% | All functions + examples |

### Methods Documentation

| Level | Criteria |
|-------|----------|
| **Insufficient** | Another researcher could NOT replicate |
| **Minimal** | Could replicate with author correspondence |
| **Adequate** | Could replicate from methods + supplementary |
| **Complete** | Could replicate from methods alone |
| **Exemplary** | Includes step-by-step protocol with troubleshooting |

### Reproducibility Checklist Pass Rates

| Rating | Checklist Items Passed |
|--------|----------------------|
| Not reproducible | <50% |
| Partially reproducible | 50-74% |
| Reproducible | 75-89% |
| Highly reproducible | 90-100% |

---

## 6. Time-Based Thresholds

### Literature Recency

| Context | "Recent" Means | "Current" Means |
|---------|---------------|-----------------|
| Fast-moving fields (ML, genomics) | ≤2 years | ≤6 months |
| Standard biomedical | ≤5 years | ≤1 year |
| Clinical guidelines | ≤3 years | ≤1 year |
| Historical context | ≤10 years | ≤5 years |
| Foundational work | Any age if seminal | N/A |

### Review Turnaround

| Review Type | Target Time | Maximum |
|-------------|-------------|---------|
| Initial read | 1-2 hours | 4 hours |
| Detailed review | 4-8 hours | 16 hours |
| Revision check | 1-2 hours | 4 hours |
| Statistical review | 2-4 hours | 8 hours |

---

## 7. Quality Scoring Rubrics

### Figure Quality Rubric (10-point scale)

| Criterion | 0 | 1 | 2 |
|-----------|---|---|---|
| Resolution | Pixelated/blurry | Acceptable | Crisp at 100% zoom |
| Labels | Missing/unreadable | Present but issues | Clear, consistent |
| Color | Inaccessible | Mostly accessible | Fully colorblind-safe |
| Message | Unclear | Understandable | Immediately obvious |
| Professional | Amateur appearance | Acceptable | Publication-ready |

**Scoring:** 0-4 = Needs revision, 5-7 = Acceptable, 8-10 = Excellent

### Schematic Quality Rubric (10-point scale)

| Criterion | 0 | 1 | 2 |
|-----------|---|---|---|
| Accuracy | Factual errors | Minor issues | Scientifically accurate |
| Clarity | Confusing layout | Understandable | Intuitive flow |
| Completeness | Missing key elements | Most elements | All elements present |
| Accessibility | Fails contrast | Partial compliance | Full WCAG compliance |
| Style | Inconsistent | Mostly consistent | Cohesive design |

### Writing Quality Rubric (10-point scale)

| Criterion | 0 | 1 | 2 |
|-----------|---|---|---|
| Clarity | Ambiguous/confusing | Mostly clear | Unambiguous |
| Concision | Verbose/redundant | Appropriate length | Precise |
| Flow | Disjointed | Reasonable transitions | Smooth narrative |
| Grammar | Many errors (>5/page) | Few errors (<2/page) | Error-free |
| Jargon | Unexplained terms | Most defined | Appropriate for audience |

---

## 8. Iteration & Stopping Criteria

### When to Stop Iterating

| Process | Stop When |
|---------|-----------|
| Literature search | <5% new relevant papers in final iteration |
| Figure revision | Quality score ≥8/10 OR 3 iterations completed |
| Schematic generation | Quality score ≥8/10 OR 5 iterations completed |
| Writing revision | Passes all checklist items OR editor satisfied |
| Statistical analysis | Pre-registered analysis complete |

### Maximum Iteration Limits

| Task | Soft Limit | Hard Limit |
|------|------------|------------|
| AI image generation | 3 iterations | 5 iterations |
| Figure revision | 2 rounds | 4 rounds |
| Manuscript revision | 2 major | 3 major |
| Peer review response | 1 round | 2 rounds |

---

## Cross-Reference Guide

When a skill uses these terms, apply the corresponding threshold:

| Qualitative Term | See Section | Key Metric |
|------------------|-------------|------------|
| "comprehensive literature review" | §1 | ≥50 papers, ≥3 databases |
| "sufficient replication" | §3 | ≥3 bio reps, ≥80% power |
| "appropriate statistical test" | §4, statistical-analysis | Decision matrix criteria |
| "high quality figure" | §2, §7 | ≥8/10 rubric score |
| "clear message" | §7 | 5-second test pass |
| "adequate documentation" | §5 | 80% coverage |
| "recent literature" | §6 | ≤5 years (field-dependent) |
| "major issue" | §4 | Affects interpretation |
| "minor issue" | §4 | Affects clarity only |

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | 2025-12-27 | Initial release with 8 threshold categories |

---

*This document provides the quantitative backbone for deterministic skill behavior. All skills should reference these thresholds rather than using vague qualitative terms.*
