# Skills Folder Analysis Report

> Generated: 2025-12-26
> **Updated: 2025-12-28** - Applied consolidation recommendations
> Scope: Review of `/skills/` folder for redundancies and consolidation opportunities

## Executive Summary

The skills folder contains **21 top-level skills** (plus 4 document sub-skills) that cover the complete research lifecycle from literature discovery through publication and dissemination.

**Changes Applied (2025-12-28):**
- ✅ Merged `scientific-critical-thinking` into `peer-review` (reduced skill count by 1)
- ✅ Added `SKILL.md` to `document-skills/` container
- ✅ Updated cross-references in affected skills
- ✅ Updated SKILL_ROUTER.md routing logic

---

## Skills Inventory

### By Category

| Category | Skills | Count |
|----------|--------|-------|
| **Writing** | `scientific-writing`, `literature-review`, `hypothesis-generation` | 3 |
| **Presentations** | `scientific-slides`, `latex-posters`, `pptx-posters` | 3 |
| **Research** | `research-lookup`, `citation-management`, `peer-review`, `scholar-evaluation` | 4 |
| **Visuals** | `scientific-schematics`, `generate-image`, `plotting-libraries`, `visual-design` | 4 |
| **Documents** | `markitdown`, `venue-templates`, `document-skills` (container with 4 sub-skills) | 3 |
| **Analysis** | `statistical-analysis`, `code-documentation`, `reproducible-research` | 3 |
| **Conversion** | `paper-2-web` | 1 |

**Total: 21 top-level skills + 4 document sub-skills (docx, pdf, pptx, xlsx)**

### Structure Summary

| Skill | SKILL.md | scripts/ | references/ | assets/ | Complexity |
|-------|:--------:|:--------:|:-----------:|:-------:|------------|
| scientific-writing | ✓ | | ✓ | | Medium |
| literature-review | ✓ | ✓ | ✓ | ✓ | Medium-High |
| hypothesis-generation | ✓ | ✓ | ✓ | ✓ | High |
| scientific-slides | ✓ | ✓ | ✓ | ✓ | Very High |
| latex-posters | ✓ | ✓ | ✓ | ✓ | High |
| pptx-posters | ✓ | | ✓ | ✓ | Medium |
| research-lookup | ✓ | ✓ | ✓ | | Medium-High |
| citation-management | ✓ | ✓ | ✓ | ✓ | Medium |
| peer-review | ✓ | | ✓ | | **Very High** (includes critical analysis) |
| scholar-evaluation | ✓ | ✓ | | | Medium-High |
| scientific-schematics | ✓ | ✓ | ✓ | | Very High |
| generate-image | ✓ | ✓ | | | Medium |
| plotting-libraries | ✓ | | ✓ | | Medium |
| visual-design | ✓ | | ✓ | | Medium |
| markitdown | ✓ | ✓ | ✓ | ✓ | Medium-High |
| paper-2-web | ✓ | | ✓ | | Very High |
| venue-templates | ✓ | ✓ | ✓ | ✓ | High |
| statistical-analysis | ✓ | | ✓ | | High |
| code-documentation | ✓ | | ✓ | | Medium |
| reproducible-research | ✓ | | ✓ | | Medium-High |
| document-skills | ✓ | | | | Container |

---

## Completed Consolidations

### 1. Scientific Critical Thinking → Peer Review ✅

**Date:** 2025-12-28

The `scientific-critical-thinking` skill has been merged into `peer-review` as a "Critical Analysis Framework" section.

**What was merged:**
- Evidence Quality Assessment (GRADE framework)
- Logical Fallacy Identification
- Claim Evaluation process
- Bias Detection taxonomy

**Reference files moved to peer-review/references/:**
- `common_biases.md`
- `evidence_hierarchy.md`
- `experimental_design.md`
- `logical_fallacies.md`
- `scientific_method.md`
- `statistical_pitfalls.md`

**Files updated:**
- `skills/peer-review/SKILL.md` - Added Critical Analysis Framework section
- `skills/statistical-analysis/SKILL.md` - Updated cross-references
- `skills/scholar-evaluation/SKILL.md` - Updated cross-references
- `skills/SKILL_ROUTER.md` - Updated routing logic
- `skills/SKILL_TESTS.md` - Updated test scenarios

| Before | After |
|--------|-------|
| 3 review skills | 2 review skills |
| `peer-review` | `peer-review` (enhanced with critical analysis) |
| `scholar-evaluation` | `scholar-evaluation` (unchanged) |
| `scientific-critical-thinking` | (absorbed into peer-review) |

### 2. Document Skills Container ✅

**Date:** 2025-12-28

Added `SKILL.md` to `document-skills/` explaining the container structure and routing logic.

**Container structure:**
```
document-skills/
├── SKILL.md           ← NEW: Container routing logic
├── docx/SKILL.md      → Word documents
├── pdf/SKILL.md       → PDF documents
├── pptx/SKILL.md      → PowerPoint presentations
└── xlsx/SKILL.md      → Excel spreadsheets
```

---

## Non-Redundant Skill Groups

### Visual Generation (3 skills) ✓ Keep Separate

| Skill | Use Case | Backend | Output |
|-------|----------|---------|--------|
| `scientific-schematics` | Technical diagrams (flowcharts, circuits, pathways) | Nano Banana Pro + Gemini review | Publication-quality diagrams |
| `generate-image` | Artistic/photorealistic imagery | FLUX.2 Pro, Gemini 3 Pro | PNG images |
| `plotting-libraries` | Data-driven charts/plots | matplotlib/seaborn | Python code |

**Verdict**: Keep separate - clear purpose boundaries, different tooling, different quality control mechanisms.

### Poster Skills (2 skills) ✓ Keep Separate

| Skill | Format | Use Case |
|-------|--------|----------|
| `latex-posters` | PDF via LaTeX | High control, academic standard, precise layout |
| `pptx-posters` | PowerPoint | Easier editing, collaboration, iterative design |

**Verdict**: Keep separate - different output formats serve different user workflows and preferences.

### Design Guidance (Layered Architecture) ✓ Keep Separate

| Skill | Level | Scope |
|-------|-------|-------|
| `visual-design` | Philosophy/principles | All visual outputs |
| `plotting-libraries` | Decision framework | Data visualization |
| `scientific-slides` | Implementation | Presentations |
| `latex-posters` | Implementation | Posters |

**Verdict**: Intentional layering - `visual-design` provides principles, others implement for specific formats.

---

## Dependency Patterns

### High-Dependency Skills

These skills are required or referenced by many others:

#### 1. `scientific-schematics` (Most Referenced)

Marked as **MANDATORY** in:
- `hypothesis-generation` — 1-2 figures minimum required
- `latex-posters` — 2-3 figures required (40-50% of poster area)
- `scientific-slides` — Essential for technical diagrams
- `paper-2-web` — Methodology/pipeline diagrams

Also referenced by: `markitdown`, `peer-review`, `venue-templates`

#### 2. `research-lookup` (Critical Infrastructure)

Required by:
- `scientific-slides` — 8-15 papers for citations
- `literature-review` — Primary tool for finding sources
- `hypothesis-generation` — Literature synthesis phase
- `scientific-writing` — Background/context research

#### 3. `visual-design` (Philosophy Layer)

Referenced by:
- `plotting-libraries` — Styling integration
- `scientific-slides` — Design philosophy
- `latex-posters` — Design principles
- `paper-2-web` — Visual requirements

### Workflow Chains

```
Research Initiation:
research-lookup → literature-review → hypothesis-generation → scientific-writing

Presentation:
scientific-writing → scientific-slides → [scientific-schematics + research-lookup]

Publication:
venue-templates → scientific-writing → peer-review → citation-management

Dissemination:
scientific-writing → paper-2-web → [web/video/poster outputs]
```

---

## Structural Improvements Made

### 1. ✅ Document Skills Container (Fixed)

`document-skills/` now has a `SKILL.md` explaining:
- Container structure and routing logic
- Decision tree for format selection
- Links to sub-skills

### 2. Complexity Imbalance (Noted)

Some skills are significantly more complex than others:

| High Complexity | Reason |
|-----------------|--------|
| `scientific-slides` | Covers 3 distinct workflows (Nano Banana Pro, PowerPoint, LaTeX Beamer) |
| `latex-posters` | Multiple LaTeX packages, extensive templates, quality control workflow |
| `scientific-schematics` | Two-stage AI pipeline with quality thresholds |
| `paper-2-web` | Autonomous multi-output pipeline |
| `peer-review` | Now includes 7-stage review + critical analysis framework |

**Consideration**: `scientific-slides` could potentially be split into 3 sub-skills by workflow type.

### 3. MANDATORY Requirements

Figure requirements are documented in SKILL_ROUTER.md and individual skills:

- `hypothesis-generation`: "1-2 figures minimum"
- `latex-posters`: "2-3 figures, 40-50% of area"
- `scientific-slides`: "figures/diagrams essential"

---

## Current State Summary

### Skill Count

| State | Count |
|-------|-------|
| Top-level skills | 21 |
| Document sub-skills | 4 |
| **Total skills** | **25** |

### Files in skills/ (Non-Skill)

| File | Purpose |
|------|---------|
| `SKILL_ROUTER.md` | Decision trees for skill selection |
| `SKILL_TEMPLATE.md` | Template for creating new skills |
| `SKILL_TESTS.md` | Validation test scenarios |
| `QUANTIFICATION_THRESHOLDS.md` | Shared quality thresholds |

---

## Appendix: Skill Purposes

| Skill | One-Line Purpose |
|-------|------------------|
| `scientific-writing` | Write research papers (IMRaD structure) |
| `literature-review` | Conduct systematic literature reviews |
| `hypothesis-generation` | Develop testable scientific hypotheses |
| `scientific-slides` | Create research presentations |
| `latex-posters` | Create LaTeX-based academic posters |
| `pptx-posters` | Create PowerPoint-based posters |
| `research-lookup` | Real-time research lookup via Perplexity |
| `citation-management` | Organize and format citations |
| `peer-review` | Systematic manuscript/grant evaluation + critical analysis |
| `scholar-evaluation` | 8-dimension scholarly work assessment |
| `scientific-schematics` | AI-powered technical diagram generation |
| `generate-image` | General-purpose AI image generation |
| `plotting-libraries` | Python plotting library decision framework |
| `visual-design` | Design philosophy for scientific visuals |
| `markitdown` | Convert documents to Markdown |
| `paper-2-web` | Transform papers to web/video/poster |
| `venue-templates` | Publication venue templates and formatting |
| `statistical-analysis` | Statistical method selection and reporting |
| `code-documentation` | Code documentation standards and templates |
| `reproducible-research` | FAIR principles and reproducibility workflows |
| `document-skills` | Document format processing (container) |

---

*Analysis conducted by Claude Opus 4.5*
*Updated: 2025-12-28 - Consolidation recommendations applied*
