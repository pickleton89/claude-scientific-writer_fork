# Research Paper Summarizer Integration Plan

> Implementation plan for integrating the research-paper-summarizer skill into the skill library
> Created: 2025-12-30
> Status: **IN PROGRESS** (Phases 0-3 complete, Phase 4 pending)

## Overview

This document outlines the step-by-step implementation plan to fully integrate the `research-paper-summarizer` skill with the claude-scientific-writer skill library. The plan addresses brand alignment, router integration, cross-references, and full XML structure compliance.

### Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Brand System | **Oligon brand** | Project-wide consistency |
| BRAND_COLORS.md | **Delete** | Contradicts actual implementation |
| PDF Generator | **Keep separate** | Specialized figure extraction justifies it |
| Skill Structure | **Full XML restructure** | Complete compliance with skill standards |

### Scope Summary

- **Files to modify:** 8-12 files
- **Files to delete:** 1 file
- **Files to create:** 2-3 files
- **Estimated effort:** 4-6 hours

---

## Phase 0: Quick Wins

### Task 0.1: Add to SKILL_ROUTER.md

**File:** `skills/SKILL_ROUTER.md`

**Action:** Add new decision tree section for paper summarization

**Insert after "Decision Tree: Research Workflow" section:**

```markdown
---

## Decision Tree: Paper Summarization

```
User wants to summarize or analyze a research paper
│
├─ Is it a SINGLE paper needing comprehensive summary?
│  │
│  └─ YES → research-paper-summarizer
│            │
│            ├─ What article type?
│            │  ├─ General research → general_researcher_summarizer
│            │  ├─ Review/meta-analysis → review_article_summarizer
│            │  ├─ Computational/bioinformatics → compbio_bioinformatic_summarizer
│            │  └─ Cell & molecular biology → cellmolbio_summarizer
│            │
│            └─ What output format?
│               ├─ Markdown only → _summary.md
│               ├─ Interactive HTML → html-report
│               ├─ PDF summary → yaml-for-pdf → generate_summary_pdf.py
│               └─ SVG infographic → svg-infographic
│
├─ Is it a COLLECTION of papers for systematic review?
│  │
│  └─ YES → literature-review
│            └─ Then optionally → research-paper-summarizer (for each key paper)
│
├─ Is it a quick lookup for citation/abstract only?
│  │
│  └─ YES → research-lookup
│
└─ Is it evaluating paper QUALITY (not summarizing)?
   │
   └─ YES → peer-review
```

### Paper Summarization: Quick Selection

| Task | Primary Skill | Output |
|------|---------------|--------|
| Deep analysis of single paper | research-paper-summarizer | Markdown + optional visual |
| Quick abstract/citation lookup | research-lookup | Text response |
| Systematic review of many papers | literature-review | Review document |
| Evaluate paper methodology | peer-review | Evaluation report |
| Summarize paper for presentation | research-paper-summarizer → scientific-slides | Slides |
```

**Also update "Quick Reference Matrix" table to include:**

```markdown
| **Paper Analysis** | | |
| Summarize single paper | research-paper-summarizer | + oligon-brand (for visual outputs) |
| Quick paper lookup | research-lookup | — |
| Systematic literature review | literature-review | + citation-management |
```

**Verification:**
- [ ] Decision tree added after Research Workflow section
- [ ] Quick Selection table added
- [ ] Quick Reference Matrix updated
- [ ] Version number incremented in frontmatter

---

### Task 0.2: Update Related Skills' Cross-References

**File:** `skills/literature-review/SKILL.md`

**Action:** Add reference to research-paper-summarizer in cross-references section

**Add to cross-references table:**

```markdown
| `research-paper-summarizer` | Summarize individual papers found during review; use for key papers requiring deep analysis |
```

**File:** `skills/peer-review/SKILL.md`

**Action:** Clarify distinction from summarization in when_to_use or cross-references

**Add clarification:**

```markdown
| `research-paper-summarizer` | Summarization (extracting content) vs evaluation (assessing quality); use peer-review for critical assessment |
```

**File:** `skills/research-lookup/SKILL.md`

**Action:** Add reference for deeper analysis

**Add:**

```markdown
| `research-paper-summarizer` | For comprehensive paper analysis beyond quick lookup |
```

**Verification:**
- [ ] literature-review cross-reference added
- [ ] peer-review cross-reference added
- [ ] research-lookup cross-reference added

---

## Phase 1: Brand Alignment

### Task 1.1: Delete Inconsistent Brand File

**File to delete:** `skills/research-paper-summarizer/brand/BRAND_COLORS.md`

**Action:** Remove the file (1022 lines of "Scientific Minimalist" brand that contradicts the actual implementation)

```bash
rm skills/research-paper-summarizer/brand/BRAND_COLORS.md
rmdir skills/research-paper-summarizer/brand/  # Remove empty directory
```

**Verification:**
- [ ] BRAND_COLORS.md deleted
- [ ] brand/ directory removed
- [ ] No broken references in SKILL.md

---

### Task 1.2: Update SKILL.md Brand References

**File:** `skills/research-paper-summarizer/SKILL.md`

**Action:** Update brand color section to reference oligon-brand skill

**Find and replace the "Brand Visual Identity" section (around line 299-310):**

**Current:**
```markdown
### Brand Visual Identity

All visual outputs follow Scientific Brand Visual Identity v4.0:

```
Primary Highlight:  #2DB2E8 (Brand Blue) - Use sparingly for PRIMARY finding only
Contrast/Alert:     #E8622D (Orange) - Opposing effects, warnings, limitations
...
```
```

**Replace with:**
```markdown
### Brand Visual Identity

All visual outputs follow the **Oligon Scientific Brand** (see `oligon-brand` skill).

**Quick Reference:**

| Role | Color | HEX |
|------|-------|-----|
| Primary Highlight | Brand Blue | `#2DB2E8` |
| Contrast/Alert | Orange | `#E8622D` |
| Primary Data | Dark Gray | `#222222` |
| Secondary Data | Medium Gray | `#666666` |
| Tertiary Data | Muted Gray | `#999999` |
| Background | White | `#FFFFFF` |

For complete brand specification, color cycles, and accessibility guidelines, see:
- `skills/oligon-brand/SKILL.md`
- `skills/oligon-brand/references/brand-colors-full.md`
```

**Verification:**
- [ ] Brand section updated to reference oligon-brand
- [ ] Quick reference table matches Oligon colors
- [ ] Link to oligon-brand skill included

---

### Task 1.3: Update PDF Generator Import Comments

**File:** `skills/research-paper-summarizer/scripts/generate_summary_pdf.py`

**Action:** Update header comments to reference Oligon brand

**Find (around line 8-9):**
```python
- Brand-compliant styling (Scientific Brand Visual Identity v4.0)
```

**Replace with:**
```python
- Brand-compliant styling (Oligon Scientific Brand v4.0)
```

**Also update the comment block around line 61:**
```python
# =============================================================================
# BRAND COLOR PALETTE (Oligon Scientific Brand v4.0)
# =============================================================================
```

**Verification:**
- [ ] Header comment updated
- [ ] Color palette comment updated
- [ ] No functional code changes needed (colors already correct)

---

## Phase 2: Full XML Structure Compliance

### Task 2.1: Update Frontmatter

**File:** `skills/research-paper-summarizer/SKILL.md`

**Action:** Add missing frontmatter fields

**Current:**
```yaml
---
name: research-paper-summarizer
description: Summarize scientific research papers into multiple output formats...
---
```

**Replace with:**
```yaml
---
name: research-paper-summarizer
version: 1.0.0
description: "Summarize scientific research papers into multiple output formats including structured markdown, interactive HTML presentations, brand-compliant PDF summaries, and visual SVG infographics. Use when the user provides a research paper (PDF, text, DOI, or URL) and wants it summarized, analyzed, or visualized."
when_to_use: "When user uploads a research paper (PDF, text, abstract) or provides a DOI/URL and wants comprehensive summarization with optional visual outputs"
allowed-tools: [Read, Write, Edit, Bash, Glob, Grep, AskUserQuestion, Task]
brand-reference: "../oligon-brand/"
---
```

**Verification:**
- [x] version field added
- [x] when_to_use field added
- [x] allowed-tools field added
- [x] brand-reference field added

---

### Task 2.2: Wrap Existing Content in XML Sections

**File:** `skills/research-paper-summarizer/SKILL.md`

**Action:** Add XML section wrappers around existing content. The content largely stays the same, just wrapped.

**Structure to implement:**

```markdown
# Research Paper Summarizer

<overview>
[Move "Create rigorous, quantitative summaries..." paragraph here]
[Move Quick Start section here]
</overview>

<when_to_use>
[Keep existing "When to Use This Skill" content]
</when_to_use>

<decision_framework>
[Move "Processing Modes" section here]
[Move "Phase 1: Article Type Selection" here]
</decision_framework>

<workflow>
[Move "Workflow Overview" section here]
[Move "Chunked Mode Workflow" section here]
[Move "Phase 2: Visual Output Format Selection" section here]
</workflow>

<output_formats>
[Move format-specific instructions here: HTML, PDF, SVG]
</output_formats>

<core_principles>
[Move "Core Principles (All Formats)" section here]
[Move "Quantitative Preservation" here]
[Move "Scientific Rigor" here]
</core_principles>

<article_types>
[Move "Article Type Details" section here]
[Move "Special Handling Notes" here]
</article_types>

<error_handling>
[Move "Troubleshooting & Error Recovery" section here]
</error_handling>

<anti_patterns>
## Common Pitfalls

### 1. Skipping Article Type Selection
**Problem:** Using general summarizer for specialized content (compbio, cellmolbio).
**Solution:** Always ask user to confirm article type before processing.

### 2. Ignoring Page Count Threshold
**Problem:** Using Standard Mode for papers >12 pages, causing context exhaustion.
**Solution:** Let the skill auto-detect; trust the 12-page threshold.

### 3. Summarizing Away Statistics
**Problem:** Replacing "p = 0.003" with "statistically significant."
**Solution:** Preserve ALL numerical data exactly as reported.

### 4. Missing Figure References
**Problem:** Summary mentions findings without linking to source figures.
**Solution:** Include figure_ref for every finding with supporting data.

### 5. Incomplete Chunked Processing
**Problem:** Stopping when one subagent fails.
**Solution:** Use recovery workflow; check for PENDING markers and resume.
</anti_patterns>

<success_criteria>
## Success Criteria

### Markdown Summary Quality
- [ ] All statistics preserved with exact values
- [ ] Every finding linked to figure/table reference
- [ ] Limitations section includes both author-stated and additional
- [ ] Technical terminology preserved (no paraphrasing jargon)
- [ ] Conditional language preserved (authors' hedging)

### Visual Output Quality
- [ ] Brand colors applied correctly (Oligon brand)
- [ ] White background, black text
- [ ] Statistics formatted consistently (p = X.XX, 95% CI: X–X)
- [ ] Figures extracted and matched to references (PDF output)
</success_criteria>

<cross_references>
## Related Skills

| Skill | Relationship |
|-------|-------------|
| `oligon-brand` | **Brand colors** — all visual outputs (HTML, PDF, SVG) use Oligon brand styling |
| `literature-review` | **Finding papers** — use literature-review to find papers, then summarize key ones |
| `research-lookup` | **Quick lookup** — for simple queries; use research-paper-summarizer for deep analysis |
| `peer-review` | **Evaluation vs summary** — peer-review assesses quality; this skill extracts content |
| `markdown-to-pdf` | **Alternative PDF** — for template-based documents; this skill has specialized figure extraction |
| `scientific-critical-thinking` | **Analysis framework** — provides structure for the critical analysis sections |
| `plotting-libraries` | **Figure creation** — if creating new figures from extracted data |
| `visual-design` | **Design principles** — for infographic composition guidance |
| `markitdown` | **PDF conversion** — convert paper PDFs to markdown before summarizing (optional) |

**Skill Selection:** See `SKILL_ROUTER.md` for decision trees when multiple skills may apply.

**Typical Workflow:**
1. Receive paper PDF from user
2. `research-paper-summarizer` → generate markdown summary
3. (Optional) `oligon-brand` colors applied → visual output
4. (Optional) `scientific-critical-thinking` → deeper analysis
</cross_references>

<references>
[Keep existing "Supporting Files" section, reorganized]
</references>
```

**Verification:**
- [x] `<overview>` section added
- [x] `<when_to_use>` section wrapped
- [x] `<decision_framework>` section added
- [x] `<workflow>` section added
- [x] `<output_formats>` section added
- [x] `<core_principles>` section added
- [x] `<article_types>` section added
- [x] `<error_handling>` section wrapped
- [x] `<anti_patterns>` section added (new content)
- [x] `<success_criteria>` section added (new content)
- [x] `<cross_references>` section added (new content)
- [x] `<references>` section wrapped

---

## Phase 3: Documentation

### Task 3.1: Create Subagent Architecture Document

**File to create:** `skills/research-paper-summarizer/references/subagent_architecture.md`

**Content:**

```markdown
# Subagent Pipeline Architecture

> Technical reference for the chunked processing subagent system
> Version: 1.0.0

## Overview

The research-paper-summarizer uses a multi-agent pipeline for processing papers >12 pages. This document describes the architecture, coordination protocol, and error handling.

## Processing Mode Selection

```
Paper received
│
├─ Count pages
│
├─ pages ≤ 12?
│  └─ YES → Standard Mode (single-pass, legacy prompts)
│
└─ pages > 12?
   └─ YES → Chunked Mode (subagent pipeline)
```

## Chunked Mode Architecture

### Pipeline Stages

```
┌─────────────────────────────────────────────────────────────┐
│                    SETUP PHASE                               │
├─────────────────────────────────────────────────────────────┤
│  1. Section Detection                                        │
│     └─ Input: First 5-10 pages                              │
│     └─ Output: JSON mapping {section → page_range}          │
│                                                              │
│  2. Skeleton Creation                                        │
│     └─ Input: summary-skeleton.md + skeleton-config.md      │
│     └─ Output: {filename}_summary.md with PENDING markers   │
└─────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────┐
│                   DISPATCH PHASE                             │
├─────────────────────────────────────────────────────────────┤
│  For each subagent in sequence:                             │
│                                                              │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐     │
│  │  overview   │ →  │   methods   │ →  │   results   │     │
│  │   agent     │    │    agent    │    │    agent    │     │
│  └─────────────┘    └─────────────┘    └─────────────┘     │
│         │                  │                  │              │
│         ▼                  ▼                  ▼              │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐     │
│  │  critique   │ →  │   context   │ →  │  synthesis  │     │
│  │   agent     │    │    agent    │    │    agent    │     │
│  └─────────────┘    └─────────────┘    └─────────────┘     │
│                                                              │
│  Each agent:                                                 │
│  1. Reads assigned PDF pages                                │
│  2. Reads current summary file                              │
│  3. Writes output to summary file                           │
│  4. Updates section marker: PENDING → COMPLETE              │
└─────────────────────────────────────────────────────────────┘
```

### Section Marker Protocol

The summary file uses HTML comment markers to track progress:

```html
<!-- SECTION: executive_summary PENDING -->
[placeholder text]
<!-- END SECTION: executive_summary -->
```

After subagent writes:

```html
<!-- SECTION: executive_summary COMPLETE -->
[actual content from subagent]
<!-- END SECTION: executive_summary -->
```

### Marker States

| State | Meaning |
|-------|---------|
| `PENDING` | Section not yet processed |
| `COMPLETE` | Section successfully written |
| `SKIPPED` | Section not applicable for article type |
| `FAILED` | Subagent failed after retry |

## Article-Type-Specific Sequences

### General Research
```
overview → methods → results → critique → context → synthesis
```

### Review Article
```
overview → evidence → landscape → critique → references → synthesis
```

### Cell & Molecular Biology
```
overview → models → mechanisms → results → critique → translational → synthesis
```

### Computational Biology
```
overview → data → methods → validation → critique → reproducibility → synthesis
```

## Subagent Inventory

### Shared Agents (All Article Types)

| Agent | Location | Purpose |
|-------|----------|---------|
| overview-agent | `subagents/overview-agent.md` | Executive summary, background, hypothesis |
| critique-agent | `subagents/critique-agent.md` | Critical analysis, limitations |
| synthesis-agent | `subagents/synthesis-agent.md` | Final synthesis (reads summary only) |

### Article-Type-Specific Agents

| Article Type | Directory | Agents |
|--------------|-----------|--------|
| General | `subagents/general/` | methods, results, context |
| Review | `subagents/review/` | evidence, landscape, references |
| Cell/Mol Bio | `subagents/cellmolbio/` | models, mechanisms, results, translational |
| Comp Bio | `subagents/compbio/` | data, methods, validation, reproducibility |

## Error Handling

### Retry Logic

```
Subagent fails
│
├─ Attempt 1 failed?
│  └─ Retry with same parameters
│
├─ Attempt 2 failed?
│  └─ Mark section as FAILED
│  └─ Log error details
│  └─ Continue to next agent
│
└─ Multiple failures?
   └─ Fall back to Standard Mode
   └─ Notify user of partial completion
```

### Recovery Workflow

User can resume incomplete summarization:

```
User: Continue summarizing this paper - some sections are still PENDING

Claude: [Reads summary file]
        Found 3 sections marked PENDING: critique, context, synthesis
        Resuming from critique agent...
```

## Performance Characteristics

| Paper Length | Mode | Approximate Time |
|--------------|------|------------------|
| ≤12 pages | Standard | 2-3 minutes |
| 13-24 pages | Chunked | 5-8 minutes |
| 25-40 pages | Chunked | 8-12 minutes |
| >40 pages | Chunked | 12-15 minutes |

Visual output generation adds ~1-2 minutes per format.
```

**Verification:**
- [x] Architecture document created
- [x] Pipeline diagram included
- [x] Section marker protocol documented
- [x] Article-type sequences listed
- [x] Error handling documented

---

### Task 3.2: Update CHANGELOG.md

**File:** `CHANGELOG.md`

**Action:** Add entry for integration work

**Add under appropriate date:**

```markdown
### Research Paper Summarizer Integration

- **Router integration**: Added decision tree to `SKILL_ROUTER.md` for paper summarization tasks
- **Brand alignment**: Removed conflicting `brand/BRAND_COLORS.md` (was "Scientific Minimalist" v5.0, but implementation used Oligon v4.0)
- **Skill structure**: Full XML restructure of `SKILL.md` with all standard sections
- **Cross-references**: Added bidirectional references with `literature-review`, `peer-review`, `research-lookup`, `oligon-brand`
- **Documentation**: Created `references/subagent_architecture.md` documenting the chunked processing pipeline
- **Anti-patterns**: Added common pitfalls section
- **Success criteria**: Added quality checklist for summaries
```

**Verification:**
- [x] CHANGELOG.md updated with integration details

---

## Phase 4: Verification

### Task 4.1: Verify All File References

**Action:** Ensure all referenced files exist

**Check list:**

```bash
# Prompts (Standard Mode)
test -f skills/research-paper-summarizer/prompts/general_researcher_summarizer.md
test -f skills/research-paper-summarizer/prompts/review_article_summarizer.md
test -f skills/research-paper-summarizer/prompts/compbio_bioinformatic_summarizer.md
test -f skills/research-paper-summarizer/prompts/cellmolbio_summarizer.md

# Visual output prompts
test -f skills/research-paper-summarizer/prompts/html-report.md
test -f skills/research-paper-summarizer/prompts/yaml-for-pdf.md
test -f skills/research-paper-summarizer/prompts/svg-infographic.md

# Templates
test -f skills/research-paper-summarizer/templates/summary-skeleton.md
test -f skills/research-paper-summarizer/templates/skeleton-config.md

# Scripts
test -f skills/research-paper-summarizer/scripts/generate_summary_pdf.py

# Subagents (shared)
test -f skills/research-paper-summarizer/subagents/overview-agent.md
test -f skills/research-paper-summarizer/subagents/critique-agent.md
test -f skills/research-paper-summarizer/subagents/synthesis-agent.md
```

**Verification:**
- [ ] All prompt files exist
- [ ] All template files exist
- [ ] Script exists
- [ ] Shared subagents exist

---

### Task 4.2: Test Router Integration

**Action:** Verify SKILL_ROUTER.md is valid markdown and decision trees render correctly

**Verification:**
- [ ] No markdown syntax errors
- [ ] Decision tree ASCII renders correctly
- [ ] Tables are well-formed

---

### Task 4.3: Final Review

**Checklist:**

- [ ] `brand/BRAND_COLORS.md` deleted
- [ ] `brand/` directory removed
- [ ] SKILL.md frontmatter complete
- [ ] SKILL.md has all XML sections
- [ ] SKILL_ROUTER.md updated
- [ ] Cross-references added to 3+ related skills
- [ ] `references/subagent_architecture.md` created
- [ ] CHANGELOG.md updated
- [ ] No broken file references

---

## Implementation Order

Execute tasks in this order to minimize risk:

```
Phase 0.1: Add to SKILL_ROUTER.md
    ↓
Phase 0.2: Update related skills' cross-references
    ↓
Phase 1.1: Delete brand/BRAND_COLORS.md
    ↓
Phase 1.2: Update SKILL.md brand references
    ↓
Phase 1.3: Update PDF generator comments
    ↓
Phase 2.1: Update frontmatter
    ↓
Phase 2.2: Add XML section wrappers
    ↓
Phase 3.1: Create subagent architecture doc
    ↓
Phase 3.2: Update CHANGELOG.md
    ↓
Phase 4.1-4.3: Verification
    ↓
Git commit
```

---

## Git Commit Strategy

**Option A: Single commit (recommended for cohesive change)**

```bash
git add -A
git commit -m "feat(research-paper-summarizer): full skill library integration

- Add to SKILL_ROUTER.md with paper summarization decision tree
- Remove conflicting brand/BRAND_COLORS.md (use Oligon brand)
- Restructure SKILL.md with full XML section compliance
- Add cross-references to/from related skills
- Create references/subagent_architecture.md
- Add anti-patterns and success criteria sections

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

**Option B: Multiple commits (for granular history)**

1. `feat(skill-router): add paper summarization decision tree`
2. `fix(research-paper-summarizer): remove conflicting brand file`
3. `refactor(research-paper-summarizer): XML structure compliance`
4. `docs(research-paper-summarizer): add subagent architecture reference`

---

## Post-Implementation

After completing integration:

1. **Test the skill** by summarizing a sample paper
2. **Verify visual outputs** generate with correct Oligon colors
3. **Check router** by asking ambiguous questions to see if routing works
4. **Update version** in frontmatter if significant changes

---

*Plan created: 2025-12-30*
*Target completion: [TBD]*
