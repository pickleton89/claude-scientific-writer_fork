---
name: research-paper-summarizer
description: Summarize scientific research papers into multiple output formats including structured markdown, interactive HTML presentations, brand-compliant PDF summaries, and visual SVG infographics. Use when the user provides a research paper (PDF, text, DOI, or URL) and wants it summarized, analyzed, or visualized.
---

# Research Paper Summarizer

Create rigorous, quantitative summaries of scientific research papers in multiple output formats while preserving all numerical data, statistical details, and domain terminology.

---

## Quick Start

### Basic Usage

1. **Provide a PDF**: Upload or specify path to a research paper
2. **Select article type**: General / Review / CompBio / CellMolBio
3. **Get markdown summary**: Saved as `{filename}_summary.md`
4. **Choose visual format**: HTML report, PDF, SVG infographic, or none

### Example Session

```
User: Summarize this paper [uploads smith2024_cancer.pdf]

Claude: This paper has 24 pages. Using Chunked processing mode.
        What type of article is this?
        1. General research paper
        2. Review article
        3. Computational biology
        4. Cell & molecular biology

User: 1

Claude: [Processes paper through 6 specialized agents]
        ✓ Summary complete! Saved to: smith2024_cancer_summary.md

        Which visual format would you like?
        1. Interactive HTML report
        2. PDF summary
        3. SVG infographic
        4. None

User: 1

Claude: ✓ HTML report generated: smith2024_cancer_report.html
```

### What You Get

| Output | Format | Best For |
|--------|--------|----------|
| Markdown summary | `_summary.md` | Reading, editing, version control |
| HTML report | `_report.html` | Presentations, sharing |
| PDF summary | `_summary.pdf` | Printing, archival |
| SVG infographic | `_infographic.svg` | Posters, slides |

### Timing

- **Short papers** (≤12 pages): ~2-3 minutes
- **Long papers** (>12 pages): ~5-8 minutes (chunked processing)
- Visual output generation: ~1-2 minutes per format

---

## When to Use This Skill

Activate when the user:
- Uploads or pastes a research paper (PDF, text, abstract)
- Provides a DOI, PubMed ID, or paper URL
- Asks to "summarize this paper", "create a summary", or "analyze this research"
- Requests specific outputs: "make an infographic", "create HTML report", "generate PDF summary"

---

## Processing Modes

This skill automatically selects the best processing approach based on paper length:

| Mode | Paper Length | Approach |
|------|--------------|----------|
| **Standard Mode** | ≤12 pages | Full document analysis using legacy prompts |
| **Chunked Mode** | >12 pages | Section-by-section processing via subagent pipeline |

### Why Two Modes?

- **Standard Mode**: Efficient for shorter papers that fit comfortably in context
- **Chunked Mode**: Handles large papers without context exhaustion by processing sections incrementally and writing to a shared summary file

### Processing Mode Detection

When a PDF is provided:

1. **Read PDF and count pages**
2. **Determine mode**:
   - If pages ≤ 12: Use Standard Mode (legacy prompts)
   - If pages > 12: Use Chunked Mode (subagent pipeline)
3. **Inform user**: "This paper has X pages. Using [Standard/Chunked] processing mode."

---

## Workflow Overview

This skill follows a multi-phase workflow with mode-dependent processing:

### Phase 0: Initialization (All Papers)
1. Receive PDF from user
2. Read PDF and count pages
3. Determine processing mode (Standard ≤12 pages, Chunked >12 pages)
4. Inform user of mode selection

### Phase 1: Article Type Selection & Summary Generation

**1a. Ask user to select article type** (both modes)

**1b. Generate summary** (mode-dependent):
- **Standard Mode**: Use legacy full-document prompts
- **Chunked Mode**: Use subagent pipeline (see Chunked Mode Workflow below)

**1c. Save markdown file** to same directory as source document

### Phase 2: Visual Output Generation
- Ask user which visual format(s) they want
- Generate visuals using the saved markdown summary as the source

---

## Phase 1: Article Type Selection

**ALWAYS ask the user which type of article they are summarizing:**

| Article Type | Best For | Prompt File |
|--------------|----------|-------------|
| **General Research** | Standard primary research papers | `prompts/general_researcher_summarizer.md` |
| **Review Article** | Review papers, meta-analyses, systematic reviews | `prompts/review_article_summarizer.md` |
| **Computational/Bioinformatics** | Genomics, methods, pipelines, algorithms | `prompts/compbio_bioinformatic_summarizer.md` |
| **Cell & Molecular Biology** | Mechanistic, wet-lab, cancer biology papers | `prompts/cellmolbio_summarizer.md` |

Ask: "What type of article is this?"
1. General research paper
2. Review article / meta-analysis
3. Computational biology / bioinformatics
4. Cell & molecular biology

### Generating the Markdown Summary (Standard Mode)

For papers ≤12 pages, use the standard approach:

1. Read the appropriate summarizer prompt file
2. Apply the framework to analyze the paper comprehensively
3. Generate the markdown summary following the prompt's structure
4. **Save the summary** to the same directory as the source document
   - Filename format: `{original_filename}_summary.md`
   - Example: `smith2024_cancer_cells.pdf` → `smith2024_cancer_cells_summary.md`

---

## Chunked Mode Workflow (Papers >12 Pages)

For longer papers, use the chunked subagent pipeline to avoid context exhaustion.

### Setup Phase

1. **Detect PDF sections** using `prompts/subagents/_shared/section-detection.md`
   - Read first 5-10 pages of PDF
   - Identify page ranges for: Abstract, Introduction, Methods, Results, Discussion, References
   - Output JSON mapping of sections to pages

2. **Create skeleton summary file**
   - Read `templates/summary-skeleton.md`
   - Read `templates/skeleton-config.md` for article-type-specific placeholders
   - Fill in metadata: TITLE, TIMESTAMP, FILENAME, PAGE_COUNT, ARTICLE_TYPE
   - Apply article-type-specific customizations
   - Write skeleton to `{original_filename}_summary.md`

### Subagent Dispatch Phase

Process the paper section-by-section using specialized subagents:

| Order | Agent | PDF Sections | Output Sections |
|-------|-------|--------------|-----------------|
| 1 | `overview-agent` | Abstract, Introduction | Executive Summary, Background, Hypothesis |
| 2 | `methods-agent` | Methods section | Experimental Approach |
| 3 | `results-agent` | Results section | Key Findings, Figures, Statistics |
| 4 | `critique-agent` | Discussion | Strengths, Limitations, Red Flags |
| 5 | `context-agent` | Introduction + Discussion | Bigger Picture, Implications |
| 6 | `synthesis-agent` | **Summary file only** | Creative Insights, Synthesis |

**Subagent execution pattern:**

For each subagent:
1. Extract relevant PDF pages based on section mapping
2. Read the subagent prompt from `prompts/subagents/` or `subagents/`
3. Invoke Task tool with the subagent
4. Subagent reads PDF pages + current summary file
5. Subagent writes output to summary file, updating section markers from PENDING → COMPLETE
6. Verify section was written before proceeding to next agent

### Section Markers

The summary file uses HTML comment markers to track progress:

```html
<!-- SECTION: executive_summary PENDING -->   (not yet written)
<!-- SECTION: executive_summary COMPLETE -->  (successfully written)
<!-- SECTION: methods SKIPPED -->             (not applicable for article type)
```

### Article-Type-Specific Subagent Sequences

**General Research:**
```
overview → methods → results → critique → context → synthesis
```

**Review Article:**
```
overview → evidence → landscape → critique → references → synthesis
```

**Cell & Molecular Biology:**
```
overview → models → mechanisms → results → critique → translational → synthesis
```

**Computational Biology:**
```
overview → data → methods → validation → critique → reproducibility → synthesis
```

### Error Handling

| Scenario | Recovery |
|----------|----------|
| Subagent fails to write section | Retry once, then fall back to Standard Mode |
| Section detection fails | Use page-ratio fallback estimates |
| Context exceeded in subagent | Reduce page range, retry |

### Progress Reporting

Keep user informed during chunked processing:

```
Processing: smith2024_cancer_cells.pdf (24 pages)
Mode: Chunked subagent pipeline

[■■■□□□] Step 3/6: Analyzing results...
  ✓ Overview extracted
  ✓ Methods analyzed
  ◉ Results in progress
  ○ Critical analysis pending
  ○ Context pending
  ○ Synthesis pending
```

---

## Phase 2: Visual Output Format Selection

**After markdown summary is saved, ask the user which visual format(s) they want:**

| Format | Best For | Output |
|--------|----------|--------|
| **Interactive HTML** | Presentations, sharing, detailed review | Single-file HTML with navigation |
| **PDF Summary** | Archival, printing, formal documentation | Brand-compliant PDF via YAML intermediate |
| **SVG Infographic** | Visual presentations, posters, slides | Single-page graphical summary |

Ask: "Which visual output format(s) would you like?"
- Interactive HTML report (web presentation)
- PDF summary (printable document)
- SVG infographic (visual poster)
- Multiple formats
- None (markdown summary only)

## Core Principles (All Formats)

### Quantitative Preservation (Non-Negotiable)

NEVER summarize away numerical data. Always preserve:
- **P-values**: p < 0.05, p = 0.003 (with significance context)
- **Sample sizes**: n = 24/group, N = 156 total
- **Effect sizes**: Cohen's d = 0.8, HR = 0.45, OR = 2.3
- **Confidence intervals**: 95% CI: 1.2–3.4
- **Fold changes**: 2.3-fold increase
- **Concentrations/doses**: 1 mM DFMO, 5 mg/kg
- **Timepoints**: 24h, 4 weeks post-treatment

### Scientific Rigor

- Distinguish novel findings from confirmatory/background
- Preserve conditional language (authors' hedging)
- Include negative/null results
- Note controls and comparators
- Flag limitations (both author-stated and obvious)
- Use exact technical terminology (no paraphrasing jargon)

### Brand Visual Identity

All visual outputs follow Scientific Brand Visual Identity v4.0:

```
Primary Highlight:  #2DB2E8 (Brand Blue) - Use sparingly for PRIMARY finding only
Contrast/Alert:     #E8622D (Orange) - Opposing effects, warnings, limitations
Neutral Data:       #222222 (Dark Gray) - Controls, baseline, most data
Secondary:          #666666 (Medium Gray) - Secondary comparisons
Tertiary:           #999999 (Muted Gray) - Tertiary data
Annotations:        #BDBDBD (Light Gray) - NOT for data, only gridlines/annotations
Background:         #FFFFFF (White only)
Text:               #000000 (Black)
```

## Format-Specific Instructions

**IMPORTANT**: All visual formats use the saved markdown summary file as their source. Read the markdown summary first before generating any visual output.

### 1. Interactive HTML Report

Read: `prompts/html-report.md`

**Source**: The markdown summary file generated in Phase 1

Single-file HTML document with:
- Sticky navigation bar
- Collapsible sections
- Tabbed findings (In Vitro / In Vivo / Clinical)
- P-value badges with significance styling
- Mechanism diagrams (CSS/HTML)
- Verification checklist

### 2. PDF Summary (via YAML)

Read: `prompts/yaml-for-pdf.md`

**Source**: The markdown summary file generated in Phase 1

Two-step workflow:
1. Generate YAML summary with figure references (extracted from markdown summary)
2. User runs Python script to generate PDF with automatic figure extraction

Script location: `scripts/generate_summary_pdf.py`

### 3. SVG Infographic

Read: `prompts/svg-infographic.md`

**Source**: The markdown summary file generated in Phase 1

Single-page visual summary using component library:
- DonutChart, StatCallout, HorizontalBarChart
- TimelineNode, HierarchicalTree
- Brand-compliant colors only
- White background, flat design, no gradients

## Article Type Details

Each summarizer prompt is optimized for its article type:

| Article Type | Summarizer | Key Focus Areas |
|--------------|------------|-----------------|
| **General Research** | `general_researcher_summarizer.md` | Hypothesis, experimental approach, findings, critical analysis, actionable takeaways |
| **Review Article** | `review_article_summarizer.md` | Scope/framing, field landscape, evidence synthesis, knowledge gaps, reference mining |
| **Comp Bio / Bioinformatics** | `compbio_bioinformatic_summarizer.md` | Data foundation, computational approach, validation strategy, reproducibility assessment |
| **Cell & Molecular Biology** | `cellmolbio_summarizer.md` | Model systems, experimental design, mechanistic depth, cancer biology specifics, translational assessment |

### Special Handling Notes
| Condition | Handling |
|-----------|----------|
| **Preprint** | Add prominent "Not Peer-Reviewed" warning in all outputs |
| **Meta-analysis** | Use Review Article summarizer; include forest plot data, heterogeneity stats (I², Q) |
| **Clinical Trial** | Use General Research summarizer; emphasize CONSORT flow, primary/secondary endpoints, adverse events |
| **Methods Paper** | Use Comp Bio summarizer; prioritize protocol details, validation, benchmarks |

## Quick Reference

### Statistics Formatting
```
p = 0.003          (not "p<0.05" if exact value available)
95% CI: 1.2–3.4    (en-dash, not hyphen)
n = 24/group       (specify per-group)
HR = 0.45          (include measure type)
2.3-fold increase  (include direction)
```

### Color Assignment (Visual Outputs)
```
Treatment/Intervention → Brand Blue (#2DB2E8)
Control/Baseline → Dark Gray (#222222)
Adverse/Opposing effect → Orange (#E8622D)
Secondary comparison → Medium Gray (#666666)
```

---

## Troubleshooting & Error Recovery

### Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| Summary incomplete | Subagent failed mid-process | Retry the skill; check for PENDING markers in summary file |
| Context exceeded | Paper too long for chunked mode | Reduce page ranges; contact support for very long papers |
| Section detection failed | Non-standard paper format | Uses page-ratio fallback automatically |
| Missing statistics | Data in figures, not text | Check figure captions; may require manual extraction |
| Wrong article type | Misclassified paper | Re-run with correct article type selection |

### Processing Mode Issues

**Chunked mode not triggering for long papers:**
- Verify PDF is readable (not image-only scan)
- Check page count detection worked correctly
- Force chunked mode by confirming paper >12 pages

**Standard mode producing truncated output:**
- Paper may be near the 12-page threshold
- Dense text/tables increase token count
- Consider requesting chunked mode explicitly

### Subagent Pipeline Recovery

**If a subagent fails:**

1. Check the summary file for `<!-- SECTION: xxx PENDING -->` markers
2. Identify which section(s) failed
3. Claude will automatically:
   - Retry the failed agent once
   - Fall back to Standard Mode if retry fails
   - Report which sections couldn't be completed

**Partial completion recovery:**

```
# To resume from partial completion:
User: Continue summarizing this paper - some sections are still PENDING

Claude: [Reads existing summary file]
        Found 3 sections marked PENDING. Resuming from critique agent...
```

### PDF Reading Issues

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| Empty content extracted | Image-only PDF (scanned) | Use OCR tool first, then provide text |
| Garbled text | PDF encoding issues | Re-export from source or use different PDF |
| Missing sections | Multi-column layout confusion | Provide page hints for key sections |
| Figures not described | Text extraction excludes images | Refer to figure captions in PDF |

### Visual Output Issues

**HTML report not styling correctly:**
- Check browser console for errors
- Verify brand colors applied from `brand/BRAND_COLORS.md`
- Single-file HTML should work offline

**PDF generation fails:**
- Ensure weasyprint is installed: `pip install weasyprint`
- Check YAML syntax in intermediate file
- Run: `python scripts/generate_summary_pdf.py input.yaml`

**SVG infographic rendering issues:**
- Verify viewBox dimensions match content
- Check for unsupported CSS properties
- Test in multiple browsers

### When to Fall Back to Standard Mode

Standard mode (full-document) may work better when:
- Paper is well-structured with clear sections
- Paper is 10-15 pages (borderline length)
- You need faster processing (single pass)
- Chunked mode repeatedly fails

```
User: Use standard mode for this paper instead

Claude: Switching to Standard Mode (full-document analysis)...
```

### Getting Help

If issues persist:
1. Check that all required files exist (see Supporting Files section)
2. Verify the PDF is not corrupted
3. Try with a different paper to isolate the issue
4. Report issues at https://github.com/anthropics/claude-code/issues

---

## Supporting Files

### Article Type Summarizers (Standard Mode - Phase 1)
- `prompts/general_researcher_summarizer.md` - General primary research papers
- `prompts/review_article_summarizer.md` - Review articles and meta-analyses
- `prompts/compbio_bioinformatic_summarizer.md` - Computational biology and bioinformatics
- `prompts/cellmolbio_summarizer.md` - Cell and molecular biology papers

### Chunked Mode Infrastructure
- `templates/summary-skeleton.md` - Base skeleton template with section markers
- `templates/skeleton-config.md` - Article-type-specific placeholder values
- `prompts/subagents/_shared/section-detection.md` - PDF section boundary detection

### Chunked Mode Subagent Prompts (Phase 2 onwards)
Shared prompts (all article types):
- `prompts/subagents/_shared/overview.md` - Executive summary, background, hypothesis
- `prompts/subagents/_shared/critique.md` - Critical analysis, limitations
- `prompts/subagents/_shared/synthesis.md` - Final synthesis (reads summary file only)

Article-type-specific prompts in:
- `prompts/subagents/general/` - General research variants
- `prompts/subagents/review/` - Review article variants
- `prompts/subagents/cellmolbio/` - Cell & molecular biology variants
- `prompts/subagents/compbio/` - Computational biology variants

### Subagent Configuration Files
- `subagents/overview-agent.md` - Shared overview agent config
- `subagents/critique-agent.md` - Shared critique agent config
- `subagents/synthesis-agent.md` - Shared synthesis agent config
- `subagents/{type}/*.md` - Article-type-specific agent configs

### Visual Output Prompts (Phase 2)
- `prompts/html-report.md` - Interactive HTML prompt
- `prompts/yaml-for-pdf.md` - YAML generation for PDF
- `prompts/svg-infographic.md` - Infographic prompt with component library

### Brand Guidelines
- `brand/BRAND_COLORS.md` - Complete brand color palette and usage rules

### Scripts
- `scripts/generate_summary_pdf.py` - PDF generator script (requires weasyprint)

### Examples
- `examples/example_summary.yaml` - Example YAML output for PDF generation
- `examples/example_chunked_summary.md` - Complete markdown summary from chunked mode
- `examples/example_workflow_general.md` - Step-by-step general research workflow
- `examples/example_workflow_review.md` - Step-by-step review article workflow
