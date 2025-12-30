---
name: results-agent
description: Extracts key findings, figures, and statistics from Results section
article_types: [general]
execution_order: 3
pdf_sections:
  - results
page_hints:
  start_percent: 30
  end_percent: 65
output_sections:
  - results
  - key_findings
  - key_figures
  - statistics
depends_on:
  - overview-agent
  - methods-agent
prompt_file: prompts/subagents/general/results.md
---

# Results Agent (General Research)

## Purpose

The Results Agent extracts and organizes the key findings from general research papers. It ensures all quantitative data is preserved and identifies the most important figures and statistical information.

## Input

- **PDF pages**: Results section
- **Article type**: general
- **Prior context**: Overview established the hypothesis; Methods documented the approach

## Output

The agent writes to the summary file, filling in:

1. **Results** (main section)
   - Overview of findings organized by importance

2. **Key Findings** (subsection)
   - Each major finding with effect size, p-value, robustness assessment
   - Organized by importance to conclusions, not paper order

3. **Key Figures** (subsection)
   - 2-4 most important figures/tables
   - What each shows and key data points

4. **Statistical Summary** (subsection)
   - Sample size ranges
   - Effect size summary
   - P-value distribution
   - Tests used

## Section Markers Updated

```
<!-- SECTION: results PENDING --> → COMPLETE
<!-- SUBSECTION: key_findings PENDING --> → COMPLETE
<!-- SUBSECTION: key_figures PENDING --> → COMPLETE
<!-- SUBSECTION: statistics PENDING --> → COMPLETE
```

## Dependencies

Runs AFTER:
- overview-agent (order 1)
- methods-agent (order 2)

## Invocation Pattern

```
Task tool invocation:
- subagent_type: "general-purpose"
- prompt: Read prompts/subagents/general/results.md
- Input: Results section pages from PDF
- Output: Updates summary file sections
```

## Quality Criteria

- EVERY p-value, fold change, n, CI captured
- Effect directions specified (increase/decrease)
- Units included for all measurements
- Negative/null results explicitly documented
- Novel vs confirmatory findings distinguished
- Figure descriptions include specific data points, not just "shows X"
