---
name: critique-agent
description: Performs critical analysis identifying strengths, limitations, and red flags
article_types: [general, review, cellmolbio, compbio]
execution_order: 4
pdf_sections:
  - discussion
  - conclusions
page_hints:
  start_percent: 65
  end_percent: 85
output_sections:
  - critical_analysis
  - strengths
  - limitations
  - red_flags
depends_on:
  - overview-agent
  - methods-agent
  - results-agent
prompt_file: prompts/subagents/_shared/critique.md
---

# Critique Agent

## Purpose

The Critique Agent performs rigorous critical analysis of the paper, identifying both strengths and weaknesses. It reads the Discussion section where authors interpret their findings and acknowledge limitations.

## Input

- **PDF pages**: Discussion and Conclusions sections
- **Article type**: Selected by user
- **Prior context**: Can reference what overview/methods/results agents found

## Output

The agent writes to the summary file, filling in:

1. **Strengths** - What the paper does well (methodological rigor, novel contributions, quality of evidence)
2. **Limitations & Gaps** - Missing experiments, control issues, confounding variables, generalizability concerns
3. **Red Flags** - Data presentation concerns, statistical issues, logical gaps, overclaiming

## Section Markers Updated

```
<!-- SECTION: critical_analysis PENDING --> → COMPLETE
<!-- SUBSECTION: strengths PENDING --> → COMPLETE
<!-- SUBSECTION: limitations PENDING --> → COMPLETE
<!-- SUBSECTION: red_flags PENDING --> → COMPLETE
```

## Dependencies

Runs AFTER:
- overview-agent (order 1)
- methods-agent (order 2)
- results-agent (order 3)

This allows the critique to be informed by what was found in methods and results.

## Invocation Pattern

```
Task tool invocation:
- subagent_type: "general-purpose"
- prompt: Read prompts/subagents/_shared/critique.md
- Input: Discussion/Conclusions pages from PDF
- Output: Updates summary file sections
```

## Quality Criteria

- Strengths are specific, not generic praise
- Limitations explain WHY they matter and their impact on conclusions
- Red flags are evidence-based, not speculative
- Distinguishes major vs minor concerns
- Includes author-acknowledged limitations with quotes
- "No red flags" stated explicitly if none found (rather than inventing concerns)
