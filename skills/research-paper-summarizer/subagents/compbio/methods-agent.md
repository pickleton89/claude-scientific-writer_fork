---
name: methods-agent
description: Documents computational pipeline, algorithms, and statistical framework
article_types: [compbio]
execution_order: 3
pdf_sections:
  - methods
  - supplementary_methods
page_hints:
  start_percent: 15
  end_percent: 40
output_sections:
  - methods_2
depends_on:
  - overview-agent
  - data-agent
prompt_file: prompts/subagents/compbio/methods.md
---

# Methods Agent (Computational Biology)

## Purpose

The Methods Agent documents and critically evaluates all computational methods, algorithms, and analytical approaches. It assesses whether the methods are appropriate for the question and identifies potential sources of bias.

## Input

- **PDF pages**: Methods section, Supplementary Methods
- **Article type**: compbio
- **Prior context**: Overview established question; Data documented inputs

## Output

The agent writes to the summary file, filling in:

1. **Computational Pipeline & Algorithms** (methods_2 subsection)
   - Pipeline overview with flowchart
   - Methods inventory table (tools, versions, parameters)
   - Algorithm assessments
   - Statistical framework (tests, ML approaches)
   - Parameter choices and potential biases

## Section Markers Updated

```
<!-- SUBSECTION: methods_2 PENDING --> → COMPLETE
```

## Dependencies

Runs AFTER:
- overview-agent (order 1)
- data-agent (order 2)

## Invocation Pattern

```
Task tool invocation:
- subagent_type: "general-purpose"
- prompt: Read prompts/subagents/compbio/methods.md
- Input: Methods section pages from PDF
- Output: Updates summary file sections
```

## Quality Criteria

- EVERY tool documented with version (or "not stated")
- Parameters assessed as default/optimized/arbitrary
- ML workflows checked for data leakage
- Multiple testing correction verified
- Potential biases explicitly identified
- Alternative methods noted where relevant
