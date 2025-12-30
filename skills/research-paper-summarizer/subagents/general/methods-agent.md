---
name: methods-agent
description: Extracts experimental approach, model systems, and controls from Methods section
article_types: [general]
execution_order: 2
pdf_sections:
  - methods
  - materials_and_methods
page_hints:
  start_percent: 15
  end_percent: 30
output_sections:
  - methods
  - methods_1
  - methods_2
depends_on:
  - overview-agent
prompt_file: prompts/subagents/general/methods.md
---

# Methods Agent (General Research)

## Purpose

The Methods Agent extracts detailed methodological information from general research papers. It documents study design, model systems, key techniques, and critically evaluates the control strategy.

## Input

- **PDF pages**: Methods/Materials and Methods section
- **Article type**: general
- **Prior context**: Overview agent has established the research question

## Output

The agent writes to the summary file, filling in:

1. **Experimental Approach** (main section)
   - Study type and design structure
   - Key comparisons being made

2. **Study Design** (subsection)
   - Model systems table (cell lines, animal models, patient samples)
   - Sources and catalog numbers where available

3. **Key Methods & Controls** (subsection)
   - Major techniques with purpose, parameters, validation
   - Control assessment (adequate vs missing)
   - Sample sizes and statistical approaches

## Section Markers Updated

```
<!-- SECTION: methods PENDING --> → COMPLETE
<!-- SUBSECTION: methods_1 PENDING --> → COMPLETE
<!-- SUBSECTION: methods_2 PENDING --> → COMPLETE
```

## Dependencies

Runs AFTER:
- overview-agent (order 1)

## Invocation Pattern

```
Task tool invocation:
- subagent_type: "general-purpose"
- prompt: Read prompts/subagents/general/methods.md
- Input: Methods section pages from PDF
- Output: Updates summary file sections
```

## Quality Criteria

- ALL numerical details preserved (concentrations, timepoints, n values)
- Model systems include sources when provided
- Techniques linked to the questions they answer (WHY chosen)
- Controls critically assessed, not just listed
- Missing controls explicitly identified
- Sample sizes clearly stated per group
