---
name: data-agent
description: Evaluates data sources, quality, and appropriateness for computational biology analyses
article_types: [compbio]
execution_order: 2
pdf_sections:
  - methods
  - materials_and_methods
  - data_availability
page_hints:
  start_percent: 10
  end_percent: 30
output_sections:
  - methods
  - methods_1
depends_on:
  - overview-agent
prompt_file: prompts/subagents/compbio/data.md
---

# Data Agent (Computational Biology)

## Purpose

The Data Agent critically evaluates the data foundation of computational biology papers. Data quality and appropriateness are fundamental to computational work - garbage in, garbage out.

## Input

- **PDF pages**: Methods/Data section, Data Availability statement
- **Article type**: compbio
- **Prior context**: Overview agent has established the research question

## Output

The agent writes to the summary file, filling in:

1. **Data Foundation** (main methods section)
   - Overview of all data sources
   - Data types and sizes

2. **Data Sources & Quality** (subsection)
   - Comprehensive data inventory table
   - Quality assessment per dataset
   - Batch effects and confounders

3. **Data Limitations** (subsection)
   - Acknowledged and unacknowledged limitations
   - Impact on conclusions

## Section Markers Updated

```
<!-- SECTION: methods PENDING --> → COMPLETE
<!-- SUBSECTION: methods_1 PENDING --> → COMPLETE
```

## Dependencies

Runs AFTER:
- overview-agent (order 1)

## Invocation Pattern

```
Task tool invocation:
- subagent_type: "general-purpose"
- prompt: Read prompts/subagents/compbio/data.md
- Input: Methods and Data sections from PDF
- Output: Updates summary file sections
```

## Quality Criteria

- ALL data sources documented with accession numbers
- Sample sizes critically assessed for intended analyses
- Batch effects and confounders identified
- Data accessibility assessed (can others get this data?)
- Missing data handling documented
- Limitations specific, not generic
