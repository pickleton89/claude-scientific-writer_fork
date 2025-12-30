---
name: context-agent
description: Places findings in broader scientific context and assesses translational implications
article_types: [general]
execution_order: 5
pdf_sections:
  - introduction
  - discussion
page_hints:
  introduction_pages: [1, 2, 3]
  discussion_start_percent: 65
  discussion_end_percent: 85
output_sections:
  - context
  - field_impact
  - future_directions
  - article_specific
  - clinical_relevance
  - followup_studies
depends_on:
  - overview-agent
  - methods-agent
  - results-agent
  - critique-agent
prompt_file: prompts/subagents/general/context.md
---

# Context Agent (General Research)

## Purpose

The Context Agent places research findings in their broader scientific context. It analyzes how the work fits into the field, identifies competing approaches, and assesses translational potential.

## Input

- **PDF pages**: Introduction (for field background) + Discussion (for interpretation)
- **Article type**: general
- **Prior context**: Full picture from overview, methods, results, critique agents

## Output

The agent writes to the summary file, filling in:

1. **Bigger Picture** (main section)
   - How this advances the field
   - Paradigm confirmation/extension/challenge

2. **Field Impact** (subsection)
   - Competing approaches and how this compares
   - Controversies addressed

3. **Future Directions** (subsection)
   - Immediate next steps
   - Key questions opened
   - Long-term implications

4. **Implications & Applications** (article-specific section)
   - Clinical/practical relevance
   - Suggested follow-up studies

## Section Markers Updated

```
<!-- SECTION: context PENDING --> → COMPLETE
<!-- SUBSECTION: field_impact PENDING --> → COMPLETE
<!-- SUBSECTION: future_directions PENDING --> → COMPLETE
<!-- SECTION: article_specific PENDING --> → COMPLETE
<!-- SUBSECTION: clinical_relevance PENDING --> → COMPLETE
<!-- SUBSECTION: followup_studies PENDING --> → COMPLETE
```

## Dependencies

Runs AFTER:
- overview-agent (order 1)
- methods-agent (order 2)
- results-agent (order 3)
- critique-agent (order 4)

This allows context to integrate the full picture.

## Invocation Pattern

```
Task tool invocation:
- subagent_type: "general-purpose"
- prompt: Read prompts/subagents/general/context.md
- Input: Introduction + Discussion pages from PDF
- Output: Updates summary file sections
```

## Quality Criteria

- Competing work cited with references mentioned in paper
- Assessment of contribution magnitude is honest (incremental vs transformative)
- Future directions are specific experiments, not vague "more work"
- Translational path includes concrete barriers and realistic timeline
- Distinguishes author claims from analyst assessment
