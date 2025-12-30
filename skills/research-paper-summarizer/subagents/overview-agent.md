---
name: overview-agent
description: Extracts executive summary, background, and hypothesis from paper introduction
article_types: [general, review, cellmolbio, compbio]
execution_order: 1
pdf_sections:
  - abstract
  - introduction
page_hints:
  start: 1
  end: 3
  fallback_pages: 4
output_sections:
  - executive_summary
  - background
  - knowledge_gap
  - hypothesis
  - background_quote
depends_on: []
prompt_file: prompts/subagents/_shared/overview.md
---

# Overview Agent

## Purpose

The Overview Agent is the FIRST agent in the chunked processing pipeline. It reads the Abstract and Introduction sections to establish the foundational context for the paper summary.

## Input

- **PDF pages**: Abstract + Introduction (typically pages 1-3)
- **Article type**: Selected by user
- **Research focus**: User's research area (if provided)

## Output

The agent writes to the summary file, filling in:

1. **Executive Summary** - 2-3 sentence synopsis of the most important finding and why it matters
2. **Background & Motivation** - The knowledge gap being addressed
3. **Central Hypothesis** - Clear statement of the research question
4. **Key Quote** - Verbatim quote capturing authors' motivation

## Section Markers Updated

```
<!-- SECTION: executive_summary PENDING --> → COMPLETE
<!-- SECTION: background PENDING --> → COMPLETE
<!-- SUBSECTION: knowledge_gap PENDING --> → COMPLETE
<!-- SUBSECTION: hypothesis PENDING --> → COMPLETE
<!-- SUBSECTION: background_quote PENDING --> → COMPLETE
```

## Dependencies

None - this agent runs first.

## Invocation Pattern

```
Task tool invocation:
- subagent_type: "general-purpose"
- prompt: Read prompts/subagents/_shared/overview.md
- Input: PDF pages 1-3 (or detected abstract/intro pages)
- Output: Updates summary file sections
```

## Quality Criteria

- Executive summary captures THE key finding, not a list of findings
- Knowledge gap is specific, not generic ("cancer is bad")
- Hypothesis is clearly stated or explicitly noted as implicit
- Quote uses exact author words with quotation marks
- All numerical data from abstract is preserved
