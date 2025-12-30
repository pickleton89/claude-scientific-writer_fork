---
name: synthesis-agent
description: Creates final synthesis by reading completed summary file (does NOT read PDF)
article_types: [general, review, cellmolbio, compbio]
execution_order: 6
pdf_sections: []
reads_summary_only: true
output_sections:
  - insights
  - takeaways
  - synthesis
depends_on:
  - overview-agent
  - methods-agent
  - results-agent
  - critique-agent
  - context-agent
prompt_file: prompts/subagents/_shared/synthesis.md
---

# Synthesis Agent

## Purpose

The Synthesis Agent is the FINAL agent in the pipeline. Unlike other agents, it does NOT read the PDF. Instead, it reads the completed summary file and generates integrative insights that tie everything together.

## Input

- **Summary file**: The accumulated output from all prior agents
- **Article type**: Selected by user
- **Research focus**: User's research area for tailored takeaways

**IMPORTANT**: This agent receives NO PDF content.

## Output

The agent writes to the summary file, completing:

1. **Creative Insights & Connections** - Cross-field applications, unexpected connections, speculation about follow-up findings
2. **Actionable Takeaways** - Techniques to adapt, targets to consider, pitfalls to avoid (tailored to research focus)
3. **One-Paragraph Synthesis** - 150-200 word "coffee chat" summary integrating key finding, approach, significance, and critical perspective

## Section Markers Updated

```
<!-- SECTION: insights PENDING --> → COMPLETE
<!-- SECTION: takeaways PENDING --> → COMPLETE
<!-- SECTION: synthesis PENDING --> → COMPLETE
```

## Dependencies

Runs LAST, after ALL other agents complete:
- overview-agent (order 1)
- methods-agent (order 2)
- results-agent (order 3)
- critique-agent (order 4)
- context-agent (order 5)

## Invocation Pattern

```
Task tool invocation:
- subagent_type: "general-purpose"
- prompt: Read prompts/subagents/_shared/synthesis.md
- Input: ONLY the summary file path (no PDF content)
- Output: Completes final sections of summary file
```

## Quality Criteria

- Insights are scientifically grounded, not wild speculation
- Cross-field connections are specific and plausible
- Takeaways are tailored to user's research focus
- One-paragraph synthesis reads naturally (not a bullet list in disguise)
- Synthesis integrates information from all prior sections
- No repetition of content already written by other agents
