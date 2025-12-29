# Phase 3: Skill Creation Implementation Plan

**Date:** December 28, 2025
**Goal:** Create unified `markdown-to-pdf` skill with three commands

---

## Overview

Create a skill that leverages Phase 1's `TemplateParser` and Phase 2's components to provide a complete markdown-to-PDF workflow.

---

## Architecture

```
skills/markdown-to-pdf/
├── SKILL.md              # Main skill definition
└── references/
    └── component_map.md  # Template element → PDF component mapping
```

**Dependencies:**
- `TemplateParser` from `src/oligon_reports/templates/parser.py`
- PDF components from `src/oligon_reports/components.py`
- `ReportGenerator` from `src/oligon_reports/report_generator.py`

---

## Commands (Inline in SKILL.md)

### 1. `/list-templates`

**Purpose:** Show available document templates

**Workflow:**
1. Call `TemplateParser.list_templates()`
2. Display formatted table with type, category, description

**Output:**
```
| Type              | Category        | Description                           |
|-------------------|-----------------|---------------------------------------|
| analysis-report   | scientific      | Structured analysis with findings     |
| literature-review | scientific      | Systematic literature review          |
| meeting-notes     | project-mgmt    | Meeting documentation                 |
...
```

### 2. `/new-doc <type>`

**Purpose:** Create new document from template

**Workflow:**
1. Validate `<type>` exists in template list
2. Call `TemplateParser.get_template(type)`
3. Create new file: `<type>_YYYYMMDD.md` in current directory
4. Open for editing

**Example:**
```
> /new-doc analysis-report

Created: analysis-report_20251228.md
Template type: analysis-report (scientific)
```

### 3. `/doc-to-pdf <file>`

**Purpose:** Convert markdown to branded PDF with confirmation

**Workflow:**

**Step 1: Parse & Detect**
```python
parser = TemplateParser()
tree = parser.parse(content)
```

**Step 2: Show Structure for Confirmation**
```
📄 Document: analysis-report_20251228.md
📋 Type: analysis-report (detected from frontmatter)

Structure:
├── 1. Objective
├── 2. Key Questions
├── 3. Methodology
│   └── 3.1 Data Sources
├── 4. Findings
│   ├── Finding 1: [title]
│   ├── Finding 2: [title]
│   └── Finding 3: [title]
└── 5. Recommendations

Elements Detected:
- 3 tables (→ StatusTable)
- 3 finding cards (→ FindingCard)
- 2 code blocks (→ MethodBlock)

⚠️ Validation:
- ✅ All required sections present
- ⚠️ Missing optional: executive_summary

Proceed with PDF generation? [Y/n/edit]
```

**Step 3: Generate PDF**
- Map elements to components (see component_map.md)
- Generate using `ReportGenerator`
- Output to `output/<filename>.pdf`
- Create `output/` if needed, add to `.gitignore`

**Step 4: Optional Figure Enhancement**
If user requests, invoke `scientific-schematics` skill:
```
> /doc-to-pdf analysis-report.md --figures

Generating figures for document...
[Invokes scientific-schematics for each figure placeholder]
```

---

## Component Mapping

| Markdown Element | Detection Pattern | PDF Component |
|------------------|-------------------|---------------|
| Tables with ✓/✗ | Cell contains ✓, ✗, ✔, ✘ | `StatusTable` |
| Tables with grades | Cell contains A-F or tier values | `GradedTable` |
| Regular tables | Standard markdown table | `BaseTable` |
| `#### Finding N:` | Heading pattern | `FindingCard` |
| Code blocks | Triple backticks | `MethodBlock` or code |
| Callouts `> **Type:**` | Blockquote with bold type | `CalloutBox` |
| Checklists `- [ ]` | Checkbox pattern | Checklist component |
| Bullet/numbered lists | Standard list markers | List rendering |

---

## Success Criteria

1. `/list-templates` shows all 12 template types
2. `/new-doc` creates valid templated files
3. `/doc-to-pdf` shows structure confirmation before generation
4. PDF output matches brand styling from `oligon_reports`
5. `--figures` flag triggers `scientific-schematics` integration

---

## Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Command location | Inline in SKILL.md | Single source of truth, matches document-skills pattern |
| Confirmation format | ASCII tree | Shows hierarchy, readable in terminal |
| Output location | `output/` folder | Clean separation, gitignored |
| Figure generation | Opt-in via `--figures` | User control, not automatic |

---

## File Changes

| File | Action |
|------|--------|
| `skills/markdown-to-pdf/SKILL.md` | CREATE |
| `skills/markdown-to-pdf/references/component_map.md` | CREATE |
| `.gitignore` | EDIT (add output/) |
| `INTEGRATION_ANALYSIS.md` | EDIT (mark Phase 3 complete) |
| `CHANGELOG.md` | EDIT (document Phase 3) |

---

## Next Steps

1. Review this plan
2. Create SKILL.md with decision tree and command definitions
3. Create component_map.md reference
4. Test each command manually
5. Update documentation

---

*Phase 3 of INTEGRATION_ANALYSIS.md implementation*
