# Phase 1 Foundation: Template Infrastructure

> **Status**: Approved | **Date**: 2025-12-28
> **Parent**: [INTEGRATION_ANALYSIS.md](../INTEGRATION_ANALYSIS.md)

**Goal**: Create template schemas, markdown templates, and TemplateParser for deterministic markdown-to-PDF conversion.

---

## Deliverables

| Item | Count | Location |
|------|-------|----------|
| YAML schemas | 12 | `src/oligon_reports/templates/schemas/*.yaml` |
| Markdown templates | 12 | `src/oligon_reports/templates/markdown/*.md` |
| TemplateParser class | 1 | `src/oligon_reports/templates/parser.py` |
| Module init | 2 | `templates/__init__.py`, update main `__init__.py` |

---

## Files to Create/Modify

### New Files

```
src/oligon_reports/templates/
├── __init__.py              # Export TemplateParser, DocumentTree, etc.
├── parser.py                # TemplateParser class (~300 lines)
├── schemas/
│   ├── analysis-report.yaml
│   ├── phase-plan.yaml
│   ├── data-report.yaml
│   ├── literature-review.yaml
│   ├── meeting-notes.yaml
│   ├── project-status.yaml
│   ├── technical-spec.yaml
│   ├── task-list.yaml
│   ├── standards-guide.yaml
│   ├── agent-definition.yaml
│   ├── readme.yaml
│   └── method-guide.yaml
└── markdown/
    └── [same 12 files as .md]
```

### Modify

- `pyproject.toml` - Add `python-frontmatter>=1.1.0`
- `src/oligon_reports/__init__.py` - Export template classes

---

## Schema Format

Each YAML schema defines:

```yaml
schema_version: "1.0"
document_type: "analysis-report"
category: "scientific"
description: "..."

frontmatter:
  required: [type, title, date]
  optional: [project, phase, author, version, assessment]

sections:
  - id: objective
    heading: "Objective"
    level: 2
    required: true
    subsections: [...]
    expected_elements: [tables, code_blocks]

component_mappings:
  "## N. Title": "SectionDivider"
  "> **Type:**": "CalloutBox"
  "#### Finding N:": "FindingCard"

validation:
  require_frontmatter: true
  min_sections: 3
```

---

## TemplateParser API

```python
@dataclass
class DocumentTree:
    frontmatter: FrontmatterData
    sections: list[Section]
    raw_content: str

class TemplateParser:
    def __init__(self, template_type: str | None = None): ...
    def parse(self, content: str) -> DocumentTree: ...
    def validate(self, tree: DocumentTree) -> list[ValidationError]: ...
    def detect_type(self, content: str) -> str | None: ...

    @staticmethod
    def list_templates() -> list[dict[str, str]]: ...

    @staticmethod
    def get_template(template_type: str) -> str: ...
```

---

## Implementation Order

### Step 1: Setup
- [ ] Add `python-frontmatter>=1.1.0` to pyproject.toml
- [ ] Run `uv sync`
- [ ] Create directory structure

### Step 2: Parser Foundation
- [ ] Create `templates/__init__.py` with dataclass definitions
- [ ] Create `parser.py` with:
  - `_parse_frontmatter()` using python-frontmatter
  - `detect_type()` method
  - `_load_schema()` using importlib.resources
  - Skeleton `parse()` and `validate()` methods

### Step 3: Schemas - Batch 1
- [ ] `analysis-report.yaml` (complex, establishes patterns)
- [ ] `meeting-notes.yaml` (simple, quick validation)
- [ ] `readme.yaml` (common reference type)
- [ ] `data-report.yaml` (data-driven variant)

### Step 4: Schemas - Batch 2
- [ ] `phase-plan.yaml`
- [ ] `literature-review.yaml`
- [ ] `project-status.yaml`
- [ ] `technical-spec.yaml`

### Step 5: Schemas - Batch 3
- [ ] `task-list.yaml`
- [ ] `standards-guide.yaml`
- [ ] `agent-definition.yaml`
- [ ] `method-guide.yaml`

### Step 6: Markdown Templates
- [ ] Create full-example templates for all 12 types
- [ ] Based on DOCUMENT_TEMPLATING_SYSTEM.md Section 4

### Step 7: Complete Parser
- [ ] Implement `_parse_sections()` with heading regex
- [ ] Implement `_detect_elements()` for tables, callouts, code blocks
- [ ] Implement `validate()` against schema rules
- [ ] Implement `list_templates()` and `get_template()`

### Step 8: Integration
- [ ] Update `src/oligon_reports/__init__.py` with exports
- [ ] Verify `uv run python -c "from oligon_reports import TemplateParser"`

---

## Document Types (12 total)

| Category | Types |
|----------|-------|
| Scientific | analysis-report, phase-plan, data-report |
| Project Mgmt | literature-review, meeting-notes, project-status |
| Development | technical-spec, task-list, standards-guide, agent-definition |
| Reference | readme, method-guide |

---

## Source of Truth

- Template structures: `docs/template-project/brand/DOCUMENT_TEMPLATING_SYSTEM.md`
- Component mappings: Same document, Section 3
- Existing components: `src/oligon_reports/components.py`

---

## Success Criteria

1. All 12 schemas loadable via `TemplateParser("type-name")`
2. All 12 templates retrievable via `TemplateParser.get_template("type-name")`
3. `detect_type()` correctly identifies type from any template
4. `parse()` returns valid DocumentTree for sample documents
5. `validate()` catches missing required sections/frontmatter
6. Clean import: `from oligon_reports import TemplateParser`

---

## Detailed Schema Specification

### Example: analysis-report.yaml

```yaml
schema_version: "1.0"
document_type: "analysis-report"
category: "scientific"
description: |
  Phase-style analysis documents with objective, methods, results, discussion.
  Supports optional PASS/FAIL assessment mode.

frontmatter:
  required:
    - type: { expected: "analysis-report" }
    - title: { type: string }
    - date: { type: string, format: "YYYY-MM-DD" }
  optional:
    - project: { type: string }
    - phase: { type: string }
    - focus: { type: string }
    - author: { type: string }
    - version: { type: string, default: "1.0" }
    - assessment: { type: enum, values: ["none", "pass-fail"], default: "none" }

sections:
  - id: objective
    heading: "Objective"
    level: 2
    required: true
    subsections:
      - id: key_questions
        heading: "Key Questions"
        required: false
        expected_content: { type: "numbered_list" }
      - id: critical_context
        heading: "Critical Context"
        required: false
        expected_content: { type: "callout" }

  - id: methods
    heading: "Methods"
    level: 2
    required: true
    allows_dynamic_subsections: true
    expected_elements: [tables, code_blocks, method_blocks]

  - id: results
    heading: "Results"
    level: 2
    required: true
    allows_multiple: true

  - id: summary
    heading: "Summary"
    level: 2
    required: true
    subsections:
      - id: key_findings
        heading: "Key Findings"
        expected_content: { type: "finding_cards" }

  - id: discussion
    heading: "Discussion"
    level: 2
    required: false

  - id: recommendations
    heading: "Recommendations"
    level: 2
    required: false

component_mappings:
  frontmatter: "CoverPage"
  "## N. Title": "SectionDivider"
  "> **Type:**": "CalloutBox"
  "#### Finding N:": "FindingCard"
  "table_with_status": "StatusTable"
  "- [ ]": "ChecklistItem"
  "code_block": "CodeBlock"

validation:
  require_frontmatter: true
  min_sections: 3
  allow_extra_sections: true
```

---

## Parser Dataclass Definitions

```python
from dataclasses import dataclass, field
from typing import Any

@dataclass
class FrontmatterData:
    type: str
    title: str
    date: str
    raw: dict[str, Any] = field(default_factory=dict)

@dataclass
class Element:
    type: str  # "table", "callout", "code_block", etc.
    content: str
    attributes: dict[str, Any] = field(default_factory=dict)

@dataclass
class Section:
    id: str
    heading: str
    level: int
    content: str
    subsections: list["Section"] = field(default_factory=list)
    elements: list[Element] = field(default_factory=list)

@dataclass
class DocumentTree:
    frontmatter: FrontmatterData
    sections: list[Section]
    raw_content: str

@dataclass
class ValidationError:
    level: str  # "error" | "warning"
    message: str
    location: str | None = None
```

---

## Dependencies

Add to pyproject.toml:
```toml
dependencies = [
    "python-frontmatter>=1.1.0",
]
```
