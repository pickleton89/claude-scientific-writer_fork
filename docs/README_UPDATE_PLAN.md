# README.md Update Plan

> Generated: 2025-12-29
> Status: Ready for implementation

## Executive Summary

The README.md requires significant updates to accurately reflect the current state of the project. Key changes include:

1. **Correct skill count** from 19 to 22 top-level skills
2. **Add missing skills** to the skills table (4 skills missing)
3. **Add new sections** for testing, recent changes, and architecture
4. **Update project structure** to reflect current directory layout
5. **Expand components list** from 5 to full set
6. **Add quality badges** for professional appearance

**Impact**: These changes will make the README accurate, comprehensive, and valuable for users.

---

## Analysis Summary

### Current State Issues

| Issue | Severity | Current | Actual |
|-------|----------|---------|--------|
| Skill count | Critical | 19 top-level | 22 top-level |
| Skills in table | Critical | 17 skills | 22 skills |
| Components listed | Medium | 5 | 11 |
| Testing section | Medium | Missing | 24 tests exist |
| Recent changes | Medium | Not highlighted | Major work done |

### Missing Skills in README

These skills exist but aren't in the skills table:
- `code-documentation` (bioinformatics)
- `reproducible-research` (bioinformatics)
- `statistical-analysis` (bioinformatics)
- `markdown-to-pdf` (mentioned later but not in table)

### Major Work Not Documented

1. Visual Design Skill Inheritance (Option D)
2. Skills Determinism Audit (7.5/10 → 8.5/10)
3. Bioinformatics skills implementation
4. SKILL_ROUTER.md for multi-skill routing
5. 24-test integration suite

---

## Section 1: Structural Changes

### 1.1 Add Badges Section
**Priority**: Medium
**Location**: After title, before first paragraph

```markdown
![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)
![Tests](https://img.shields.io/badge/tests-24%20passing-brightgreen.svg)
![Skills](https://img.shields.io/badge/skills-22%20active-purple.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
```

### 1.2 Reorganize Sections
**Priority**: High

**Proposed section order:**
1. Title and badges
2. Introduction (existing)
3. What This Is (existing, minor updates)
4. **NEW: What's New** - Recent improvements summary
5. Key Customizations (existing, update counts)
6. Skills (existing, major updates)
7. **NEW: Architecture** - Skill inheritance, router
8. Document Template System (existing, minor updates)
9. Components (existing, expand list)
10. Project Structure (existing, update tree)
11. Getting Started (existing)
12. **NEW: Testing** - How to run tests
13. Development (existing)
14. Documentation (existing, add new docs)
15. Upstream Reference (existing)
16. License (existing)

---

## Section 2: Content Updates

### 2.1 Update Skills Section
**Priority**: Critical

#### 2.1.1 Fix Skill Count in Header
**Location**: Line 27
**Current**: `## Skills (19 top-level + 4 document sub-skills)`
**Change to**: `## Skills (22 top-level + 4 document sub-skills)`

#### 2.1.2 Update Skills Table
**Location**: Lines 29-37

**Replace with:**
```markdown
| Category | Skills |
|----------|--------|
| **Writing** | `scientific-writing`, `literature-review`, `hypothesis-generation` |
| **Presentations** | `scientific-slides`, `latex-posters`, `pptx-posters` |
| **Research** | `research-lookup`, `citation-management`, `peer-review`, `scholar-evaluation` |
| **Visuals** | `scientific-schematics`, `generate-image`, `plotting-libraries`, `visual-design` |
| **Documents** | `markdown-to-pdf`, `markitdown`, `venue-templates`, `document-skills/*` (docx, pdf, pptx, xlsx) |
| **Analysis** | `statistical-analysis`, `reproducible-research`, `code-documentation` |
| **Conversion** | `paper-2-web` |
```

**Changes made:**
- Added `markdown-to-pdf` to Documents
- New "Analysis" category with 3 bioinformatics skills
- Removed "scientific-critical-thinking" (merged into peer-review)

### 2.2 Add "What's New" Section
**Priority**: High
**Location**: After "Key Customizations" section

```markdown
## What's New

### Visual Design Skill Inheritance
All 6 visual output skills now inherit from `visual-design` as a parent skill, providing consistent design philosophy across:
- `scientific-schematics`, `scientific-slides`, `latex-posters`, `pptx-posters`, `generate-image`, `plotting-libraries`

### Bioinformatics Skills (New)
Three new skills for computational biology workflows:
- **`statistical-analysis`** - Test selection, multiple testing correction, normalization methods
- **`reproducible-research`** - FAIR data principles, containers, environment management
- **`code-documentation`** - Docstrings, READMEs, API documentation

### Skills Determinism Audit
Library-wide quality audit improving determinism score from 7.5/10 to **8.5/10**:
- Standardized all skills with decision matrices and exit criteria
- Added `SKILL_ROUTER.md` for multi-skill task routing
- Created shared `QUANTIFICATION_THRESHOLDS.md` for measurable standards

### Integration Testing
24 integration tests covering the full template-to-PDF workflow.
```

### 2.3 Add Architecture Section
**Priority**: Medium
**Location**: After Skills section

```markdown
## Architecture

### Skill Inheritance
Visual output skills follow a parent-child architecture:

```
visual-design (parent - design philosophy)
├── scientific-schematics
├── scientific-slides
├── latex-posters
├── pptx-posters
├── generate-image
└── plotting-libraries
```

Each child skill includes `extends: visual-design` in its frontmatter and references the parent for design decisions.

### Skill Router
`skills/SKILL_ROUTER.md` provides decision trees for routing tasks to the correct skill. Covers:
- Figure & Visual Creation
- Document Creation
- Research Workflow
- Analysis & Reproducibility
- Evaluation/Review
- Conversion

### Template System
The `oligon_reports` package provides a full template pipeline:

```
Markdown → TemplateParser → DocumentTree → ReportGenerator → Branded PDF
```

- 12 YAML schemas define document structure
- 12 markdown templates provide starting points
- 11 visual components render to PDF
```

### 2.4 Update Key Customizations Table
**Priority**: High
**Location**: Lines 19-25

**Change row:**
```markdown
| **Skills Streamlined** | Removed 5 clinical/business skills; added 4 new skills (22 total) |
```

### 2.5 Update Components Section
**Priority**: Medium
**Location**: Lines 63-72

**Replace with:**
```markdown
### Components

The `oligon_reports` package provides 11 branded PDF components:

| Component | Purpose |
|-----------|---------|
| **FindingCard** | Numbered finding boxes with badge and description |
| **StatusTable** | Tables with ✓/✗ color-coded cells |
| **GradedTable** | Tables with tier-based color bands (gold/silver/bronze) |
| **MethodBlock** | What/Why/How structured blocks |
| **CalloutBox** | Styled callouts (Note, Warning, Example) |
| **MetadataHeader** | Document header with type badge and metadata grid |
| **MetricCard** | Key metric display cards |
| **Timeline** | Visual timeline components |
| And more... | See `src/oligon_reports/components.py` |
```

### 2.6 Update Project Structure
**Priority**: High
**Location**: Lines 73-95

**Replace with:**
```markdown
## Project Structure

```
claude-scientific-writer_fork/
├── skills/                  # Canonical skill definitions (22 + 4 sub-skills)
│   ├── SKILL_ROUTER.md      # Multi-skill task routing
│   ├── SKILL_TEMPLATE.md    # Standard skill template
│   ├── SKILL_TESTS.md       # Validation test scenarios
│   ├── QUANTIFICATION_THRESHOLDS.md  # Shared numeric thresholds
│   ├── markdown-to-pdf/     # Template-based PDF generation
│   ├── statistical-analysis/  # Statistical methods
│   ├── code-documentation/  # Code documentation
│   ├── reproducible-research/  # FAIR data practices
│   └── [18 more skills...]
├── .claude/
│   └── WRITER.md            # Agent system instructions
├── src/oligon_reports/      # Branded PDF generation (ReportLab)
│   ├── components.py        # 11 visual components
│   ├── report_generator.py  # PDF orchestrator
│   └── templates/           # Document templating system
│       ├── schemas/         # YAML validation schemas (12 types)
│       ├── markdown/        # Markdown templates (12 types)
│       └── parser.py        # TemplateParser class
├── tests/                   # Integration test suite (24 tests)
│   ├── conftest.py          # Pytest fixtures
│   ├── test_integration.py  # Full workflow tests
│   └── fixtures/            # Sample documents
├── docs/
│   ├── original/            # Archived upstream documentation
│   ├── template-project/    # Oligon brand standards
│   └── archived/            # Completed planning documents
├── prompts/                 # Meta-prompts for skill generation
├── scientific_writer/       # Python package (CLI/API)
├── examples/                # Working examples
├── CLAUDE.md                # Development documentation
└── CHANGELOG.md             # Fork development history
```
```

### 2.7 Add Testing Section
**Priority**: Medium
**Location**: Before or after Development section

```markdown
## Testing

The project includes a comprehensive integration test suite covering the template-to-PDF workflow.

### Running Tests

```bash
# Run all tests
PYTHONPATH=src uv run pytest tests/ -v

# Run specific test categories
PYTHONPATH=src uv run pytest tests/test_integration.py::TestTemplateDiscovery -v
PYTHONPATH=src uv run pytest tests/test_integration.py::TestPDFGeneration -v
```

### Test Coverage

| Category | Tests | Coverage |
|----------|-------|----------|
| Template Discovery | 4 | Template listing and retrieval |
| Document Parsing | 3 | Frontmatter, sections, type detection |
| Element Detection | 5 | Tables, code blocks, findings, checklists |
| Validation | 2 | Schema validation |
| Component Mapping | 5 | Element → Component mapping |
| PDF Generation | 2 | End-to-end PDF workflow |
| Cross-Skill Integration | 3 | Package exports and API |
```

### 2.8 Update Documentation Table
**Priority**: Medium
**Location**: Lines 155-161

**Replace with:**
```markdown
| Document | Description |
|----------|-------------|
| `CLAUDE.md` | Development documentation and workflow |
| `CHANGELOG.md` | Fork customization history |
| `.claude/WRITER.md` | Scientific writing agent instructions |
| `skills/SKILL_ROUTER.md` | Multi-skill task routing decision trees |
| `skills/SKILL_TEMPLATE.md` | Standard template for new skills |
| `skills/QUANTIFICATION_THRESHOLDS.md` | Shared numeric thresholds |
| `docs/archived/INTEGRATION_ANALYSIS.md` | Completed template-project merge |
```

---

## Section 3: Implementation Steps

### Step-by-Step Checklist

| Step | Priority | Action | Dependencies |
|------|----------|--------|--------------|
| 1 | Medium | Add badges after title | None |
| 2 | Critical | Fix skill count (19→22) | None |
| 3 | Critical | Update skills table | Step 2 |
| 4 | High | Update Key Customizations count | None |
| 5 | High | Add "What's New" section | None |
| 6 | Medium | Add Architecture section | Step 3 |
| 7 | Medium | Update Components section | None |
| 8 | High | Update Project Structure | None |
| 9 | Medium | Add Testing section | None |
| 10 | Medium | Update Documentation table | None |

### Recommended Order
1. Steps 2, 3 (Critical - fix skill count/table)
2. Steps 4, 5, 8 (High - counts, What's New, structure)
3. Steps 1, 6, 7, 9, 10 (Medium - polish)

---

## Validation Checklist

After implementing changes, verify:
- [ ] All 22 top-level skills are listed or mentioned
- [ ] Skill count matches actual SKILL.md files (26 total)
- [ ] Project structure tree matches actual directories
- [ ] Test commands work (`PYTHONPATH=src uv run pytest tests/ -v`)
- [ ] No broken internal links
- [ ] All code blocks have proper syntax highlighting
- [ ] Markdown renders correctly in GitHub preview

---

## Appendix: Actual Skill List

For reference, here are all 26 SKILL.md files:

**Top-level (22):**
1. citation-management
2. code-documentation
3. document-skills (container)
4. generate-image
5. hypothesis-generation
6. latex-posters
7. literature-review
8. markdown-to-pdf
9. markitdown
10. paper-2-web
11. peer-review
12. plotting-libraries
13. pptx-posters
14. reproducible-research
15. research-lookup
16. scholar-evaluation
17. scientific-schematics
18. scientific-slides
19. scientific-writing
20. statistical-analysis
21. venue-templates
22. visual-design

**Sub-skills (4):**
- document-skills/docx
- document-skills/pdf
- document-skills/pptx
- document-skills/xlsx
