# Claude Scientific Writer Fork

![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)
![Tests](https://img.shields.io/badge/tests-24%20passing-brightgreen.svg)
![Skills](https://img.shields.io/badge/skills-22%20active-purple.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

A customized fork of [K-Dense-AI/claude-scientific-writer](https://github.com/K-Dense-AI/claude-scientific-writer) (v2.10.0), adapted as a personal skill repository for scientific and biotech writing workflows.

> **Note**: This is a working project for personal/organizational use. It is **not intended to be installed as a Claude Code plugin**. For the production-ready version, use the [upstream repository](https://github.com/K-Dense-AI/claude-scientific-writer).

## What This Is

This fork serves as a **customized skill repository** built on top of the Claude Scientific Writer foundation. Rather than using it as an installable plugin, I use this as a working project directory where:

- Skills are customized for my specific scientific writing needs
- Oligon brand standards are integrated for document generation
- Clinical/business-focused skills have been removed to maintain a pure scientific focus
- New visualization and design skills have been added
- **Template-based document generation** with 12 document types and branded PDF output

### Key Customizations

| Change | Description |
|--------|-------------|
| **Skills Streamlined** | Removed 5 clinical/business skills; added 4 new skills (22 total) |
| **Brand Integration** | Added Oligon brand standards for PDF generation and visual identity |
| **Visualization Skills** | Added `plotting-libraries` and `visual-design` skills |
| **Document Templates** | Added `markdown-to-pdf` skill with 12 document types and schema validation |
| **Single Source of Truth** | Consolidated skills to one location (`skills/`) instead of three |

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

### Integration Testing
24 integration tests covering the full template-to-PDF workflow.

## Skills (22 top-level + 4 document sub-skills)

| Category | Skills |
|----------|--------|
| **Writing** | `scientific-writing`, `literature-review`, `hypothesis-generation` |
| **Presentations** | `scientific-slides`, `latex-posters`, `pptx-posters` |
| **Research** | `research-lookup`, `citation-management`, `peer-review`, `scholar-evaluation` |
| **Visuals** | `scientific-schematics`, `generate-image`, `plotting-libraries`, `visual-design` |
| **Documents** | `markdown-to-pdf`, `markitdown`, `venue-templates`, `document-skills/*` (docx, pdf, pptx, xlsx) |
| **Analysis** | `statistical-analysis`, `reproducible-research`, `code-documentation` |
| **Conversion** | `paper-2-web` |

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

## Document Template System

The `markdown-to-pdf` skill provides a template-based document generation system with **12 document types** across 4 categories:

| Category | Document Types |
|----------|----------------|
| **Scientific** | `analysis-report`, `literature-review`, `data-report` |
| **Project Management** | `meeting-notes`, `project-status`, `task-list`, `phase-plan` |
| **Development** | `technical-spec`, `agent-definition`, `readme` |
| **Reference** | `standards-guide`, `method-guide` |

### Usage

```bash
# List available templates
/list-templates

# Create a new document from template
/new-doc analysis-report

# Convert markdown to branded PDF
/doc-to-pdf my-document.md
```

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

## Project Structure

```
claude-scientific-writer_fork/
├── skills/                  # Canonical skill definitions (22 + 4 sub-skills)
│   ├── SKILL_ROUTER.md      # Multi-skill task routing
│   ├── SKILL_TEMPLATE.md    # Standard skill template
│   ├── SKILL_TESTS.md       # Validation test scenarios
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

## Getting Started

### Prerequisites

- Python 3.10+
- [uv](https://github.com/astral-sh/uv) package manager
- LaTeX distribution (MacTeX, TeX Live, or MiKTeX)

### Installation

```bash
# Clone the repository
git clone <your-fork-url>
cd claude-scientific-writer_fork

# Install dependencies with uv
uv sync
```

### Usage

This project is designed to be used as a **working directory with Claude Code**, not as an installed plugin.

1. Open this directory in Claude Code
2. The skills in `skills/` provide specialized capabilities for scientific writing
3. Invoke skills as needed during your writing workflow

For programmatic usage outside Claude Code:

```bash
# CLI interface
uv run scientific-writer

# Or use the Python API
from scientific_writer import generate_paper
```

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

## Development

### Package Manager

```bash
uv sync                    # Install dependencies
uv run python <script>     # Run scripts in correct environment
uv run ruff check .        # Lint code
uv run ty                  # Type check
```

### Editing Skills

All skill edits happen in the `skills/` directory. Each skill has:
- `SKILL.md` - Skill definition and capabilities
- `scripts/` - Optional automation scripts
- `references/` - Optional reference materials
- `assets/` - Optional static assets

### Documentation

| Document | Description |
|----------|-------------|
| `CLAUDE.md` | Development documentation and workflow |
| `CHANGELOG.md` | Fork customization history |
| `.claude/WRITER.md` | Scientific writing agent instructions |
| `skills/SKILL_ROUTER.md` | Multi-skill task routing decision trees |
| `skills/SKILL_TEMPLATE.md` | Standard template for new skills |
| `docs/archived/INTEGRATION_ANALYSIS.md` | Completed template-project merge |

## Upstream Reference

This fork is based on [K-Dense-AI/claude-scientific-writer](https://github.com/K-Dense-AI/claude-scientific-writer) v2.10.0.

For the production-ready Claude Scientific Writer plugin, visit the upstream repository.

## License

MIT (inherited from upstream)
