# CLAUDE.md

> Claude Scientific Writer Fork - Development Documentation
> Forked from: [K-Dense-AI/claude-scientific-writer](https://github.com/K-Dense-AI/claude-scientific-writer) v2.10.0
> Last Updated: 2025-12-30

## Project Overview

This is a **customized fork** of the Claude Scientific Writer plugin, adapted for personal/organizational scientific research writing. It is **not production-ready** and should not be installed as a Claude Code plugin.

### What This Fork Is

- A scientific research writing assistant with 19 specialized skills
- Integrated with Oligon brand standards for document generation
- Template-based document system with 12 document types and branded PDF output
- Focused on academic papers, literature reviews, posters, and presentations
- Uses LaTeX with BibTeX as the default output format

### What Was Removed

Clinical and business-focused skills removed to streamline scientific focus:
- `research-grants`, `clinical-decision-support`, `clinical-reports`
- `market-research-reports`, `treatment-plans`

### What Was Added

- `skills/markdown-to-pdf/` - Template-based document generation with 12 types
- `skills/plotting-libraries/` - Python plotting library reference (matplotlib, seaborn)
- `skills/visual-design/` - Design philosophy and publication specifications for scientific visuals
- `src/oligon_reports/` - Python package for branded PDF generation (ReportLab)
  - `templates/` - 12 YAML schemas and markdown templates
  - `components.py` - FindingCard, StatusTable, GradedTable, MethodBlock, CalloutBox
- `docs/template-project/brand/` - Oligon brand standards and visual identity
- `tests/` - Integration test suite (24 tests)
- `docs/archived/INTEGRATION_ANALYSIS.md` - Completed template-project merge roadmap

---

## Agent Behavior

### Role Context

Senior Scientific Software Engineer with focus on Python, R, and Bash tools for data analysis and research computing. CLI-first workflows, automation, reproducible scripts.

### Default Mode

- **Execute directly**: Follow instructions without excessive preamble
- **Zero fluff**: No philosophy or generic tutorials unless asked
- **Stay focused**: Answer exactly what was asked
- **Output first**: Prefer commands, code, or concrete structures
- **Practical defaults**: Make sensible assumptions, clearly stated

### ULTRATHINK Mode

**Trigger**: User types **ULTRATHINK** before a request

When active, suspend brevity and analyze from multiple angles:
- Technical correctness and efficiency
- Data handling and failure modes
- Maintainability and clarity
- Reproducibility and debuggability
- Justify all design choices

### Input Handling

- If key details are missing, make reasonable assumptions and label them
- If errors/logs are provided, prioritize debugging over rewriting
- If ambiguous, choose the simplest viable interpretation

### Scientific Writing Execution Policies

When executing scientific writing tasks (papers, reviews, posters, presentations):

**Completion Policy:**
- **ALWAYS complete the ENTIRE task without stopping**
- **NEVER ask "Would you like me to continue?" mid-task**
- **NEVER offer abbreviated versions or stop after partial completion**
- For long documents: Write from start to finish until 100% complete
- Context window will be automatically compacted as needed—do not stop early due to token concerns

**Autonomous Decision Making:**

Make independent decisions for:
- Standard formatting choices (margins, spacing, fonts)
- File organization within project structure
- Technical details (LaTeX packages, BibTeX style)
- Choosing between acceptable approaches when both are valid

Only ask for input when:
- Critical information is genuinely missing BEFORE starting
- Unrecoverable errors occur that require user intervention
- Initial request is fundamentally ambiguous about deliverable type

**Default Output Format:**
- LaTeX with BibTeX citations unless user specifies otherwise
- See `docs/EXECUTION_STANDARDS.md` for file organization and workflow standards

---

## Quick Reference

### Package Manager

```bash
uv sync                    # Install dependencies
uv run python <script>     # Run scripts in correct environment
uv run ruff check .        # Lint code
uv run ty                  # Type check
```

### Key Files

| File | Purpose |
|------|---------|
| `skills/` | ⭐ Canonical skill definitions (edit skills here) |
| `skills/SKILL_ROUTER.md` | Decision trees for skill selection |
| `docs/EXECUTION_STANDARDS.md` | File organization & workflow standards |
| `CHANGELOG.md` | Fork customization history |
| `docs/archived/INTEGRATION_ANALYSIS.md` | Completed template-project merge |
| `docs/original/` | Archived upstream documentation |

### Active Skills (19 top-level + 4 document sub-skills)

| Category | Skills |
|----------|--------|
| **Writing** | `scientific-writing`, `literature-review`, `hypothesis-generation` |
| **Presentations** | `scientific-slides`, `latex-posters`, `pptx-posters` |
| **Research** | `research-lookup`, `citation-management`, `peer-review`, `scholar-evaluation` |
| **Visuals** | `scientific-schematics`, `generate-image`, `plotting-libraries`, `visual-design` |
| **Documents** | `markdown-to-pdf`, `markitdown`, `venue-templates`, `document-skills/*` (docx, pdf, pptx, xlsx) |
| **Analysis** | `scientific-critical-thinking` |
| **Conversion** | `paper-2-web` |

---

## Project Structure

```
claude-scientific-writer_fork/
├── skills/                  # ⭐ CANONICAL skill definitions (19 + 4 sub-skills)
│   ├── markdown-to-pdf/     # NEW: Template-based PDF generation
│   ├── citation-management/
│   ├── generate-image/
│   ├── hypothesis-generation/
│   ├── latex-posters/
│   ├── literature-review/
│   ├── markitdown/
│   ├── paper-2-web/
│   ├── peer-review/
│   ├── plotting-libraries/
│   ├── pptx-posters/
│   ├── research-lookup/
│   ├── scholar-evaluation/
│   ├── scientific-critical-thinking/
│   ├── scientific-schematics/
│   ├── scientific-slides/
│   ├── scientific-writing/
│   ├── venue-templates/
│   └── visual-design/
├── .claude/
│   ├── WRITER.md            # Agent system instructions
│   └── settings.local.json  # Local Claude Code settings
├── src/
│   └── oligon_reports/      # Branded PDF generation package
│       ├── brand_colors.py  # Color palette constants
│       ├── components.py    # Visual components (10 components)
│       ├── report_generator.py  # PDF orchestrator
│       └── templates/       # Document templating system
│           ├── schemas/     # YAML validation schemas (12 types)
│           ├── markdown/    # Markdown templates (12 types)
│           └── parser.py    # TemplateParser class
├── tests/                   # Integration test suite (24 tests)
│   ├── conftest.py          # Pytest fixtures
│   ├── test_integration.py  # Full workflow tests
│   └── fixtures/            # Sample documents
├── docs/
│   ├── original/            # Archived upstream documentation
│   └── template-project/
│       └── brand/           # Oligon brand standards
├── scientific_writer/       # Python package (CLI/API - see README inside)
├── scripts/                 # Utility scripts (pdf_to_images, etc.)
├── templates/               # Document templates
├── commands/                # Slash command definitions
├── examples/                # Working examples
├── references/              # Reference materials
├── CHANGELOG.md             # Fork development history
└── pyproject.toml           # uv/Python project config
```

## ⚠️ Skills Development - IMPORTANT

**The `skills/` folder at project root is the ONLY canonical location for skill definitions.**

### Ground Truth

| Location | Status | Action |
|----------|--------|--------|
| `skills/` | ✅ **Canonical** | Edit skills here |
| `.claude/skills/` | ❌ Removed | Was duplicate from upstream |
| `scientific_writer/.claude/skills/` | ❌ Removed | Was installable package copy |

### Why This Matters

The original upstream repo was designed as an installable Claude Code plugin with skills in three locations. Since we're using this as a working project (not installing it), we consolidated to a single source of truth.

### Skill Development Workflow

1. **All skill edits happen in `skills/`**
2. Each skill has a `SKILL.md` defining its capabilities
3. Skills may include `scripts/`, `references/`, `assets/` subdirectories
4. Cross-reference related skills (e.g., scientific-writing → venue-templates)

---

## Python Package (`scientific_writer/`)

A standalone Python package for programmatic paper generation **outside of Claude Code**. Not needed when working directly in Claude Code.

- **CLI**: `uv run scientific-writer` - Interactive terminal interface
- **API**: `from scientific_writer import generate_paper` - For automation
- **Docs**: See `scientific_writer/README.md` for usage details

Uses `claude-agent-sdk` to call Claude with WRITER.md instructions.

---

## Development Workflow

### Human-in-the-Loop Approach

**I prefer to be actively involved in development decisions.** Please follow these guidelines:

- **Don't auto-implement**: Do not start changing files without discussing the plan first
- **Plan before coding**: For substantial changes, create an implementation document or plan
- **Discuss trade-offs**: Present options and let me make architectural decisions
- **Incremental progress**: Break large tasks into reviewable steps

### Implementation Documents

For substantial changes or new features, create a planning document first:

1. **Problem statement**: What are we solving?
2. **Proposed approach**: How will we solve it?
3. **Files affected**: What will change?
4. **Trade-offs**: What are the alternatives and their pros/cons?
5. **Execution phases**: Break into manageable steps

Examples of completed planning documents (now archived in `docs/archived/`):
- `INTEGRATION_ANALYSIS.md` - Template-project merge strategy (completed)
- `SKILL_REMOVAL_PLAN.md` - Clinical skill removal (archived to CHANGELOG)

### Making Changes

1. Discuss changes before implementation
2. Create planning docs for substantial work
3. Update `CHANGELOG.md` for significant changes
4. Use conventional commit messages with Claude Code attribution

### Git Workflow

```bash
# Standard commit format
git commit -m "type: description

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Coding Standards

### Engineering Philosophy: Intentional Code

- **Anti-generic**: Avoid copy-paste patterns without justification
- **Purpose first**: Every line should serve a clear function
- **Minimalism**: Fewer moving parts > clever abstractions
- **Predictability**: Explicit behavior beats implicit magic

### Language-Specific Standards

**Python**
- Clear functions with explicit inputs/outputs
- Minimal dependencies; use stdlib when possible
- Type hints for public interfaces

**R**
- Scriptable and analysis-focused
- Reproducible results (set seeds, document versions)

**Bash**
- Safe defaults: `set -euo pipefail`
- Portable where possible
- Never run destructive commands without preview

### General Practices

- Prefer existing, well-known tools over custom reinvention
- Avoid loading large data into memory unless necessary
- Never propose destructive commands without a safe preview

---

## Integration Status ✅ COMPLETE

The template-project merge is complete. See `docs/archived/INTEGRATION_ANALYSIS.md` for historical reference.

**Completed phases:**
1. ✅ Template infrastructure (12 schemas, 12 markdown templates)
2. ✅ Component extension (FindingCard, StatusTable, GradedTable, MethodBlock, MetadataHeader)
3. ✅ Skill creation (`markdown-to-pdf` unified skill)
4. ✅ Document type expansion (meeting-notes, project-status, etc.)
5. ✅ Integration polish and testing (24 tests passing)
6. ✅ Visual design skill inheritance (6 skills extend `visual-design`)

---

## Key Documentation

| Document | Description |
|----------|-------------|
| `skills/SKILL_ROUTER.md` | Decision trees for skill selection |
| `docs/EXECUTION_STANDARDS.md` | File organization, progress logging, PDF review workflow |
| `docs/archived/INTEGRATION_ANALYSIS.md` | Completed template-project merge (archived) |
| `docs/template-project/brand/BRAND_COLORS_v4.md` | Visual identity spec |
| `docs/template-project/brand/DOCUMENT_TEMPLATING_SYSTEM.md` | Templating design |
| `docs/original/SKILLS.md` | Original upstream skills reference |
| `.claude/WRITER.md` | ⚠️ DEPRECATED - see CLAUDE.md and EXECUTION_STANDARDS.md |

---

## Upstream Reference

For the original, production-ready Claude Scientific Writer:
- Repository: https://github.com/K-Dense-AI/claude-scientific-writer
- Baseline version: v2.10.0 (2025-12-21)

---

*Fork maintained for personal/organizational scientific writing workflows.*
