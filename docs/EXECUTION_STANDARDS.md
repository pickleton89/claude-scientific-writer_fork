# Execution Standards for Scientific Writing

> Workflow standards for scientific writing tasks
> Version: 1.0.0 | Created: 2025-12-30

These standards apply when executing scientific writing tasks (papers, literature reviews, posters, presentations). They complement the individual skill definitions in `skills/` by providing cross-cutting workflow guidance.

---

## File Organization

All scientific writing outputs follow a standardized directory structure:

```
writing_outputs/
└── YYYYMMDD_HHMMSS_<brief_description>/
    ├── progress.md          # Task progress log
    ├── SUMMARY.md           # Final deliverables summary
    ├── PEER_REVIEW.md       # Self-review before delivery
    ├── drafts/              # Versioned manuscript drafts
    │   ├── v1_draft.tex
    │   ├── v2_draft.tex
    │   └── revision_notes.md
    ├── references/          # Bibliography files
    │   └── references.bib
    ├── figures/             # Generated figures and images
    │   ├── figure_01.png
    │   └── graphical_abstract.png
    ├── data/                # Data files (csv, json, xlsx)
    ├── sources/             # Context materials (.md, .docx, .pdf)
    └── final/               # Final compiled outputs
        ├── manuscript.pdf
        └── manuscript.tex
```

### File Routing Rules

When files are provided in the `data/` folder:

| File Type | Destination | Action |
|-----------|-------------|--------|
| `.tex` files | `drafts/` | **EDITING MODE** - modify existing manuscript |
| Images (`.png`, `.jpg`, `.svg`) | `figures/` | Reference in manuscript |
| Data files (`.csv`, `.json`, `.xlsx`) | `data/` | Use for analysis/figures |
| Other files (`.md`, `.docx`, `.pdf`) | `sources/` | Context materials |

---

## Progress Logging

Track progress in `progress.md` using timestamped entries:

```markdown
# Progress Log

## Task: [Brief description]
Started: YYYY-MM-DD HH:MM:SS

### Progress

[HH:MM:SS] STARTED: Task initialization
[HH:MM:SS] RESEARCH: Found N papers on [topic]
[HH:MM:SS] COMPLETED: Introduction section - 450 words, 8 citations
[HH:MM:SS] GENERATED: Graphical abstract for paper summary
[HH:MM:SS] VERIFIED: Smith2024 citation confirmed via DOI
[HH:MM:SS] COMPLETED: Methods section - 600 words, 5 citations
...
```

### Log Entry Types

| Prefix | Usage |
|--------|-------|
| `STARTED` | Task or section initiation |
| `RESEARCH` | Literature search results |
| `COMPLETED` | Section finished with word count and citations |
| `GENERATED` | Figure or schematic created |
| `VERIFIED` | Citation or fact verified |
| `COMPILED` | LaTeX compilation step |
| `REVIEWED` | Quality review checkpoint |

---

## Version Management

**Never overwrite previous versions.** Always increment version numbers when editing:

| Version | File Name | Purpose |
|---------|-----------|---------|
| Initial | `v1_draft.tex` | First complete draft |
| Revision 1 | `v2_draft.tex` | After first round of edits |
| Revision 2 | `v3_draft.tex` | After second round of edits |
| ... | `vN_draft.tex` | Continue incrementing |

### Revision Notes

Document changes in `drafts/revision_notes.md`:

```markdown
# Revision Notes

## v1 → v2 (YYYY-MM-DD)
- Expanded Methods section with protocol details
- Added 3 new citations to Discussion
- Fixed Figure 2 axis labels

## v2 → v3 (YYYY-MM-DD)
- Addressed reviewer comment on statistical methods
- Added supplementary table S1
- Shortened Introduction by 150 words
```

---

## PDF Review Workflow

After compiling any PDF, review it visually before delivery:

### Step 1: Convert to Images

**Never read PDFs directly.** Convert to images for visual inspection:

```bash
python scripts/pdf_to_images.py document.pdf review/page --dpi 150
```

### Step 2: Inspect Each Page

Check for common issues:

| Issue | What to Look For |
|-------|------------------|
| Text overlaps | Text running into margins or figures |
| Figure placement | Figures appearing on wrong pages |
| Table breaks | Tables split awkwardly across pages |
| Margins | Content too close to edges |
| Page breaks | Orphan/widow lines, section breaks |
| Captions | Caption spacing, alignment with figures |
| Bibliography | Consistent formatting, no broken entries |

### Step 3: Fix and Recompile

If issues found:
1. Edit the `.tex` source
2. Recompile: `pdflatex → bibtex → pdflatex → pdflatex`
3. Re-convert to images and inspect
4. Maximum 3 iterations (then deliver with notes)

### Step 4: Clean Up

```bash
rm -rf review/
```

---

## Quality Assurance

Before marking any scientific writing task complete:

### Completion Checklist

- [ ] All files created and properly formatted
- [ ] Version numbers incremented (if editing existing work)
- [ ] 100% of citations are REAL papers (verified via research-lookup)
- [ ] All citation metadata verified with DOIs
- [ ] Graphical abstract generated (for papers/reviews)
- [ ] Minimum figure count met (see scientific-writing skill)
- [ ] Figures properly integrated with captions and references
- [ ] `progress.md` complete with all actions logged
- [ ] `SUMMARY.md` lists all deliverables
- [ ] `PEER_REVIEW.md` completed (use peer-review skill)
- [ ] PDF formatting review passed (no visual issues)

### SUMMARY.md Template

```markdown
# Deliverables Summary

## Task
[Brief description of what was requested]

## Completed: YYYY-MM-DD HH:MM:SS

## Files Produced

| File | Description |
|------|-------------|
| `final/manuscript.pdf` | Compiled manuscript |
| `final/manuscript.tex` | LaTeX source |
| `references/references.bib` | Bibliography (N entries) |
| `figures/graphical_abstract.png` | Graphical abstract |
| `figures/figure_01.png` | [Description] |

## Metrics

- Total words: N
- Total citations: N
- Total figures: N

## Usage

To recompile:
\`\`\`bash
cd final/
pdflatex manuscript.tex
bibtex manuscript
pdflatex manuscript.tex
pdflatex manuscript.tex
\`\`\`
```

---

## Cross-References

| Document | Purpose |
|----------|---------|
| `CLAUDE.md` | Agent behavior policies (completion, decision making) |
| `skills/SKILL_ROUTER.md` | Skill selection decision trees |
| `skills/scientific-writing/SKILL.md` | Manuscript writing guidance |
| `skills/citation-management/SKILL.md` | Citation verification |
| `skills/peer-review/SKILL.md` | Self-review framework |

---

*These standards ensure consistent, high-quality scientific writing outputs across all document types.*
