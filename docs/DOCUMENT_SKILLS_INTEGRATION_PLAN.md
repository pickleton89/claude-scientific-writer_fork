# Document Skills Integration Plan

> Implementation plan for properly integrating `document-skills` into the skill system workflow
> Created: 2025-12-29

---

## Problem Statement

The `document-skills` container skill (docx, pdf, pptx, xlsx) is not fully integrated into the overall skill routing system. While individual skills like `scientific-slides` correctly reference `document-skills/pptx`, there are gaps in:

1. The central SKILL_ROUTER decision trees
2. Cross-references from content-creation skills to output-format skills
3. The complete workflow from research → writing → formatted output

---

## Tasks Overview

| # | Task | File(s) to Modify | Priority |
|---|------|-------------------|----------|
| 1 | Add "Document Export" decision tree to SKILL_ROUTER | `skills/SKILL_ROUTER.md` | High |
| 2 | Add document-skills to Quick Reference Matrix | `skills/SKILL_ROUTER.md` | High |
| 3 | Add cross-reference in scientific-writing | `skills/scientific-writing/SKILL.md` | Medium |
| 4 | Add cross-reference in venue-templates | `skills/venue-templates/SKILL.md` | Medium |
| 5 | Update Multi-Skill Workflows in SKILL_ROUTER | `skills/SKILL_ROUTER.md` | Medium |
| 6 | Verify all document-skills sub-skills have proper cross-refs | `skills/document-skills/*/SKILL.md` | Low |

---

## Task 1: Add "Document Export" Decision Tree to SKILL_ROUTER

**File:** `skills/SKILL_ROUTER.md`

**Location:** After "Decision Tree: Conversion & Transformation" section (around line 277)

**Content to Add:**

~~~markdown
---

## Decision Tree: Document Export & Output Formats

```
User has content ready for final document format
│
├─ Need Word document (.docx)?
│  │
│  └─ YES → document-skills/docx
│            │
│            ├─ From scientific-writing content
│            │  └─ Export manuscript to Word format
│            │
│            ├─ With venue-templates formatting
│            │  └─ Apply journal-specific styles
│            │
│            ├─ Need tracked changes/redlining?
│            │  └─ Use OOXML redline workflow
│            │
│            └─ Edit existing DOCX?
│               └─ Raw XML or docx-js manipulation
│
├─ Need PDF output?
│  │
│  ├─ Branded/templated PDF from Markdown?
│  │  └─ YES → markdown-to-pdf
│  │            └─ Uses Oligon templates & brand colors
│  │
│  └─ PDF manipulation (forms, merge, extract)?
│     └─ YES → document-skills/pdf
│              │
│              ├─ Fill PDF forms
│              ├─ Extract text/tables
│              ├─ Merge/split PDFs
│              └─ Metadata handling
│
├─ Need PowerPoint (.pptx)?
│  │
│  ├─ Creating presentation from scratch?
│  │  └─ YES → scientific-slides → document-skills/pptx
│  │            │
│  │            ├─ Plan content with scientific-slides
│  │            └─ Build PPTX with document-skills/pptx
│  │
│  └─ Editing existing presentation?
│     └─ YES → document-skills/pptx
│              │
│              ├─ OOXML manipulation
│              ├─ Template-based creation
│              └─ html2pptx workflow
│
└─ Need Excel (.xlsx)?
   │
   └─ YES → document-skills/xlsx
            │
            ├─ Create spreadsheets with formulas
            ├─ Data analysis and manipulation
            ├─ Chart generation
            └─ Recalculate formulas (recalc.py)
```

### Document Export: Quick Selection

| Content Source | Output Format | Primary Skill | Supporting Skill |
|----------------|---------------|---------------|------------------|
| Manuscript text | Word (.docx) | document-skills/docx | venue-templates |
| Markdown content | Branded PDF | markdown-to-pdf | visual-design |
| Existing PDF | Modified PDF | document-skills/pdf | — |
| Slide plan | PowerPoint | document-skills/pptx | scientific-slides |
| Data/analysis | Excel | document-skills/xlsx | — |
~~~

---

## Task 2: Add Document-Skills to Quick Reference Matrix

**File:** `skills/SKILL_ROUTER.md`

**Location:** In the "Quick Reference Matrix" section, add new category after "Conversion"

**Content to Add:**

```markdown
| **Document Output** | | |
| Export to Word | document-skills/docx | + venue-templates (for journal styles) |
| Create branded PDF | markdown-to-pdf | + visual-design |
| PDF forms/manipulation | document-skills/pdf | — |
| Create/edit PowerPoint | document-skills/pptx | + scientific-slides |
| Create/edit Excel | document-skills/xlsx | — |
```

---

## Task 3: Add Cross-Reference in scientific-writing

**File:** `skills/scientific-writing/SKILL.md`

**Location:** In the `## Integration with Other Skills` section (around line 279-291)

**Current Content (lines 283-289):**
```markdown
**Related skills:**
- **research-lookup**: Literature search and citation finding
- **scientific-schematics/generate-image**: Figure creation and visual elements
- **venue-templates**: Journal-specific styles and formatting requirements
- **statistical-analysis**: Test selection, statistical reporting in Methods, interpreting results
- **reproducible-research**: Data Availability statements, environment documentation, FAIR compliance
- **code-documentation**: Documenting analysis code referenced in Methods sections
```

**Updated Content:**
```markdown
**Related skills:**
- **research-lookup**: Literature search and citation finding
- **scientific-schematics/generate-image**: Figure creation and visual elements
- **venue-templates**: Journal-specific styles and formatting requirements
- **statistical-analysis**: Test selection, statistical reporting in Methods, interpreting results
- **reproducible-research**: Data Availability statements, environment documentation, FAIR compliance
- **code-documentation**: Documenting analysis code referenced in Methods sections
- **document-skills/docx**: Export manuscript to Word format for journal submission
- **document-skills/pdf**: PDF output for preprints or supplementary materials

**Output Workflow**: After drafting content with this skill, use `venue-templates` for journal-specific formatting, then export via `document-skills/docx` for Word submissions or `markdown-to-pdf` for branded PDF output.
```

---

## Task 4: Add Cross-Reference in venue-templates

**File:** `skills/venue-templates/SKILL.md`

**Location:** In the `## Integration with Other Skills` section (around line 319-344)

**Action:** Add a new subsection after the existing subsections (Scientific Writing, Literature Review, Peer Review, LaTeX Posters)

**Content to Add:**
```markdown
### Document Export
- Use **document-skills/docx** to export formatted manuscript to Word for journal submission systems
- Use **document-skills/pdf** to generate PDF version for preprint servers (arXiv, bioRxiv)
- Use **markdown-to-pdf** for branded PDF output with Oligon templates

**Export Workflow**: After applying venue-specific formatting with this skill, use `document-skills/docx` to generate the final Word document for submission.
```

---

## Task 5: Update Multi-Skill Workflows in SKILL_ROUTER

**File:** `skills/SKILL_ROUTER.md`

**Location:** In "Multi-Skill Workflows" section, update Workflow 1

**Current Content:**
~~~markdown
### Workflow 1: Complete Research Paper

```
1. literature-review      → Systematic search & screening
2. citation-management    → Organize references
3. hypothesis-generation  → Formulate testable hypotheses
4. statistical-analysis   → Plan analysis approach
5. scientific-writing     → Draft manuscript
6. venue-templates        → Format for target journal
7. plotting-libraries     → Create data figures
8. scientific-schematics  → Create method diagrams
9. visual-design          → Polish all visuals
10. peer-review           → Self-review before submission
```
~~~

**Updated Content:**
~~~markdown
### Workflow 1: Complete Research Paper

```
1. literature-review      → Systematic search & screening
2. citation-management    → Organize references
3. hypothesis-generation  → Formulate testable hypotheses
4. statistical-analysis   → Plan analysis approach
5. scientific-writing     → Draft manuscript
6. venue-templates        → Format for target journal
7. plotting-libraries     → Create data figures
8. scientific-schematics  → Create method diagrams
9. visual-design          → Polish all visuals
10. peer-review           → Self-review before submission
11. document-skills/docx  → Export to Word for submission
```
~~~

**Also Add New Workflow:**
~~~markdown
### Workflow 5: PDF Report from Markdown

```
1. scientific-writing     → Draft content in Markdown
2. visual-design          → Establish design specifications
3. plotting-libraries     → Create data figures
4. markdown-to-pdf        → Apply branded template
   OR
   document-skills/pdf    → Advanced PDF manipulation
```

### Workflow 6: Report to Presentation Conversion

```
1. markitdown             → Extract PDF/DOCX to Markdown
2. visual-design          → Load brand specifications
3. scientific-slides      → Plan slide structure
4. document-skills/pptx   → Build PowerPoint file
5. visual-design          → Validate accessibility
```
~~~

---

## Task 6: Verify Document-Skills Sub-Skill Cross-References

**Files to Check:**
- `skills/document-skills/docx/SKILL.md`
- `skills/document-skills/pdf/SKILL.md`
- `skills/document-skills/pptx/SKILL.md`
- `skills/document-skills/xlsx/SKILL.md`

**Verification Checklist:**

| Sub-Skill | Has cross_references section? | References SKILL_ROUTER? | References related content skills? |
|-----------|------------------------------|--------------------------|-----------------------------------|
| docx | Check | Should reference | scientific-writing, venue-templates |
| pdf | Check | Should reference | markdown-to-pdf, markitdown |
| pptx | Check | Should reference | scientific-slides |
| xlsx | Check | Should reference | statistical-analysis |

**Action:** For each sub-skill, ensure the cross_references section includes:
```markdown
**Skill Selection:** See `SKILL_ROUTER.md` for decision trees when multiple skills may apply.
```

---

## Implementation Order

```
Phase 1: Core Routing (High Priority)
├── Task 1: Add Document Export decision tree
└── Task 2: Add to Quick Reference Matrix

Phase 2: Cross-References (Medium Priority)
├── Task 3: Update scientific-writing
├── Task 4: Update venue-templates
└── Task 5: Update Multi-Skill Workflows

Phase 3: Verification (Low Priority)
└── Task 6: Verify sub-skill cross-references
```

---

## Validation

After implementation, verify:

1. **SKILL_ROUTER completeness:**
   ```bash
   grep -c "document-skills" skills/SKILL_ROUTER.md
   # Should return 10+ matches
   ```

2. **Cross-reference consistency:**
   ```bash
   grep -l "document-skills" skills/*/SKILL.md
   # Should include: scientific-writing, venue-templates, markdown-to-pdf, scientific-slides
   ```

3. **Decision tree coverage:**
   - [ ] All 4 document formats (docx, pdf, pptx, xlsx) appear in SKILL_ROUTER
   - [ ] Each format has clear routing logic
   - [ ] Export workflows connect content skills to output skills

---

## Success Criteria

- [ ] SKILL_ROUTER has "Document Export" decision tree
- [ ] Quick Reference Matrix includes document output category
- [ ] scientific-writing references document-skills/docx
- [ ] venue-templates references document-skills/docx
- [ ] Multi-Skill Workflows include document export steps
- [ ] All document-skills sub-skills reference SKILL_ROUTER

---

*Plan created: 2025-12-29*
*Status: Ready for implementation*
