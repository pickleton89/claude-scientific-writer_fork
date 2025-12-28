# Skills Library Determinism Audit

> Comprehensive evaluation of deterministic patterns across the scientific writing skills library
> Generated: 2025-12-27
> Overall Score: 7.5/10 (potential to reach 8.5/10 with targeted improvements)

---

## Executive Summary

This audit evaluates all 22 skills in the library across 8 determinism dimensions. The library shows **strong foundations** with several exemplary skills that demonstrate best practices. Newer bioinformatics-focused skills (Phase 1-4 implementations) show more deterministic patterns than inherited upstream skills.

### Scoring Rubric

| Score | Rating | Criteria |
|-------|--------|----------|
| 9.0-10 | Highly Deterministic | XML tags, decision tables, numbered workflows, explicit templates, quantitative thresholds |
| 7.0-8.9 | Strong Determinism | Clear workflows, good decision frameworks, some quantitative criteria |
| 5.5-6.9 | Moderate Determinism | Basic structure but missing decision trees, implicit criteria |
| <5.5 | Needs Improvement | Prose-heavy, ambiguous routing, no explicit success criteria |

### Key Findings

**Strengths:**
- Three exemplary skills (statistical-analysis, reproducible-research, code-documentation) provide templates for improvement
- Newer bioinformatics additions use consistent XML semantic tags
- Cross-references between skills are well-documented

**Gaps:**
- Inconsistent structure patterns across older skills
- Missing quantitative thresholds for qualitative terms
- Undefined workflow transition criteria
- Incomplete skill orchestration rules

---

## Tier 1: Highly Deterministic (8.5-10)

### statistical-analysis (9.5/10)

**Strengths:**
- XML semantic tags (`<decision_framework>`, `<workflow>`, `<scope>`)
- Decision matrices with explicit condition→test mappings
- Numbered 6-step workflow with transition criteria
- 10 numbered pitfalls with explicit solutions
- Reporting templates with placeholder strings
- Quantitative thresholds (e.g., "n<30", "p<0.05", "≥80% power")

**Improvements Needed:**
- Add effect size interpretation thresholds (small/medium/large by test type)

---

### reproducible-research (9.0/10)

**Strengths:**
- XML semantic tags throughout
- FAIR principles decomposed into actionable items
- 3-level environment specification hierarchy with explicit criteria
- 10 numbered reproducibility errors with anti-patterns
- Decision matrix mapping data types to repositories
- Template strings for Data Availability statements

**Improvements Needed:**
- Add workflow transition checklist between specification levels
- Quantify "minimal" vs "standard" vs "complete" with concrete metrics

---

### code-documentation (8.5/10)

**Strengths:**
- Decision flowchart for style selection (ASCII art)
- Style comparison table with 5 dimensions
- README template with explicit sections
- Docstring patterns by language (Python, R, Bash)
- Anti-patterns section with explicit "don't do this" examples

**Improvements Needed:**
- Add quantitative thresholds for "well-documented" (e.g., "≥80% public functions have docstrings")
- Include workflow for adding documentation to existing code

---

### plotting-libraries (8.5/10)

**Strengths:**
- Python vs R decision matrix with 6 evaluation factors
- Tool routing within each ecosystem (seaborn vs matplotlib, ggplot2 vs Bioconductor)
- Export specifications with exact values (300 DPI, vector formats)
- Bioinformatics visualization patterns with tool recommendations
- Success criteria checklist (7 items)

**Improvements Needed:**
- Add decision tree for choosing between 6 bioinformatics plot types
- Quantify "publication quality" criteria

---

## Tier 2: Strong Determinism (7.0-8.4)

### scientific-writing (8.0/10)

**Strengths:**
- Two-stage writing process (outline → prose) explicitly defined
- IMRAD section guidance with specific sentence counts
- Figure requirements table by document type
- Pre-submission diagnostic tests (3 specific checks)
- Cross-references to related skills

**Improvements Needed:**
- Add decision tree for section order based on paper type
- Quantify "comprehensive" literature review (target citation count)
- Include word count targets per section

---

### peer-review (7.5/10)

**Strengths:**
- 7-stage systematic workflow
- Presentation review protocol (separate from manuscript)
- Reviewer checklist with 12 specific items
- Response template with required sections

**Improvements Needed:**
- Add severity scoring rubric for issues (critical/major/minor)
- Decision tree for when to recommend accept/revise/reject
- Time allocation per review stage

---

### research-lookup (7.5/10)

**Strengths:**
- Dual-model selection (Sonar Pro vs Sonar Reasoning Pro)
- Complexity scoring algorithm with explicit thresholds (≥3 points → reasoning)
- Keyword triggers documented
- Query structure template

**Improvements Needed:**
- Add decision tree diagram for model selection
- Quantify expected response quality by model
- Include fallback strategies when API fails

---

### scientific-critical-thinking (7.5/10)

**Strengths:**
- Comprehensive evaluation framework
- Bias taxonomy with examples
- Evidence quality hierarchy
- Logical fallacy catalog

**Improvements Needed:**
- Add decision tree for selecting appropriate analysis framework
- Numbered workflow steps (currently prose-heavy)
- Scoring rubric for evidence quality levels

---

### hypothesis-generation (7.0/10)

**Strengths:**
- 8-step workflow with numbered stages
- LaTeX template with overflow prevention rules
- Box type specifications (hypothesisbox1-5, predictionbox, etc.)
- Page break strategy documented
- Citation targets (10-15 main text, 50+ total)

**Improvements Needed:**
- Add decision tree for hypothesis type selection
- Quantify "viable" hypothesis criteria
- Include workflow transition checklists

---

### venue-templates (7.0/10)

**Strengths:**
- Journal-specific formatting guides
- Style comparison tables
- Reviewer expectation documentation
- Word limit specifications

**Improvements Needed:**
- Add decision tree for journal selection based on paper characteristics
- Consolidate scattered venue guides into unified decision matrix
- Include rejection rate and review timeline data

---

## Tier 3: Moderate Determinism (5.5-6.9)

### literature-review (6.5/10)

**Strengths:**
- 7-phase systematic methodology
- PICO framework for structuring questions
- Multi-database integration documented
- PRISMA flow diagram requirement

**Improvements Needed:**
- Add decision tree for database selection by field
- Quantify screening criteria (e.g., "title screen ≤30 seconds per paper")
- Include quality assessment tool selection matrix
- Define "sufficient" search saturation criteria

---

### citation-management (6.5/10)

**Strengths:**
- Style-specific formatting rules
- DOI verification emphasis
- Reference manager integration notes

**Improvements Needed:**
- Add decision tree for citation style selection by venue
- Workflow for bulk citation cleanup
- Quantify "complete" citation (required vs optional fields)
- Include common error patterns with fixes

---

### scientific-schematics (6.5/10)

**Strengths:**
- AI generation workflow documented
- Accessibility requirements (colorblind-safe, high contrast)
- Output format specifications

**Improvements Needed:**
- Add decision tree for schematic type selection
- Quantify iteration cycles ("review and refine" → "≤3 iterations")
- Include quality checklist with pass/fail criteria
- Define style consistency rules

---

### visual-design (6.5/10)

**Strengths:**
- Four-question design thinking framework
- Accessibility requirements with WCAG ratios
- Chart selection anti-patterns
- Publication requirements table

**Improvements Needed:**
- Convert prose sections to decision matrices
- Add workflow steps with transition criteria
- Quantify "memorable" and "distinctive" criteria
- Include brand compliance checklist

---

### scientific-slides (6.0/10)

**Strengths:**
- Comprehensive content (~1150 lines)
- Multiple presentation type guides
- Visual hierarchy principles

**Improvements Needed:**
- **Critical**: Split into focused sub-documents (currently overwhelming)
- Add decision tree for slide type selection
- Create numbered workflow (currently prose paragraphs)
- Define explicit success criteria per slide type
- Quantify "appropriate" animation and transition use

---

### latex-posters (6.0/10)

**Strengths:**
- Template specifications
- Size and layout options
- Font and color guidelines

**Improvements Needed:**
- Add decision tree for poster format selection
- Create content density guidelines (text per section)
- Include quality checklist
- Define accessibility requirements

---

### pptx-posters (6.0/10)

**Strengths:**
- Tool-specific guidance (PowerPoint focus)
- Size specifications
- Export instructions

**Improvements Needed:**
- Add decision tree for PPTX vs LaTeX selection
- Create content flow template
- Include accessibility checklist
- Define print-ready criteria

---

## Tier 4: Needs Improvement (<5.5)

### generate-image (5.5/10)

**Strengths:**
- Photorealistic illustration focus
- Medical/anatomical emphasis

**Improvements Needed:**
- Add decision tree for image generation approach
- Define prompt engineering patterns
- Create quality validation checklist
- Include iteration workflow with exit criteria

---

### scholar-evaluation (5.0/10)

**Strengths:**
- Metric definitions (h-index, citations)

**Improvements Needed:**
- Add evaluation rubric with scoring
- Create decision tree for metric interpretation
- Define field-specific normalization rules
- Include workflow for comprehensive evaluation

---

### markitdown (5.0/10)

**Strengths:**
- Conversion focus
- Format compatibility notes

**Improvements Needed:**
- Add decision tree for conversion path selection
- Create quality validation checklist
- Define edge case handling
- Include workflow with pre/post-processing steps

---

### paper-2-web (5.0/10)

**Strengths:**
- Web output focus
- Interactive element mentions

**Improvements Needed:**
- Add decision tree for web format selection
- Create accessibility requirements
- Define SEO considerations
- Include validation checklist

---

## Cross-Cutting Issues

### 1. Inconsistent Structure Patterns

| Pattern | Best Skills | Missing In |
|---------|-------------|------------|
| XML semantic tags | statistical-analysis, reproducible-research | scientific-slides, visual-design |
| Decision matrices | statistical-analysis, plotting-libraries | literature-review, citation-management |
| Numbered workflows | hypothesis-generation, code-documentation | generate-image, markitdown |
| Quantitative thresholds | reproducible-research | most older skills |

### 2. Missing Workflow Transition Criteria

Most skills define phases but lack explicit "exit criteria" for moving between stages:

**Example Gap:**
```
❌ Current: "After screening, proceed to data extraction"
✅ Improved: "Proceed to data extraction when:
   □ All duplicates removed (0 remaining)
   □ Title/abstract screening complete (100% reviewed)
   □ Exclusion reasons documented for all rejected papers
   □ Inter-rater agreement ≥0.8 (if dual screening)"
```

### 3. Qualitative Terms Without Quantification

| Qualitative Term | Appears In | Suggested Quantification |
|-----------------|------------|-------------------------|
| "comprehensive" | literature-review, scientific-writing | ≥50 papers, ≥3 databases |
| "high quality" | peer-review, visual-design | score ≥8/10 on rubric |
| "appropriate" | statistical-analysis, venue-templates | meets 5/5 selection criteria |
| "sufficient" | hypothesis-generation | ≥3 competing hypotheses |

### 4. Incomplete Skill Orchestration

When multiple skills apply, routing logic is implicit:

**Example Gap:**
- User wants to create a figure → visual-design? plotting-libraries? scientific-schematics?
- No explicit decision tree for skill selection

---

## Implementation Plan

### Phase 1: Structure Standardization (High Priority)

**Goal:** Establish consistent XML-based structure across all skills

**Duration:** 1-2 sessions per skill

#### Step 1.1: Create Structure Template

Create a standard SKILL.md template with required sections:

```markdown
---
name: skill-name
version: x.y.z
description: "..."
allowed-tools: [...]
---

# Skill Name

<overview>
Brief description of what this skill does and when to use it.
</overview>

<when_to_use>
- Trigger condition 1
- Trigger condition 2
</when_to_use>

<decision_framework>
## Decision Matrix
| Condition | → Action |
|-----------|----------|
| ... | ... |
</decision_framework>

<workflow>
## Workflow
1. **Stage 1: Name**
   - Step details
   - Exit criteria: [checklist]

2. **Stage 2: Name**
   ...
</workflow>

<success_criteria>
## Success Criteria
□ Criterion 1 (quantified)
□ Criterion 2 (quantified)
</success_criteria>

<scope>
## Scope
**In Scope:** ...
**Out of Scope:** ...
</scope>

<anti_patterns>
## Common Pitfalls
1. **Pitfall Name**: Description → Solution
</anti_patterns>

<cross_references>
## Related Skills
- skill-name: relationship description
</cross_references>
```

#### Step 1.2: Apply Template to Tier 4 Skills

Priority order:
1. `scholar-evaluation` - Add evaluation rubric and scoring
2. `markitdown` - Add conversion decision tree
3. `paper-2-web` - Add web format selection criteria
4. `generate-image` - Add prompt engineering patterns

**Per-skill tasks:**
- [ ] Add XML semantic tags
- [ ] Create decision matrix
- [ ] Add numbered workflow
- [ ] Define success criteria
- [ ] Add scope boundaries

#### Step 1.3: Apply Template to Tier 3 Skills

Priority order:
1. `scientific-slides` - Split into sub-documents, add decision tree
2. `latex-posters` - Add format selection criteria
3. `pptx-posters` - Add PPTX vs LaTeX decision
4. `literature-review` - Add database selection matrix
5. `citation-management` - Add style selection decision tree
6. `scientific-schematics` - Add type selection criteria
7. `visual-design` - Convert prose to matrices

---

### Phase 2: Quantification (Medium Priority)

**Goal:** Replace qualitative terms with measurable criteria

**Duration:** 30-60 minutes per skill

#### Step 2.1: Audit Qualitative Terms

For each skill, identify and replace:

| Original | Quantified Replacement |
|----------|----------------------|
| "comprehensive review" | "≥50 papers from ≥3 databases" |
| "high quality figure" | "meets 8/10 quality checklist items" |
| "appropriate statistical test" | "passes 5/5 decision tree criteria" |
| "sufficient sample size" | "≥80% power at α=0.05" |
| "thorough documentation" | "≥90% public API coverage" |

#### Step 2.2: Add Threshold Tables

Create reference tables for common thresholds:

```markdown
## Quality Thresholds

| Metric | Minimum | Target | Excellent |
|--------|---------|--------|-----------|
| Citation count | 20 | 50 | 100+ |
| Figure resolution | 150 DPI | 300 DPI | 600 DPI |
| Code coverage | 60% | 80% | 95% |
| Inter-rater agreement | 0.6 | 0.8 | 0.9 |
```

---

### Phase 3: Workflow Transitions (Medium Priority)

**Goal:** Add explicit exit criteria between workflow stages

**Duration:** 20-30 minutes per skill

#### Step 3.1: Identify Stage Boundaries

For each skill with numbered workflows, add transition checklists:

```markdown
### Stage 2 → Stage 3 Transition

**Exit Criteria (all must be true):**
□ All Stage 2 outputs generated
□ Quality check passed (score ≥7/10)
□ No blocking issues identified
□ Required approvals obtained

**Entry Requirements for Stage 3:**
□ Stage 2 artifacts available
□ Dependencies resolved
□ Resources allocated
```

#### Step 3.2: Add Progress Indicators

Define measurable progress within each stage:

```markdown
### Stage 2: Data Collection

**Progress Milestones:**
- 25%: Query strategy finalized
- 50%: Primary database search complete
- 75%: Secondary databases searched
- 100%: Deduplication complete, screening list ready
```

---

### Phase 4: Skill Orchestration (High Priority)

**Goal:** Create explicit routing logic for multi-skill scenarios

**Duration:** 1 session

#### Step 4.1: Create Skill Router Document

Add `SKILL_ROUTER.md` at skills root:

```markdown
# Skill Router

## Decision Tree: Figure Creation

```
User wants to create a figure
│
├─ Is it a data visualization (plots, charts)?
│  │
│  ├─ Yes → plotting-libraries
│  │         ├─ Statistical comparison → seaborn/ggplot2
│  │         ├─ Time series → matplotlib/ggplot2
│  │         └─ Bioinformatics → Bioconductor
│  │
│  └─ No → Is it a schematic/diagram?
│          │
│          ├─ Yes → scientific-schematics
│          │
│          └─ No → Is it photorealistic?
│                   │
│                   ├─ Yes → generate-image
│                   │
│                   └─ No → visual-design (for guidance)
```

## Decision Tree: Document Creation

```
User wants to create a document
│
├─ Is it a research paper?
│  │
│  ├─ Yes → scientific-writing
│  │        └─ Then → venue-templates (for journal-specific formatting)
│  │
│  └─ No → Is it a presentation?
│          │
│          ├─ Yes → scientific-slides
│          │
│          └─ No → Is it a poster?
│                   │
│                   ├─ Yes → latex-posters OR pptx-posters
│                   │        (decision: LaTeX if complex equations,
│                   │         PPTX if quick turnaround)
│                   │
│                   └─ No → [route by document type]
```
```

#### Step 4.2: Add Router References

In each skill's cross-references section, add:

```markdown
**Skill Selection:**
See `SKILL_ROUTER.md` for decision trees when multiple skills may apply.
```

---

### Phase 5: Validation (Low Priority)

**Goal:** Verify improvements through testing

**Duration:** Ongoing

#### Step 5.1: Create Test Scenarios

For each skill, define 3-5 test scenarios:

```markdown
## Test Scenarios

### Scenario 1: Basic Usage
**Input:** [description]
**Expected Behavior:** [deterministic output]
**Success Criteria:** [measurable]

### Scenario 2: Edge Case
**Input:** [edge case description]
**Expected Behavior:** [how skill handles it]
**Success Criteria:** [measurable]
```

#### Step 5.2: Conduct Periodic Audits

Schedule quarterly reviews:
- Re-score each skill on 8 dimensions
- Track improvement trends
- Identify new patterns to propagate

---

## Exemplary Patterns to Propagate

### Pattern 1: Decision Matrix (from statistical-analysis)

```markdown
## Two Groups Comparison

| Data Type | Distribution | Paired? | → Test |
|-----------|--------------|---------|--------|
| Continuous | Normal | No | Independent t-test |
| Continuous | Normal | Yes | Paired t-test |
| Continuous | Non-normal | No | Mann-Whitney U |
| Continuous | Non-normal | Yes | Wilcoxon signed-rank |
| Categorical | N/A | No | Chi-square |
| Categorical | N/A | Yes | McNemar's |
```

### Pattern 2: Numbered Pitfalls (from reproducible-research)

```markdown
## Common Pitfalls

1. **Unpinned Dependencies**
   - Anti-pattern: `pandas` (no version)
   - Solution: `pandas==2.0.3` (exact pin)

2. **Missing Random Seeds**
   - Anti-pattern: Results vary between runs
   - Solution: `np.random.seed(42)` at script start
```

### Pattern 3: Exit Criteria Checklist (from code-documentation)

```markdown
### Documentation Complete Checklist

□ README.md has all 7 required sections
□ ≥80% of public functions have docstrings
□ All parameters have type hints
□ Examples compile and run successfully
□ API reference generated without errors
```

### Pattern 4: Scope Boundaries (from statistical-analysis)

```markdown
## Scope

**In Scope:**
- Frequentist hypothesis testing
- Common parametric and non-parametric tests
- Multiple testing correction
- Effect size calculation

**Out of Scope (use specialized resources):**
- Bayesian inference
- Time series analysis
- Causal inference methods
- Deep learning
```

### Pattern 5: Template Strings (from hypothesis-generation)

```markdown
## Reporting Template

```latex
\begin{hypothesisbox1}[Hypothesis 1: {{TITLE}}]
\textbf{Mechanistic Explanation:}
{{1-2 paragraphs, 6-10 sentences max}}

\textbf{Key Supporting Evidence:}
\begin{itemize}
\item {{Evidence 1 with citation}}
\item {{Evidence 2 with citation}}
\end{itemize}

\textbf{Core Assumptions:}
\begin{itemize}
\item {{Assumption 1}}
\end{itemize}
\end{hypothesisbox1}
```
```

---

## Summary

### Current State
- **Overall Score:** 7.5/10
- **Top Performers:** statistical-analysis (9.5), reproducible-research (9.0), code-documentation (8.5)
- **Needs Work:** scholar-evaluation (5.0), markitdown (5.0), paper-2-web (5.0)

### Target State
- **Target Score:** 8.5/10
- **Timeline:** 4-6 focused sessions
- **Key Improvements:** Structure standardization, quantification, workflow transitions, skill orchestration

### Priority Actions

| Priority | Action | Impact | Effort |
|----------|--------|--------|--------|
| 🔴 High | Create SKILL_ROUTER.md | Solves multi-skill confusion | 1 session |
| 🔴 High | Apply template to Tier 4 skills | Brings 4 skills to baseline | 2 sessions |
| 🟡 Medium | Quantify qualitative terms | Improves all skills | 2-3 hours |
| 🟡 Medium | Add workflow transition criteria | Improves 12 skills | 2 hours |
| 🟢 Low | Create test scenarios | Validates improvements | Ongoing |

---

*This audit provides a roadmap for improving skill determinism. Each phase can be executed independently, and improvements compound across the library.*
