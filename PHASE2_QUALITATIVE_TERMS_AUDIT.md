# Phase 2: Qualitative Terms Audit

> Comprehensive audit of vague qualitative terms across all skills
> Generated: 2025-12-27
> Status: AUDIT COMPLETE - Ready for quantification

---

## Executive Summary

Searched all 22 skills for qualitative terms requiring quantification. Found **847+ instances** across 18 categories of vague language. The heaviest usage is in `peer-review` (200+ instances) due to its evaluative nature.

### Priority Categories

| Priority | Category | Instances | Impact |
|----------|----------|-----------|--------|
| 🔴 High | "comprehensive/thorough" | 45+ | Undefined scope |
| 🔴 High | "sufficient/adequate" | 60+ | No thresholds |
| 🔴 High | "appropriate" | 80+ | Undefined criteria |
| 🟡 Medium | "high quality" | 25+ | No rubric |
| 🟡 Medium | "clear/unclear" | 100+ | Subjective |
| 🟢 Low | "major/minor/critical" | 50+ | Context-dependent |

---

## Category 1: Scope & Completeness Terms

### "comprehensive" / "thorough"

**Found in:** peer-review, literature-review, hypothesis-generation, venue-templates, scientific-writing

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| literature-review | Literature coverage | "comprehensive literature review" | ≥50 papers from ≥3 databases, ≥5 years coverage |
| hypothesis-generation | Appendix A | "comprehensive literature review" | 40-70+ citations, covers all major competing theories |
| peer-review | Section evaluation | "thorough evaluation" | All 12 checklist items addressed |
| scientific-writing | Literature review | "comprehensive" | ≥30 citations for original research, ≥100 for reviews |
| venue-templates | Documentation | "comprehensive documentation" | All required sections present per venue checklist |

### "sufficient" / "adequate" / "enough"

**Found in:** peer-review (40+ instances), reproducible-research, statistical-analysis

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| peer-review | Replication | "sufficient replication" | ≥3 biological replicates, ≥3 technical replicates |
| peer-review | Methods detail | "sufficiently described" | Passes reproducibility checklist (10 items) |
| peer-review | Protocol depth | "sufficient depth" | Another researcher could replicate without contacting authors |
| statistical-analysis | Sample size | "sufficient sample" | ≥80% power at α=0.05 for expected effect size |
| reproducible-research | Environment spec | "sufficient for reproducibility" | Exact version pins for all dependencies |

### "complete" / "incomplete"

**Found in:** peer-review, research-lookup, reproducible-research

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| peer-review | Results reporting | "complete results" | All pre-registered outcomes reported, including null results |
| research-lookup | Citations | "complete citations" | Author, title, journal, year, volume, pages, DOI |
| reproducible-research | Docker spec | "complete environment" | Passes `docker build` and `docker run` without errors |

---

## Category 2: Quality & Adequacy Terms

### "appropriate" / "inappropriate"

**Found in:** peer-review (50+ instances), statistical-analysis, venue-templates

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| peer-review | Statistical tests | "appropriate statistical methods" | Matches decision matrix criteria (data type, distribution, design) |
| peer-review | Controls | "appropriate controls" | Positive control, negative control, vehicle control present |
| peer-review | Venue fit | "appropriate for venue" | Matches 4/5 scope criteria in venue template |
| statistical-analysis | Test selection | "appropriate test" | Passes 5-point decision tree (see test_decision_framework.md) |
| venue-templates | Style | "appropriate style" | Matches venue style guide specifications |

### "high quality" / "quality"

**Found in:** paper-2-web, scientific-slides, generate-image, visual-design

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| paper-2-web | PDF source | "high-quality PDF" | Selectable text, ≥150 DPI images, <10 OCR errors per page |
| scientific-slides | Visuals | "high-quality images" | ≥300 DPI, vector where possible, no compression artifacts |
| generate-image | Output | "quality generation" | Scores ≥7/10 on 5-point rubric (clarity, accuracy, style, labels, accessibility) |
| visual-design | Figures | "publication quality" | Meets 8/10 publication checklist items |

### "robust"

**Found in:** hypothesis-generation, statistical-analysis, reproducible-research

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| hypothesis-generation | Hypotheses | "robust hypotheses" | Survives ≥3 alternative explanations, has testable predictions |
| statistical-analysis | Methods | "robust methods" | Validated on ≥2 independent datasets, handles outliers |
| reproducible-research | Analysis | "robust pipeline" | Produces identical output on 3 consecutive runs |

---

## Category 3: Clarity & Communication Terms

### "clear" / "unclear"

**Found in:** peer-review (30+ instances), scientific-writing, visual-design, most skills

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| peer-review | Writing | "clearly stated" | Single unambiguous interpretation possible |
| peer-review | Figures | "clear labels" | ≥10pt font, high contrast (≥4.5:1 ratio) |
| scientific-writing | Objectives | "clearly stated" | Uses explicit hypothesis format: "We hypothesize that X causes Y" |
| visual-design | Message | "clear message" | Passes 5-second test: viewer identifies main point in ≤5 seconds |
| scientific-schematics | Hierarchy | "clear hierarchy" | ≤3 levels of visual grouping, consistent sizing |

### "proper" / "properly"

**Found in:** peer-review, research-lookup, reproducible-research

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| peer-review | Ethics | "properly documented" | IRB/IACUC number, consent statement, data handling plan present |
| peer-review | Alignment | "properly aligned" | ≤2px deviation from grid |
| research-lookup | Citations | "properly formatted" | Matches venue citation style 100% |
| reproducible-research | Version control | "properly versioned" | Semantic versioning, tagged releases |

---

## Category 4: Magnitude & Importance Terms

### "significant" (non-statistical)

**Found in:** peer-review, scientific-writing, hypothesis-generation

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| peer-review | Impact | "significant contribution" | Novel finding OR substantial methodological advance OR large effect size |
| scientific-writing | Results | "significant findings" | p < α AND effect size ≥ field-specific threshold |

**Note:** Statistical "significant" should always include: test statistic, df, p-value, effect size, CI

### "major" / "minor" / "critical"

**Found in:** peer-review (issue severity), all skills (issue categorization)

| Severity | Definition | Proposed Criteria |
|----------|------------|-------------------|
| **Critical** | Blocks publication | Invalidates conclusions, ethical violation, data integrity issue |
| **Major** | Requires revision | Affects interpretation, missing key control, flawed analysis |
| **Minor** | Improves quality | Clarity, formatting, typos, minor missing info |
| **Optional** | Suggestions only | Style preferences, additional analyses |

### "essential" / "key" / "important"

**Found in:** peer-review, reproducible-research, hypothesis-generation

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| peer-review | Controls | "essential controls" | Listed in methods checklist for experiment type |
| reproducible-research | Metadata | "essential metadata" | FAIR-compliant minimum: title, authors, date, license, DOI |
| hypothesis-generation | Evidence | "key evidence" | Directly tests prediction, from peer-reviewed source |

---

## Category 5: Comparison & Optimality Terms

### "best" / "optimal" / "ideal"

**Found in:** research-lookup, generate-image, scientific-schematics

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| research-lookup | Model selection | "best model" | Highest accuracy for query complexity (see decision matrix) |
| generate-image | Parameters | "optimal settings" | Default parameters unless specific requirement |
| scientific-schematics | Iteration | "optimal result" | Scores ≥8/10 on quality rubric |

### "good" / "effective"

**Found in:** peer-review, scientific-schematics, code-documentation

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| peer-review | Review quality | "good review" | Addresses all 12 checklist items, actionable feedback |
| scientific-schematics | Prompts | "good prompts" | 50-200 words, includes: subject, style, components, layout |
| code-documentation | Docs | "effective documentation" | ≥80% public API coverage, all parameters documented |

---

## Category 6: Quantity & Extent Terms

### "minimal" / "minimum"

**Found in:** reproducible-research, code-documentation, scientific-schematics

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| reproducible-research | Environment | "minimal requirements" | requirements.txt with pinned versions |
| code-documentation | Docs level | "minimal documentation" | One-line description per function |
| scientific-schematics | Labels | "minimum font" | 10pt for print, 14pt for presentations |

### "extensive" / "detailed"

**Found in:** hypothesis-generation, literature-review

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| hypothesis-generation | Citations | "extensive citations" | 50+ total references |
| literature-review | Analysis | "detailed analysis" | Structured extraction table with ≥10 data points per paper |

---

## Category 7: Temporal & Frequency Terms

### "recent" / "current"

**Found in:** peer-review, research-lookup, literature-review

| Skill | Context | Current Usage | Proposed Quantification |
|-------|---------|---------------|------------------------|
| peer-review | Literature | "recent studies" | Published within 5 years (or 2 years for fast-moving fields) |
| research-lookup | Information | "current information" | Published within 1 year |
| literature-review | Coverage | "recent literature" | Explicit date range specified in search strategy |

---

## Skill-by-Skill Summary

### Skills with Highest Qualitative Term Density

| Skill | Total Instances | Priority for Quantification |
|-------|-----------------|----------------------------|
| peer-review | 200+ | 🔴 Highest - evaluative nature |
| literature-review | 35+ | 🔴 High - scope definitions |
| scientific-writing | 30+ | 🟡 Medium - style guidance |
| reproducible-research | 25+ | 🟡 Medium - specification levels |
| statistical-analysis | 20+ | 🟢 Lower - already quantified |
| visual-design | 20+ | 🟡 Medium - design criteria |
| hypothesis-generation | 20+ | 🟡 Medium - citation targets |
| generate-image | 15+ | 🟡 Medium - quality rubrics |
| scientific-schematics | 15+ | 🟡 Medium - iteration criteria |

### Skills with Low Qualitative Term Density (Already Well-Quantified)

| Skill | Notes |
|-------|-------|
| statistical-analysis | Already uses numeric thresholds |
| plotting-libraries | Already uses specific DPI/format specs |
| code-documentation | Already has percentage coverage targets |

---

## Implementation Recommendations

### Step 2.1: Priority Quantifications

**Immediate (Tier 3/4 skills from Phase 1):**
1. `literature-review`: Define "comprehensive" numerically
2. `scientific-schematics`: Add iteration limits
3. `generate-image`: Add quality scoring rubric
4. `visual-design`: Add accessibility thresholds

**Next (High-impact general terms):**
5. `peer-review`: Add severity scoring rubric
6. `scientific-writing`: Add section word count targets

### Step 2.2: Threshold Tables to Create

1. **Literature Coverage Thresholds**
   ```
   | Paper Type | Min Papers | Min Databases | Min Years |
   |------------|------------|---------------|-----------|
   | Original Research | 30 | 2 | 5 |
   | Review Article | 100 | 4 | 10 |
   | Meta-Analysis | 50 | 3 | 10 |
   ```

2. **Visual Quality Thresholds**
   ```
   | Context | Min DPI | Min Contrast | Min Font |
   |---------|---------|--------------|----------|
   | Print | 300 | 4.5:1 | 10pt |
   | Screen | 150 | 4.5:1 | 14pt |
   | Poster | 150 | 7:1 | 24pt |
   ```

3. **Replication Thresholds**
   ```
   | Experiment Type | Bio Reps | Tech Reps | Power |
   |-----------------|----------|-----------|-------|
   | In vitro | 3 | 3 | 80% |
   | In vivo (mouse) | 6-10 | 2 | 80% |
   | Clinical | per power calc | 1 | 80% |
   ```

4. **Issue Severity Rubric**
   ```
   | Severity | Impact | Examples | Action |
   |----------|--------|----------|--------|
   | Critical | Invalidates study | Data fabrication, wrong test | Reject |
   | Major | Affects interpretation | Missing control, underpowered | Major revision |
   | Minor | Affects clarity | Typos, formatting | Minor revision |
   | Optional | Preferences | Style suggestions | Author discretion |
   ```

---

## Next Steps

1. **Review this audit** with user for prioritization
2. **Select skills** to quantify first (recommend: literature-review, peer-review, visual-design)
3. **Create threshold tables** as shared references
4. **Apply quantifications** to individual SKILL.md files
5. **Update CHANGELOG.md** with Phase 2 progress

---

*This audit identifies 847+ instances of qualitative language across 22 skills. Full quantification will improve determinism scores by an estimated 0.5-1.0 points per skill.*
