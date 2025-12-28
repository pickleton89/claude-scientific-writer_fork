---
name: literature-review
version: 2.1.0
description: "Conduct systematic literature reviews using multiple academic databases with quantified coverage thresholds. Provides structured methodology for search, screening, synthesis, and citation verification with professional output in markdown and PDF formats."
allowed-tools: [Read, Write, Edit, Bash, Glob, Grep, WebFetch, WebSearch]
shared-thresholds: "../QUANTIFICATION_THRESHOLDS.md"
---

# Literature Review

> **Quantified Thresholds:** This skill references shared thresholds from [`QUANTIFICATION_THRESHOLDS.md`](../QUANTIFICATION_THRESHOLDS.md) §1 (Literature Coverage) and §6 (Time-Based Thresholds).

<overview>
Conduct systematic literature reviews with quantified coverage requirements: ≥3 databases, 30-100+ papers depending on review type, and documented search strategies achieving <5% new papers in final iteration. This skill provides structured workflows for multi-database searching, screening, thematic synthesis, and citation verification with professional output generation in markdown and PDF formats.
</overview>

<when_to_use>
## Trigger Conditions

Use this skill when:
- Conducting a systematic literature review for research or publication
- Synthesizing current knowledge on a topic across multiple sources
- Performing meta-analysis, scoping reviews, or rapid reviews
- Writing the literature review section of a paper or thesis
- Investigating state of the art in a research domain
- Identifying research gaps and future directions

Do NOT use this skill when:
- Looking up a single paper or citation → use `research-lookup`
- Searching a single database for quick reference → use `research-lookup`
- Formatting existing citations → use `citation-management`
- Creating figures for the review → use `scientific-schematics` or `plotting-libraries`
</when_to_use>

<decision_framework>
## Database Selection Matrix

```
What is the primary research domain?
│
├─ Biomedical / Life Sciences
│  │
│  ├─ Clinical / Medical → PubMed + Cochrane + CINAHL + Embase
│  │
│  ├─ Basic Biology → PubMed + bioRxiv + Web of Science
│  │
│  └─ Genomics/Bioinformatics → PubMed + bioRxiv + GEO + UniProt
│
├─ Physical Sciences / Engineering
│  │
│  ├─ Physics / Math → arXiv + Web of Science + ADS
│  │
│  ├─ Chemistry → PubMed + SciFinder + Web of Science
│  │
│  └─ Engineering → IEEE Xplore + Scopus + Web of Science
│
├─ Computer Science / AI
│  │
│  ├─ ML/AI → arXiv + Semantic Scholar + ACM DL + IEEE
│  │
│  └─ General CS → ACM DL + IEEE Xplore + arXiv + DBLP
│
├─ Social Sciences
│  │
│  └─ Psychology / Sociology → PsycINFO + Scopus + SSRN
│
└─ Multidisciplinary
   │
   └─ Cross-field → Semantic Scholar + Web of Science + Scopus + Google Scholar
```

## Minimum Database Requirements

| Review Type | Minimum Databases | Recommended |
|-------------|------------------|-------------|
| Systematic review | 3 | 4-5 |
| Scoping review | 2 | 3-4 |
| Narrative review | 2 | 3 |
| Rapid review | 2 | 2-3 |

## Review Type Selection

| Factor | Systematic | Scoping | Narrative | Rapid |
|--------|-----------|---------|-----------|-------|
| Timeline | 6-18 months | 3-6 months | 1-3 months | 2-5 weeks |
| Question type | Focused, specific | Broad, exploratory | Broad overview | Urgent, focused |
| Methodology | Rigorous protocol | Flexible framework | Expert synthesis | Streamlined |
| Quality assessment | Required | Optional | Informal | Limited |
| PRISMA required | Yes | PRISMA-ScR | No | Modified |
| Output papers | ≥30 typical | Variable | Variable | 10-30 typical |

## Quality Assessment Tool Selection

| Study Design | Assessment Tool | Domains Evaluated |
|-------------|----------------|-------------------|
| Randomized trials | Cochrane RoB 2.0 | 5 domains, 7 items |
| Non-randomized interventions | ROBINS-I | 7 domains |
| Observational cohort | Newcastle-Ottawa Scale | 3 domains, 8 items |
| Cross-sectional | JBI Checklist | 8 items |
| Diagnostic accuracy | QUADAS-2 | 4 domains |
| Qualitative studies | CASP Checklist | 10 items |
| Systematic reviews | AMSTAR 2 | 16 items, 7 critical |
| Animal studies | SYRCLE RoB | 10 items |

## Literature Coverage Thresholds

> Reference: [`QUANTIFICATION_THRESHOLDS.md`](../QUANTIFICATION_THRESHOLDS.md) §1

| Review Type | Min Papers | Min Databases | Time Range | Citation Density |
|-------------|------------|---------------|------------|------------------|
| Original Research | 30 | 2 | 5-10 years | 2-4 per paragraph |
| Review Article | 100 | 4 | 10+ years | 4-6 per paragraph |
| Systematic Review | 50+ | 3 | Per PRISMA | All relevant |
| Meta-Analysis | 10+ studies | 3 | 10+ years | All included |
| Methods Paper | 20 | 2 | 5 years | Key comparisons |
| Commentary/Letter | 10-15 | 1 | 2-3 years | Selective |

**"Comprehensive" Defined:**

A literature search is considered **comprehensive** when ALL criteria are met:
- [ ] ≥3 databases searched with documented strategy
- [ ] Forward/backward citation search performed
- [ ] Gray literature checked (preprints, theses, conference proceedings)
- [ ] ≤5% new relevant papers found in final search iteration
- [ ] Date range explicitly specified and justified

## Literature Recency Thresholds

> Reference: [`QUANTIFICATION_THRESHOLDS.md`](../QUANTIFICATION_THRESHOLDS.md) §6

| Field Type | "Recent" Means | "Current" Means |
|------------|---------------|-----------------|
| Fast-moving (ML, genomics) | ≤2 years | ≤6 months |
| Standard biomedical | ≤5 years | ≤1 year |
| Clinical guidelines | ≤3 years | ≤1 year |
| Historical context | ≤10 years | ≤5 years |
| Foundational work | Any age if seminal | N/A |

</decision_framework>

<workflow>
## Workflow

### Stage 1: Planning and Protocol

**Objective:** Define scope, research questions, and methodology

**Steps:**
1. Define research question using PICO/PEO/SPIDER framework:
   - **PICO**: Population, Intervention, Comparison, Outcome (clinical)
   - **PEO**: Population, Exposure, Outcome (observational)
   - **SPIDER**: Sample, Phenomenon, Design, Evaluation, Research type (qualitative)

2. Register protocol (for systematic reviews):
   - PROSPERO (health): https://www.crd.york.ac.uk/prospero/
   - OSF Registries (general): https://osf.io/registries

3. Define inclusion/exclusion criteria:

| Criterion | Include | Exclude |
|-----------|---------|---------|
| Date range | {{start}}-{{end}} | Before {{start}} |
| Language | English | Non-English |
| Publication type | Peer-reviewed, preprints | Conference abstracts only |
| Study design | {{designs}} | {{excluded designs}} |
| Population | {{target}} | {{excluded populations}} |

4. Select databases (minimum 3 for systematic reviews)
5. Define search strategy with Boolean operators

**Exit Criteria:**
- [ ] Research question documented (PICO/PEO/SPIDER format)
- [ ] Protocol registered (if systematic review)
- [ ] Inclusion/exclusion criteria table completed
- [ ] ≥3 databases selected (see Database Selection Matrix)
- [ ] Search strategy drafted

---

### Stage 2: Systematic Search

**Objective:** Execute documented database searches meeting coverage thresholds (≥3 databases, <5% new papers in final iteration)

**Steps:**
1. Construct search strings for each database:
   ```
   Concept 1 terms: (term1 OR synonym1 OR synonym2)
   AND
   Concept 2 terms: (term2 OR synonym3 OR synonym4)
   AND
   Concept 3 terms: (term3 OR synonym5)
   ```

2. Execute searches and document:

| Database | Date Searched | Search String | Results |
|----------|--------------|---------------|---------|
| PubMed | YYYY-MM-DD | [full query] | N papers |
| Scopus | YYYY-MM-DD | [full query] | N papers |
| ... | ... | ... | ... |

3. Export results in structured format (RIS, BibTeX, JSON)
4. Search citation networks:
   - Forward citations: Papers citing key papers
   - Backward citations: References from included papers
5. Search grey literature if appropriate (theses, reports, preprints)

**Search Saturation Criteria:**
- [ ] No new relevant papers in last 50 screened titles
- [ ] Citation chaining yields <5% new papers
- [ ] Expert consultation confirms no major gaps

**Exit Criteria:**
- [ ] All selected databases searched
- [ ] Search strategy documented verbatim for each database
- [ ] Search dates recorded
- [ ] All results exported and combined
- [ ] Total unique results counted after initial deduplication

---

### Stage 3: Screening and Selection

**Objective:** Systematically screen papers against criteria

**Screening Time Budgets:**

| Stage | Time per Paper | Action |
|-------|---------------|--------|
| Title screening | ≤15 seconds | Include/Exclude/Maybe |
| Abstract screening | ≤2 minutes | Include/Exclude with reason |
| Full-text screening | 10-20 minutes | Final inclusion decision |

**Steps:**
1. Remove duplicates:
   - Primary: Match by DOI
   - Secondary: Match by title (normalized, ≥90% similarity)
   - Document: N duplicates removed

2. Title screening:
   - Apply inclusion criteria to titles only
   - When uncertain, include for abstract review
   - Document exclusion count

3. Abstract screening:
   - Read abstracts of remaining papers
   - Apply all inclusion/exclusion criteria
   - Record exclusion reason for each paper

4. Full-text screening:
   - Obtain full texts of remaining papers
   - Apply all criteria with detailed review
   - Record specific exclusion reasons

5. Create PRISMA flow diagram:
```
Records identified (n = X)
├── Database 1 (n = a)
├── Database 2 (n = b)
└── Other sources (n = c)
│
▼
Records after duplicates removed (n = Y)
│
▼
Records screened (n = Y)
├── Excluded at title (n = d)
│
▼
Abstracts assessed (n = Z)
├── Excluded at abstract (n = e)
│   ├── Reason 1 (n = e1)
│   ├── Reason 2 (n = e2)
│   └── Reason 3 (n = e3)
│
▼
Full-text assessed (n = A)
├── Excluded (n = f)
│   ├── Reason 1 (n = f1)
│   ├── Reason 2 (n = f2)
│   └── Reason 3 (n = f3)
│
▼
Studies included (n = B)
├── In quantitative synthesis (n = B1)
└── In qualitative synthesis (n = B2)
```

**Exit Criteria:**
- [ ] Deduplication complete (documented count)
- [ ] All titles screened (≤15 sec each average)
- [ ] All abstracts screened (exclusion reasons recorded)
- [ ] All full texts assessed
- [ ] PRISMA flow diagram complete with numbers
- [ ] Inter-rater agreement ≥0.8 (if dual screening)

---

### Stage 4: Data Extraction

**Objective:** Extract standardized data from included studies

**Extraction Template:**

| Field | Description |
|-------|-------------|
| Study ID | First author, year |
| Country | Study location |
| Design | RCT, cohort, case-control, etc. |
| Sample size | N participants/samples |
| Population | Demographics, inclusion criteria |
| Intervention/Exposure | Details of what was studied |
| Comparator | Control condition |
| Outcomes | Primary and secondary outcomes |
| Key findings | Main results with effect sizes |
| Limitations | Author-reported limitations |
| Quality score | From assessment tool |
| Funding | Funding sources |
| Conflicts | Declared conflicts of interest |

**Steps:**
1. Create extraction spreadsheet with all fields
2. Pilot extraction on 3-5 papers
3. Refine extraction form based on pilot
4. Extract data from all included studies
5. Verify extraction accuracy (dual extraction for ≥20%)

**Exit Criteria:**
- [ ] Extraction form finalized after pilot
- [ ] Data extracted from 100% of included studies
- [ ] Quality assessment completed for all studies
- [ ] Extraction verified (dual extraction ≥20%)
- [ ] Missing data documented as "NR" (not reported)

---

### Stage 5: Quality Assessment

**Objective:** Evaluate risk of bias and study quality

**Steps:**
1. Select appropriate tool (see Quality Assessment Tool Selection matrix)
2. Assess each study independently
3. Create summary table:

| Study | Domain 1 | Domain 2 | Domain 3 | Overall |
|-------|---------|---------|---------|---------|
| Smith 2020 | Low | Low | High | Moderate |
| Jones 2021 | Low | Low | Low | Low |
| ... | ... | ... | ... | ... |

4. Create risk of bias visualization (traffic light plot)
5. Consider excluding very low quality studies

**Quality Rating Thresholds:**

| Overall Rating | Criteria |
|---------------|----------|
| High quality (Low RoB) | All domains low or some concerns |
| Moderate quality | One high-risk domain, others low |
| Low quality (High RoB) | ≥2 high-risk domains |
| Very low quality | ≥3 high-risk or critical flaw |

**Exit Criteria:**
- [ ] All studies assessed with selected tool
- [ ] Summary table complete
- [ ] Risk of bias plot generated
- [ ] Quality ratings assigned (High/Moderate/Low/Very Low)
- [ ] Decision made on study exclusions (if any)

---

### Stage 6: Synthesis and Analysis

**Objective:** Synthesize findings thematically across studies

**Synthesis Approaches:**

| Data Type | Approach | Output |
|-----------|----------|--------|
| Quantitative, similar outcomes | Meta-analysis | Forest plot, pooled effect |
| Quantitative, diverse outcomes | Narrative synthesis | Comparison tables |
| Qualitative | Thematic analysis | Theme framework |
| Mixed | Integrated synthesis | Mixed methods matrix |

**Steps:**
1. Identify 3-7 major themes from extracted data
2. Group studies by theme (studies may appear in multiple themes)
3. For each theme, synthesize across studies:
   - Consensus areas (≥80% agreement)
   - Conflicting findings (document discrepancies)
   - Gaps (what's not studied)

4. Write thematic synthesis (NOT study-by-study):
   ```markdown
   #### Theme 1: [Theme Name]

   [Synthesis paragraph integrating findings from multiple studies]
   Multiple studies (n=X) found that... [citations].
   However, conflicting evidence exists regarding...
   [citations]. A notable gap is...
   ```

   **Citation Density Requirements** (see Literature Coverage Thresholds):
   - Review articles: 4-6 citations per paragraph
   - Systematic reviews: All relevant studies cited per theme
   - Original research: 2-4 citations per paragraph

5. Create summary tables per theme
6. If meta-analysis: calculate pooled effects, heterogeneity (I²)

**Exit Criteria:**
- [ ] 3-7 themes identified and documented
- [ ] All studies mapped to themes
- [ ] Thematic synthesis written (not study-by-study)
- [ ] Consensus/conflict/gaps identified per theme
- [ ] Summary tables created
- [ ] Meta-analysis complete (if applicable)

---

### Stage 7: Citation Verification

**Objective:** Verify accuracy of all citations

**Steps:**
1. Extract all DOIs and citations from document
2. Verify each DOI resolves correctly
3. Cross-check metadata (authors, title, year, journal)
4. Format citations consistently (see Citation Style Matrix)

**Citation Verification Checklist:**

- [ ] All DOIs tested and resolve
- [ ] Author names match source
- [ ] Publication year correct
- [ ] Journal name correct (abbreviated or full per style)
- [ ] Page numbers/article IDs correct
- [ ] In-text citations match reference list

**Citation Style Quick Reference:**

| Style | In-text | Reference Format |
|-------|---------|------------------|
| APA 7 | (Smith et al., 2023) | Smith, J. D., Johnson, M. L., & Williams, K. (2023). Title. *Journal*, *22*(4), 301-318. |
| Nature | Superscript^1,2^ | 1. Smith, J. D. et al. Title. *J. Name* **22**, 301-318 (2023). |
| Vancouver | Superscript^1,2^ | 1. Smith JD, Johnson ML. Title. J Name. 2023;22(4):301-18. |
| IEEE | Bracketed [1], [2] | [1] J. D. Smith et al., "Title," *J. Name*, vol. 22, no. 4, pp. 301-318, 2023. |

**Exit Criteria:**
- [ ] 100% of DOIs verified
- [ ] All metadata cross-checked
- [ ] Citation style consistent throughout
- [ ] No [?] or broken references
- [ ] Reference list formatted correctly

---

### Stage 8: Document Generation

**Objective:** Produce final review document

**Steps:**
1. Compile review from template:
   - Abstract (structured: Background, Methods, Results, Conclusions)
   - Introduction (context, rationale, objectives)
   - Methods (search strategy, selection, assessment)
   - Results (PRISMA, synthesis by theme, quality summary)
   - Discussion (interpretation, limitations, implications)
   - Conclusions (key findings, future directions)

2. Generate PDF with citations formatted
3. Review and proofread

**Document Quality Checklist:**

- [ ] Abstract covers all sections (≤350 words for structured)
- [ ] PRISMA flow diagram included
- [ ] Search strategy fully reproducible
- [ ] Results organized thematically (not study-by-study)
- [ ] Quality assessment summarized
- [ ] Limitations section included
- [ ] All references verified
- [ ] PDF generates without errors

**Exit Criteria:**
- [ ] All sections complete
- [ ] PRISMA diagram included
- [ ] Citations verified and formatted
- [ ] PDF generated successfully
- [ ] Proofread for errors

</workflow>

<success_criteria>
## Success Criteria

**Quantitative Thresholds:**

| Metric | Minimum | Target | Systematic Review |
|--------|---------|--------|-------------------|
| Databases searched | 2 | 3 | ≥4 |
| Search documentation | Dates only | Full queries | Reproducible queries |
| Included studies | ≥10 | ≥20 | ≥30 |
| Quality assessment | Informal | Standardized tool | Dual assessment |
| Citation verification | Spot-check | 100% DOIs | 100% + metadata |
| PRISMA compliance | None | Partial | Full checklist |

**Completion Checklist:**
- [ ] Research question in structured format (PICO/PEO/SPIDER)
- [ ] ≥3 databases searched with documented queries
- [ ] PRISMA flow diagram with complete numbers
- [ ] Quality assessment using standardized tool
- [ ] Thematic synthesis (not study-by-study)
- [ ] 100% of DOIs verified and resolving
- [ ] PDF generates without citation errors
- [ ] All PRISMA checklist items addressed (for systematic reviews)

</success_criteria>

<scope>
## Scope

**In Scope:**
- Systematic, scoping, narrative, and rapid reviews
- Multi-database literature searching
- Screening and selection workflows
- Quality assessment guidance
- Thematic synthesis methodology
- Citation verification
- PRISMA-compliant documentation
- PDF generation

**Out of Scope** (use specialized resources):
- Meta-analysis statistical calculations → use `statistical-analysis`
- Single paper lookups → use `research-lookup`
- Citation formatting only → use `citation-management`
- Figure creation → use `plotting-libraries` or `scientific-schematics`
- Full-text retrieval → manual or institutional access

</scope>

<anti_patterns>
## Common Pitfalls

### 1. Single Database Search

**Anti-pattern:**
```
Searched only PubMed for a systematic review.
Missed 40% of relevant papers indexed only
in Scopus, Web of Science, or preprint servers.
```

**Solution:**
```
Minimum 3 databases for systematic reviews.
Select based on domain (see Database Selection Matrix).
Include preprint servers for recent findings.
Document all searches for reproducibility.
```

---

### 2. Study-by-Study Summaries

**Anti-pattern:**
```
Results section:
"Smith (2020) found X. Jones (2021) found Y.
Brown (2022) found Z..."
No synthesis, just serial summaries.
```

**Solution:**
```
Organize by themes, synthesize across studies:
"Multiple studies (n=8) found consistent evidence
that [finding], though effect sizes varied (d=0.3-0.8).
Conflicting results emerged regarding [topic]..."
```

---

### 3. Undocumented Search Strategy

**Anti-pattern:**
```
Methods: "We searched PubMed and Google Scholar
for relevant articles."
Not reproducible. No dates, no search terms.
```

**Solution:**
```
Document verbatim for each database:
- Database name
- Date searched
- Exact search string used
- Number of results
- Filters applied (date range, language, etc.)
```

---

### 4. No Quality Assessment

**Anti-pattern:**
```
Included all studies regardless of quality.
Weighted a flawed pilot study equally with
a rigorous multi-center RCT.
```

**Solution:**
```
Use standardized quality assessment tool:
- RCTs: Cochrane RoB 2.0
- Observational: Newcastle-Ottawa Scale
- Systematic reviews: AMSTAR 2
Report quality in results, consider in synthesis.
```

---

### 5. Unverified Citations

**Anti-pattern:**
```
Copy-pasted citations from Google Scholar
without verification. Final document has
3 broken DOIs and 2 wrong publication years.
```

**Solution:**
```
Verify 100% of citations:
1. Test all DOIs resolve
2. Cross-check author names, titles, years
3. Use reference manager for consistency
4. Run verification script before submission
```

---

### 6. Missing Exclusion Reasons

**Anti-pattern:**
```
PRISMA diagram shows "312 excluded" at
abstract screening with no breakdown.
Not reproducible or auditable.
```

**Solution:**
```
Track exclusion reasons at each stage:
- Wrong population (n=45)
- Wrong study design (n=82)
- Wrong outcome (n=31)
- Not English (n=12)
- Conference abstract only (n=142)
```

</anti_patterns>

<templates>
## Output Templates

### PRISMA Flow Diagram Template

```
┌─────────────────────────────────────────────────────────┐
│                     IDENTIFICATION                       │
├─────────────────────────────────────────────────────────┤
│  Records identified through          Additional records  │
│  database searching                  from other sources  │
│  (n = {{db_total}})                  (n = {{other}})     │
└─────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────┐
│                      SCREENING                           │
├─────────────────────────────────────────────────────────┤
│  Records after duplicates removed                        │
│  (n = {{after_dedup}})                                   │
│                           │                              │
│                           ▼                              │
│  Records screened ──────────────► Records excluded       │
│  (n = {{screened}})                (n = {{title_excl}})  │
│                           │                              │
│                           ▼                              │
│  Full-text assessed ────────────► Excluded with reasons  │
│  (n = {{ft_assessed}})            (n = {{ft_excl}})      │
│                                   • Reason 1 (n = )      │
│                                   • Reason 2 (n = )      │
│                                   • Reason 3 (n = )      │
└─────────────────────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────┐
│                       INCLUDED                           │
├─────────────────────────────────────────────────────────┤
│  Studies included in synthesis                           │
│  (n = {{included}})                                      │
│  • Quantitative synthesis (n = {{quant}})                │
│  • Qualitative synthesis (n = {{qual}})                  │
└─────────────────────────────────────────────────────────┘
```

### Search Documentation Template

```markdown
## Search Strategy

### Database: {{database_name}}

**Date searched:** {{YYYY-MM-DD}}

**Search query:**
```
{{full_search_string}}
```

**Filters applied:**
- Date range: {{start}} to {{end}}
- Language: {{language}}
- Publication type: {{types}}

**Results:** {{N}} records

---

### Database: {{database_name_2}}

[Repeat for each database]
```

### Data Extraction Template

```markdown
## Study: {{First Author}} {{Year}}

| Field | Value |
|-------|-------|
| Full citation | {{citation}} |
| Country | {{country}} |
| Study design | {{design}} |
| Sample size | n = {{N}} |
| Population | {{demographics}} |
| Intervention | {{intervention}} |
| Comparator | {{comparator}} |
| Primary outcome | {{outcome}} |
| Key finding | {{finding with effect size}} |
| Quality rating | {{High/Moderate/Low}} |
| Limitations | {{limitations}} |
| Funding | {{funding}} |
```

</templates>

<cross_references>
## Related Skills

| Skill | Relationship |
|-------|-------------|
| `research-lookup` | Use for single paper lookups and quick searches |
| `citation-management` | Use for citation formatting after extraction |
| `scientific-schematics` | Generate PRISMA diagrams and synthesis figures |
| `plotting-libraries` | Create forest plots and summary visualizations |
| `statistical-analysis` | Perform meta-analysis calculations |
| `scientific-writing` | Write manuscript sections from review |
| `venue-templates` | Format for specific journal requirements |

**Skill Selection:**
See `SKILL_ROUTER.md` for decision trees when multiple skills may apply.

</cross_references>

<references>
## Reference Documents

| Document | Purpose |
|----------|---------|
| `references/citation_styles.md` | Detailed formatting for APA, Nature, Vancouver, IEEE |
| `references/database_strategies.md` | Database-specific search optimization |
| `references/quality_assessment_tools.md` | Tool selection and scoring guidance |
| `references/prisma_checklist.md` | Full PRISMA 2020 checklist |

## External Resources

| Resource | URL |
|----------|-----|
| PRISMA Statement | http://www.prisma-statement.org/ |
| Cochrane Handbook | https://training.cochrane.org/handbook |
| PROSPERO Registry | https://www.crd.york.ac.uk/prospero/ |
| MeSH Browser | https://meshb.nlm.nih.gov/search |
| AMSTAR 2 | https://amstar.ca/ |

</references>
