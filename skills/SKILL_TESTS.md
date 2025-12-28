# Skill Test Scenarios

> Phase 5 of Skills Determinism Audit: Validation
> Created: 2025-12-28
> Purpose: Verify skill determinism through reproducible test scenarios

---

## Overview

This document defines test scenarios for each skill in the library. Each skill has:
- **Scenario 1: Basic Usage** — Typical use case with expected deterministic behavior
- **Scenario 2: Edge Case** — Boundary condition testing skill routing and error handling

**Test Execution:** Run scenarios manually and verify outputs match expected behavior.

---

## statistical-analysis

### Scenario 1: Basic Usage

**Input:** "I have gene expression data from 3 treatment groups (Control, Drug A, Drug B) with 6 biological replicates each. Data appears normally distributed. Which test should I use to compare groups?"

**Expected Behavior:**
1. Skill identifies: 3 groups, continuous data, normal distribution, unpaired
2. Recommends: One-way ANOVA with Tukey HSD post-hoc
3. Provides: Code snippet in Python (scipy) or R
4. Notes: Effect size (η²) and power considerations

**Success Criteria:**
- [ ] Test selection matches decision matrix (One-way ANOVA)
- [ ] Post-hoc test specified (Tukey HSD)
- [ ] Code provided is syntactically correct
- [ ] Reporting template referenced

---

### Scenario 2: Edge Case

**Input:** "I tested 15,000 genes for differential expression between 2 conditions using t-tests. I found 800 genes with p < 0.05. Are these all significant?"

**Expected Behavior:**
1. Skill identifies: Multiple testing problem (15,000 tests)
2. Calculates expected false positives: ~750 at α=0.05
3. Recommends: Benjamini-Hochberg FDR correction
4. Explains: 800 raw hits likely includes ~750 false positives
5. Provides: Code for FDR correction

**Success Criteria:**
- [ ] Multiple testing issue identified immediately
- [ ] Expected false positive count stated (~750)
- [ ] FDR correction recommended (not Bonferroni for discovery)
- [ ] Pitfall #1 (Multiple Testing Without Correction) explicitly referenced

---

## scientific-writing

### Scenario 1: Basic Usage

**Input:** "I need to write the Methods section for a randomized controlled trial comparing two drug treatments. The study had 120 participants, 60 per arm."

**Expected Behavior:**
1. Skill applies IMRAD structure guidance
2. References CONSORT reporting guideline
3. Provides Methods section template with subsections:
   - Study design, Participants, Interventions, Outcomes, Sample size, Randomization, Blinding, Statistical analysis
4. Suggests word count target: 900-1400 words (standard paper)
5. Recommends citation count: 3-8 for Methods

**Success Criteria:**
- [ ] CONSORT guideline explicitly referenced
- [ ] All required Methods subsections listed
- [ ] Word count guidance provided
- [ ] Two-stage process mentioned (outline → prose)

---

### Scenario 2: Edge Case

**Input:** "Here's my Results section draft with bullet points:
• Treatment A showed 45% improvement
• Treatment B showed 32% improvement
• p-value was 0.03
• Side effects were minimal"

**Expected Behavior:**
1. Skill identifies: Bullet points not acceptable for final manuscript
2. References: "Critical Principle: Always write in full paragraphs"
3. Provides: Conversion to prose example
4. Notes missing elements: effect sizes, confidence intervals, exact statistics (t, df)
5. Suggests proper reporting format from templates

**Success Criteria:**
- [ ] Bullet points flagged as inappropriate for submission
- [ ] Prose conversion demonstrated
- [ ] Missing statistical details identified (CI, effect size, test statistic)
- [ ] Results Reporting Format template referenced

---

## generate-image

### Scenario 1: Basic Usage

**Input:** "I need a scientific illustration of a DNA double helix with a highlighted mutation site for my paper on BRCA1 variants."

**Expected Behavior:**
1. Skill confirms: This is generate-image scope (illustration, not diagram)
2. Constructs prompt using template: [SUBJECT] + [STYLE] + [MEDIUM] + [TECHNICAL] + [QUALITY]
3. Suggests: "DNA double helix structure with highlighted mutation site, scientific illustration style, clean white background, high detail, professional molecular visualization, suitable for academic publication"
4. Recommends model: gemini-3-pro-image-preview or flux.2-pro
5. Provides command with output path

**Success Criteria:**
- [ ] Prompt follows 5-component structure
- [ ] Model selection justified from matrix
- [ ] Command syntax correct
- [ ] Quality validation checklist referenced

---

### Scenario 2: Edge Case

**Input:** "Generate a flowchart showing the CRISPR gene editing process with boxes and arrows."

**Expected Behavior:**
1. Skill identifies: This is NOT generate-image scope
2. Routes to: `scientific-schematics` skill
3. Explains: Flowcharts require structured diagrams, not AI image generation
4. Decision tree path: Technical diagram? → Flowchart → scientific-schematics

**Success Criteria:**
- [ ] Skill correctly declines the request
- [ ] Routes to scientific-schematics explicitly
- [ ] Decision tree logic explained
- [ ] No image generation attempted

---

## literature-review

### Scenario 1: Basic Usage

**Input:** "I want to conduct a systematic review on the effectiveness of immunotherapy for melanoma. I plan to search PubMed only."

**Expected Behavior:**
1. Skill identifies: Systematic review requires ≥3 databases
2. Recommends databases from matrix: PubMed + Cochrane + Embase + Scopus
3. Provides PICO framework:
   - P: Melanoma patients
   - I: Immunotherapy
   - C: Standard treatment/placebo
   - O: Effectiveness (survival, response rate)
4. References coverage thresholds: ≥50 papers, ≥3 databases
5. Provides Stage 1 (Planning) checklist

**Success Criteria:**
- [ ] Single database flagged as insufficient
- [ ] ≥3 databases recommended
- [ ] PICO framework applied
- [ ] Coverage thresholds from QUANTIFICATION_THRESHOLDS.md cited

---

### Scenario 2: Edge Case

**Input:** "I've found 25 papers on my topic. Here's my synthesis:
- Smith 2020 found X
- Jones 2021 found Y
- Brown 2022 found Z
- [continues for all 25 papers]"

**Expected Behavior:**
1. Skill identifies: Study-by-study summary anti-pattern
2. References: Anti-pattern #2 "Study-by-Study Summaries"
3. Provides thematic synthesis guidance:
   - Group by 3-7 themes
   - Synthesize across studies per theme
   - Note consensus, conflicts, gaps
4. Shows correct format: "Multiple studies (n=X) found..."
5. References citation density: 4-6 citations per paragraph for reviews

**Success Criteria:**
- [ ] Anti-pattern #2 explicitly identified
- [ ] Thematic synthesis structure explained
- [ ] Example of correct synthesis format provided
- [ ] Citation density requirements stated

---

## citation-management

### Scenario 1: Basic Usage

**Input:** "I have a DOI: 10.1038/s41586-024-07560-8. Convert it to BibTeX format for my manuscript."

**Expected Behavior:**
1. Skill identifies DOI format
2. Queries CrossRef API as primary source
3. Generates BibTeX with citation key format: AuthorYearkeyword
4. Includes all required fields: author, title, journal, year, volume, pages, doi
5. Applies title case protection for proper nouns

**Success Criteria:**
- [ ] DOI resolution verified via doi.org
- [ ] BibTeX entry type correct (@article)
- [ ] Citation key follows AuthorYearkeyword format
- [ ] All required fields present per entry type table

---

### Scenario 2: Edge Case

**Input:** "I'm citing this arXiv preprint: 2401.12345. But I think it was published in a journal now."

**Expected Behavior:**
1. Skill identifies: arXiv ID format
2. Recognizes potential preprint → published transition
3. Anti-pattern #5 referenced: "Citing Outdated Preprints"
4. Suggests checking for published version via CrossRef or preprint page
5. If published version exists, recommends updating to journal citation

**Success Criteria:**
- [ ] Preprint update check recommended
- [ ] Anti-pattern #5 explicitly referenced
- [ ] Both @misc (preprint) and @article (published) formats shown
- [ ] DOI for published version requested if available

---

## code-documentation

### Scenario 1: Basic Usage

**Input:** "Write NumPy-style docstrings for these Python functions in my RNA-seq analysis script."

**Expected Behavior:**
1. Skill identifies: Scientific Python context → NumPy style
2. Uses style selection decision tree
3. Generates docstrings with: Summary, Parameters, Returns, Examples
4. Includes type hints where applicable
5. Provides runnable example in docstring

**Success Criteria:**
- [ ] NumPy style selected (not Google) for scientific Python
- [ ] All four sections present: Summary, Parameters, Returns, Examples
- [ ] Example code is syntactically correct
- [ ] Documentation hierarchy level 2 (docstrings) applied

---

### Scenario 2: Edge Case

**Input:** "Here's my code with comments:
```python
i += 1  # increment i
data = load_data()  # load the data
```
Review my documentation."

**Expected Behavior:**
1. Skill identifies: Over-commenting anti-pattern
2. References: "Anti-pattern: Commenting obvious code"
3. Explains comments should explain "why" not "what"
4. Suggests removing self-evident comments
5. Recommends adding docstrings to functions instead

**Success Criteria:**
- [ ] Over-commenting anti-pattern identified
- [ ] "Why not what" principle explained
- [ ] Specific examples flagged (`# increment i`)
- [ ] Documentation hierarchy referenced (docstrings > inline comments)

---

## hypothesis-generation

### Scenario 1: Basic Usage

**Input:** "We observed that mice fed a high-fat diet develop cognitive impairment. I need to generate hypotheses to explain this phenomenon."

**Expected Behavior:**
1. Skill follows 8-stage workflow
2. Stage 2: Searches ≥3 databases, retrieves ≥40 papers
3. Stage 4: Generates 3-5 competing hypotheses with different mechanisms
4. Stage 5: Scores each hypothesis on 10-point quality rubric
5. Stage 6: Provides experimental design for each viable hypothesis
6. References QUANTIFICATION_THRESHOLDS.md for literature coverage

**Success Criteria:**
- [ ] 3-5 distinct hypotheses generated
- [ ] Each hypothesis scored on 5 criteria (testability, falsifiability, parsimony, explanatory power, distinctiveness)
- [ ] ≥1 hypothesis scores ≥8/10 (strong)
- [ ] Experimental design includes controls and sample size rationale

---

### Scenario 2: Edge Case

**Input:** "Here's my hypothesis: 'Something is causing the cognitive impairment.' Please evaluate it."

**Expected Behavior:**
1. Skill identifies: Hypothesis lacks testability and mechanistic explanation
2. Quality rubric applied: Likely scores 1-2/10
3. Explains hypothesis must be mechanistic (explain HOW and WHY)
4. References hypothesis quality criteria
5. Offers to help reformulate into testable hypothesis

**Success Criteria:**
- [ ] Low quality score assigned (<5/10)
- [ ] Testability and mechanistic criteria explicitly failed
- [ ] Quality rubric criteria breakdown provided
- [ ] Weak hypothesis appropriately deprioritized

---

## latex-posters

### Scenario 1: Basic Usage

**Input:** "I need to create an A0 portrait poster for a neuroscience conference. I'm familiar with Beamer."

**Expected Behavior:**
1. Skill uses decision framework: Knows Beamer → beamerposter
2. Confirms A0 portrait dimensions (841 × 1189 mm)
3. Mandates 2-3 AI-generated figures via scientific-schematics
4. Targets 40-50% visual content, 300-800 words
5. Provides beamerposter setup code with correct page size

**Success Criteria:**
- [ ] beamerposter package selected (not tikzposter/baposter)
- [ ] Poster size configured correctly (A0, portrait)
- [ ] 2-3 figures requirement stated
- [ ] Word count target (300-800) provided

---

### Scenario 2: Edge Case

**Input:** "Here's my poster draft. It has 1500 words, full paragraphs, and 1 small figure."

**Expected Behavior:**
1. Skill identifies: Text-heavy poster anti-pattern
2. Flags: >1000 words, paragraphs instead of bullets, <30% visual content
3. References: Anti-pattern #1 "Text-Heavy Poster"
4. Recommends: Cut to 300-800 words, use bullet points, add 2-3 AI schematics
5. Cites visual content threshold: 40-50%

**Success Criteria:**
- [ ] Text-heavy anti-pattern identified
- [ ] Word count (1500) flagged as exceeding 800 max
- [ ] Full paragraphs flagged (bullets required)
- [ ] Missing figures identified (need 2-3, has 1)

---

## markitdown

### Scenario 1: Basic Usage

**Input:** "Convert this PDF research paper to Markdown for processing: paper.pdf"

**Expected Behavior:**
1. Skill identifies: PDF input type
2. Uses decision tree: Simple text PDF → Basic MarkItDown
3. Executes: `markitdown paper.pdf -o output.md`
4. Validates output: non-empty, valid Markdown, tables preserved
5. Adds metadata header with source and conversion date

**Success Criteria:**
- [ ] Correct conversion method selected from decision tree
- [ ] markitdown command syntax correct
- [ ] Output validation performed (>100 characters)
- [ ] Markdown structure verified

---

### Scenario 2: Edge Case

**Input:** "Convert this scanned PDF to Markdown. I don't have tesseract installed."

**Expected Behavior:**
1. Skill identifies: Scanned PDF requires OCR
2. Checks dependencies: tesseract missing
3. References: Anti-pattern #1 "Missing Dependencies"
4. Provides installation command before conversion
5. Explains fallback options (Azure Document Intelligence)

**Success Criteria:**
- [ ] OCR requirement identified for scanned PDF
- [ ] Missing dependency detected
- [ ] Installation instructions provided
- [ ] Alternative (Azure DocIntel) mentioned

---

## paper-2-web

### Scenario 1: Basic Usage

**Input:** "Convert my LaTeX manuscript to an interactive HTML page for sharing."

**Expected Behavior:**
1. Skill identifies: LaTeX → HTML conversion
2. Selects appropriate conversion pipeline
3. Preserves equations, figures, and references
4. Generates responsive, accessible HTML
5. Includes proper citation formatting

**Success Criteria:**
- [ ] LaTeX source correctly processed
- [ ] Equations render properly (MathJax/KaTeX)
- [ ] Figures embedded with captions
- [ ] HTML passes basic accessibility checks

---

### Scenario 2: Edge Case

**Input:** "I want to convert this paper to a website, but I also need to extract its figures separately."

**Expected Behavior:**
1. Skill identifies: Multi-output requirement
2. Differentiates from markitdown (input→markdown) - this is output direction
3. Routes figure extraction to appropriate tool
4. Clarifies: paper-2-web for web output, scientific-schematics for creating new figures

**Success Criteria:**
- [ ] Conversion vs. extraction differentiated
- [ ] Correct skill routing explained
- [ ] Multi-step workflow suggested if needed

---

## peer-review

### Scenario 1: Basic Usage

**Input:** "Review this manuscript for a Nature-style journal. Focus on methodology and statistical rigor."

**Expected Behavior:**
1. Skill applies 7-stage workflow
2. Stage 3: Checks sample sizes against replication thresholds (≥3 biological replicates)
3. References QUANTIFICATION_THRESHOLDS.md for severity classification
4. Uses venue-templates `reviewer_expectations.md` for Nature standards
5. Structures report: Summary, Major Comments, Minor Comments

**Success Criteria:**
- [ ] All 7 review stages applied
- [ ] Sample size thresholds checked (n≥3 biological, n≥3 technical for in vitro)
- [ ] Issues classified by severity (Critical/Major/Minor)
- [ ] Constructive, actionable recommendations provided

---

### Scenario 2: Edge Case

**Input:** "Review this PowerPoint presentation (slides.pdf) for a conference talk."

**Expected Behavior:**
1. Skill identifies: Presentation PDF → must convert to images first
2. CRITICAL: Never reads PDF directly for presentations
3. Uses: `python skills/scientific-slides/scripts/pdf_to_images.py slides.pdf review/`
4. Reviews each slide image for layout, overflow, contrast, font size
5. Documents issues by slide number

**Success Criteria:**
- [ ] PDF-to-images conversion mandated (not direct PDF reading)
- [ ] Image-based review workflow followed
- [ ] Presentation-specific criteria applied (font ≥18pt, contrast ≥4.5:1)
- [ ] Issues documented by slide number

---

## plotting-libraries

### Scenario 1: Basic Usage

**Input:** "Create a publication-quality volcano plot of my differential expression results using Python."

**Expected Behavior:**
1. Skill identifies: Volcano plot = matplotlib/seaborn territory
2. Provides code pattern with proper styling
3. Applies publication standards: ≥300 DPI, proper font sizes, colorblind-safe
4. Includes significance thresholds and fold-change cutoffs
5. References visual-design for styling principles

**Success Criteria:**
- [ ] Appropriate library selected (matplotlib/seaborn)
- [ ] Code syntax is correct and runnable
- [ ] Publication standards applied (resolution, fonts)
- [ ] Colorblind-safe palette used

---

### Scenario 2: Edge Case

**Input:** "Create a neural network architecture diagram using matplotlib."

**Expected Behavior:**
1. Skill identifies: This is NOT a data visualization
2. Routes to: `scientific-schematics` skill
3. Explains: matplotlib is for plots/charts, not diagrams
4. Decision tree: Neural network architecture → scientific-schematics

**Success Criteria:**
- [ ] Request correctly declined for plotting-libraries scope
- [ ] Routes to scientific-schematics explicitly
- [ ] Scope boundary explained (plots vs. diagrams)

---

## pptx-posters

### Scenario 1: Basic Usage

**Input:** "Create a 36×48 inch research poster in PowerPoint for a chemistry conference."

**Expected Behavior:**
1. Skill identifies: PowerPoint poster format
2. Confirms 36×48" dimensions (US standard)
3. Similar to latex-posters: 40-50% visual content, 300-800 words
4. Mandates 2-3 figures via scientific-schematics
5. Provides PowerPoint-specific layout guidance

**Success Criteria:**
- [ ] PowerPoint format confirmed (not LaTeX)
- [ ] Poster size correctly configured
- [ ] Visual content requirements stated
- [ ] Figure generation requirement included

---

### Scenario 2: Edge Case

**Input:** "I need version control and equations in my poster."

**Expected Behavior:**
1. Skill identifies: Version control and equations favor LaTeX
2. Routes to: `latex-posters` skill
3. Uses decision tree: Need equations/version control → LaTeX posters
4. Explains trade-offs between PowerPoint and LaTeX approaches

**Success Criteria:**
- [ ] LaTeX route recommended for equations/version control
- [ ] Decision tree logic explained
- [ ] Trade-offs articulated

---

## reproducible-research

### Scenario 1: Basic Usage

**Input:** "I'm preparing my RNA-seq project for publication. What do I need for reproducibility?"

**Expected Behavior:**
1. Skill applies FAIR principles
2. Recommends: environment.yml with pinned versions (≥75% checklist)
3. Guides data deposition: SRA for FASTQ, GEO for processed
4. Provides Data Availability statement template
5. References reproducibility checklist (20 items, need ≥15)

**Success Criteria:**
- [ ] FAIR principles referenced
- [ ] Correct repositories recommended (SRA/GEO for RNA-seq)
- [ ] Environment specification level determined (Level 2+)
- [ ] Reproducibility score threshold stated (≥75%)

---

### Scenario 2: Edge Case

**Input:** "Here's my requirements.txt:
```
numpy
pandas
scipy
```
Is this good enough for reproducibility?"

**Expected Behavior:**
1. Skill identifies: Unpinned dependencies (critical error)
2. References: Common pitfall "Unpinned dependencies" (-20% score)
3. Explains: Level 1 (Minimal) = 25-49% reproducibility score
4. Recommends: Pin all versions with `==`, upgrade to environment.yml
5. Target: Level 2 (Standard) for publication-ready

**Success Criteria:**
- [ ] Unpinned dependencies flagged as critical
- [ ] Reproducibility level identified (Level 1 = 25-49%)
- [ ] Specific fix provided (pin versions with ==)
- [ ] Target level stated (Level 2+ for publication)

---

## research-lookup

### Scenario 1: Basic Usage

**Input:** "Find recent papers on CRISPR base editing for sickle cell disease treatment."

**Expected Behavior:**
1. Skill uses automatic model selection
2. Query classified: Straightforward lookup → Sonar Pro Search
3. Returns: Summary, 3-5 relevant papers with citations
4. Includes: DOIs, publication years, key findings
5. Prioritizes peer-reviewed sources from past 2-3 years

**Success Criteria:**
- [ ] Sonar Pro Search model selected (not Reasoning Pro)
- [ ] 3-5 papers returned with complete citations
- [ ] Recent papers prioritized (2022-2024)
- [ ] DOIs included for verification

---

### Scenario 2: Edge Case

**Input:** "Compare and contrast the mechanisms of base editing versus prime editing for therapeutic applications."

**Expected Behavior:**
1. Skill identifies: Complex analytical query (compare, contrast, mechanisms)
2. Triggers Sonar Reasoning Pro (reasoning keywords detected)
3. Provides detailed comparative analysis
4. Synthesizes across multiple studies
5. Discusses trade-offs and ongoing debates

**Success Criteria:**
- [ ] Sonar Reasoning Pro model selected (complexity ≥3 points)
- [ ] Reasoning keywords detected (compare, contrast, mechanisms)
- [ ] Comparative analysis provided (not just listing)
- [ ] Trade-offs and nuances included

---

## scholar-evaluation

### Scenario 1: Basic Usage

**Input:** "Evaluate this research manuscript for publication readiness using the ScholarEval framework."

**Expected Behavior:**
1. Skill applies 8-dimension framework
2. Scores each dimension on 1-5 scale with evidence
3. Calculates weighted aggregate score
4. Assigns quality tier (A-F)
5. Provides prioritized recommendations (P1-P4)

**Success Criteria:**
- [ ] All 8 dimensions assessed
- [ ] Scores include specific evidence with page references
- [ ] Weighted aggregate calculated correctly
- [ ] Quality tier assigned based on score ranges

---

### Scenario 2: Edge Case

**Input:** "Evaluate this paper. The methodology score is 3."

**Expected Behavior:**
1. Skill identifies: Subjective scoring anti-pattern
2. References: Anti-pattern #1 "Subjective Scoring"
3. Requires: Specific evidence, strengths, weaknesses for each score
4. Demonstrates proper format with page/section references

**Success Criteria:**
- [ ] Subjective scoring anti-pattern identified
- [ ] Evidence requirement emphasized
- [ ] Proper scoring format demonstrated
- [ ] Strengths AND weaknesses required (balanced feedback)

---

## peer-review: Critical Analysis Framework

> **Note:** These scenarios test the Critical Analysis Framework section of peer-review, which absorbed the former scientific-critical-thinking skill.

### Scenario 1: Basic Usage (Critical Analysis)

**Input:** "Evaluate this claim: 'Our treatment showed a 50% improvement (p=0.049) in patient outcomes.'"

**Expected Behavior:**
1. Skill applies Critical Analysis Framework (GRADE, Claim Evaluation, Bias Detection)
2. Identifies: p-value barely significant (0.049 ≈ 0.05)
3. Questions: Effect size? Confidence interval? Multiple comparisons?
4. Notes potential issues: p-hacking risk, clinical vs. statistical significance
5. Suggests what additional information is needed

**Success Criteria:**
- [ ] Borderline p-value flagged
- [ ] Missing information identified (CI, effect size)
- [ ] Statistical vs. clinical significance differentiated
- [ ] Potential biases/issues raised

---

### Scenario 2: Edge Case (Critical Analysis)

**Input:** "This paper claims their drug cures cancer based on a study with 5 patients."

**Expected Behavior:**
1. Skill identifies: Multiple critical issues via Critical Analysis Framework
2. Flags: Small sample size (n=5), overstated conclusions ("cures")
3. References: Correlation vs. causation, replication requirements, logical fallacies
4. Notes: Cannot establish causation from single small study
5. Recommends caution in interpretation

**Success Criteria:**
- [ ] Sample size inadequacy identified
- [ ] Overstated claim flagged ("cures")
- [ ] Causal inference limitations explained
- [ ] Need for replication emphasized

---

## scientific-schematics

### Scenario 1: Basic Usage

**Input:** "Create a CONSORT participant flow diagram for my clinical trial with 500 screened, 150 excluded, 350 randomized into 2 arms."

**Expected Behavior:**
1. Skill identifies: CONSORT flowchart type
2. Uses decision matrix: Study Design → Clinical trial → CONSORT flowchart
3. Constructs detailed prompt with all numbers and exclusion reasons
4. Sets quality threshold based on document type
5. Generates via Nano Banana Pro with Gemini 3 Pro review

**Success Criteria:**
- [ ] CONSORT flowchart type selected from decision matrix
- [ ] Prompt includes specific numbers (500, 150, 350)
- [ ] Command syntax correct with doc-type flag
- [ ] Quality threshold appropriate for document type

---

### Scenario 2: Edge Case

**Input:** "Generate a bar chart showing my gene expression results."

**Expected Behavior:**
1. Skill identifies: This is NOT a schematic/diagram
2. Bar chart = data visualization → `plotting-libraries` scope
3. Routes to plotting-libraries skill
4. Explains scope boundary: schematics = diagrams, not plots

**Success Criteria:**
- [ ] Request correctly declined for scientific-schematics scope
- [ ] Routes to plotting-libraries explicitly
- [ ] Scope boundary explained (diagrams vs. data visualizations)

---

## scientific-slides

### Scenario 1: Basic Usage

**Input:** "Create a 15-minute conference presentation for my machine learning paper."

**Expected Behavior:**
1. Skill calculates: ~15 slides for 15 minutes (~1 slide/minute)
2. Structures: Title → Background (3-4) → Methods (3-4) → Results (5-6) → Conclusions (1-2)
3. Mandates: 40-50% visual content
4. Font requirements: ≥24pt body, ≥32pt headers
5. References scientific-schematics for diagram generation

**Success Criteria:**
- [ ] Slide count matches duration (~15 slides)
- [ ] Standard structure applied (intro/methods/results/conclusions)
- [ ] Font size requirements stated (≥24pt body)
- [ ] Visual content requirement included

---

### Scenario 2: Edge Case

**Input:** "Here's my talk with 45 slides for a 15-minute presentation."

**Expected Behavior:**
1. Skill identifies: Slide count mismatched (45 vs ~15)
2. Flags as Major Issue: "Slide count drastically mismatched to duration"
3. Recommends: Reduce to 12-18 slides
4. Explains: 1 slide per minute guideline

**Success Criteria:**
- [ ] Mismatch identified (45 slides too many)
- [ ] 1 slide/minute guideline cited
- [ ] Specific reduction recommended (to ~15)

---

## venue-templates

### Scenario 1: Basic Usage

**Input:** "I'm submitting to Nature. What are the formatting requirements?"

**Expected Behavior:**
1. Skill provides Nature-specific guidelines
2. Includes: Word limits, figure requirements, reference style
3. References Nature/Science formatting details
4. Covers: Abstract length, Methods section format, supplementary info
5. Cites reviewer_expectations.md for what reviewers look for

**Success Criteria:**
- [ ] Nature-specific requirements provided
- [ ] Word limits stated
- [ ] Reference style specified
- [ ] Reviewer expectations referenced

---

### Scenario 2: Edge Case

**Input:** "I'm writing a bioinformatics methods paper. Which journal should I submit to?"

**Expected Behavior:**
1. Skill identifies: Journal selection query
2. Provides bioinformatics journal options from bioinformatics_journals.md
3. Lists: Nucleic Acids Research, Bioinformatics (Oxford), Genome Biology, Nature Methods
4. Compares: Scope, requirements, data deposition policies
5. Offers decision guidance based on paper type

**Success Criteria:**
- [ ] Multiple relevant journals listed
- [ ] Comparison criteria provided (scope, requirements)
- [ ] bioinformatics_journals.md referenced
- [ ] Decision guidance offered

---

## visual-design

### Scenario 1: Basic Usage

**Input:** "What are the design principles for creating effective scientific figures?"

**Expected Behavior:**
1. Skill provides design philosophy guidance
2. Covers: Hierarchy, contrast, alignment, proximity
3. Includes: Color theory for scientific visualization
4. Addresses: Accessibility (colorblind-safe palettes)
5. References publication specifications (DPI, font sizes)

**Success Criteria:**
- [ ] Core design principles explained
- [ ] Accessibility requirements included
- [ ] Publication specifications provided
- [ ] Colorblind-safe guidance included

---

### Scenario 2: Edge Case

**Input:** "Create a matplotlib heatmap for my gene expression data."

**Expected Behavior:**
1. Skill identifies: Code generation request, not design principles
2. Routes to: `plotting-libraries` skill
3. Explains: visual-design provides philosophy/principles, plotting-libraries provides code
4. May offer design guidance to apply after routing

**Success Criteria:**
- [ ] Scope boundary recognized (principles vs. code)
- [ ] Routes to plotting-libraries for implementation
- [ ] Relationship between skills clarified

---
