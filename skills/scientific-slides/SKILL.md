---
name: scientific-slides
version: 2.0.0
description: "Build slide decks and presentations for research talks. Use for PowerPoint slides, conference presentations, seminar talks, research presentations, thesis defense slides, or any scientific talk."
allowed-tools: [Read, Write, Edit, Bash]
---

# Scientific Slides

<overview>
Create effective scientific presentations for conferences, seminars, defenses, and professional talks. This skill provides decision frameworks for selecting presentation format, structured workflows for development, and integration with AI slide generation tools. Supports PDF (via Nano Banana Pro), PowerPoint (via PPTX skill), and LaTeX Beamer outputs.

**Core Philosophy**: Visual-first design with research-backed context. Every slide needs strong visuals; text supports but doesn't replace them.
</overview>

<when_to_use>
## Trigger Conditions

Use this skill when:
- Preparing conference presentations (5-20 minutes)
- Developing academic seminars (45-60 minutes)
- Creating thesis or dissertation defense presentations
- Designing grant pitch presentations
- Preparing journal club presentations
- Giving research talks at institutions or companies

Do NOT use this skill when:
- Creating posters → use `latex-posters` or `pptx-posters`
- Writing papers → use `scientific-writing`
- Creating standalone figures → use `plotting-libraries` or `scientific-schematics`
- Need only data visualizations → use `plotting-libraries`
</when_to_use>

<decision_framework>
## Talk Type Selection

```
What type of presentation?
│
├─ Conference talk (5-20 min)
│   └─ Focus: 1-2 key findings, minimal methods, memorable impact
│       Slides: 8-20, ~1 min/slide
│       See: references/talk_types_guide.md
│
├─ Academic seminar (45-60 min)
│   └─ Focus: Comprehensive coverage, detailed methods, multiple studies
│       Slides: 40-60, allow for questions
│       See: references/talk_types_guide.md
│
├─ Thesis defense (45-60 min)
│   └─ Focus: Complete dissertation, demonstrate mastery
│       Slides: 45-60 + backup slides
│       See: references/talk_types_guide.md
│
├─ Grant pitch (10-20 min)
│   └─ Focus: Significance, feasibility, impact, preliminary data
│       Slides: 10-20, persuasive structure
│       See: references/talk_types_guide.md
│
└─ Journal club (20-45 min)
    └─ Focus: Critical analysis of published work
        Slides: 20-35, discussion-oriented
        See: references/talk_types_guide.md
```

## Implementation Method Selection

```
Choose output format
│
├─ Need stunning visuals, fast creation, non-editable OK?
│   └─ YES → Nano Banana Pro PDF (DEFAULT, RECOMMENDED)
│       - Generate each slide as complete image
│       - Combine into PDF
│       - Best visual quality
│       - See: Stage 3A below
│
├─ Need editable slides, company templates, animations?
│   └─ YES → PowerPoint via PPTX Skill
│       - Generate visuals with --visual-only
│       - Build PPTX with text separately
│       - See: document-skills/pptx/SKILL.md
│
└─ Need equations, version control, reproducibility?
    └─ YES → LaTeX Beamer
        - Use templates in assets/
        - Compile to PDF
        - See: references/beamer_guide.md
```

## Time Allocation Matrix

| Talk Duration | Intro | Methods | Results | Discussion | Conclusion |
|---------------|-------|---------|---------|------------|------------|
| 10 min        | 15%   | 10%     | 50%     | 20%        | 5%         |
| 15 min        | 15%   | 15%     | 45%     | 20%        | 5%         |
| 30 min        | 15%   | 20%     | 40%     | 20%        | 5%         |
| 45 min        | 10%   | 20%     | 45%     | 20%        | 5%         |
| 60 min        | 10%   | 20%     | 45%     | 20%        | 5%         |

</decision_framework>

<workflow>
## Workflow

### Stage 1: Research & Planning

**Objective:** Define scope, gather citations, create content outline

**Steps:**
1. Determine talk type using decision tree above
2. Use `research-lookup` skill to find 8-15 relevant papers for citations
3. Identify 1-3 core messages (max 3 for short talks)
4. Select key figures/results (3-6 for 15-min talk)
5. Create detailed slide-by-slide plan:
   ```markdown
   ## Slide Plan

   ### Slide 1: Title
   - Title text, speaker, affiliation
   - Visual: Modern abstract background

   ### Slide 2: Hook
   - Compelling problem statement
   - Citations: (Smith et al., 2023; Jones et al., 2024)
   - Visual: Relevant image or diagram
   ...
   ```

**Exit Criteria:**
- [ ] Talk type determined (conference/seminar/defense/grant/journal club)
- [ ] 8-15 papers identified for citations
- [ ] 1-3 core messages defined
- [ ] Slide-by-slide plan complete with visual elements noted
- [ ] Target slide count set (duration × 1 slide/min ± 20%)

---

### Stage 2: Implementation Method Selection

**Objective:** Choose and configure output format

**Steps:**
1. Use Implementation Method decision tree above
2. For Nano Banana Pro PDF (default):
   - Ensure OPENROUTER_API_KEY is set
   - Create slides/ output directory
3. For PowerPoint:
   - Read document-skills/pptx/SKILL.md
   - Prepare template or HTML source
4. For LaTeX Beamer:
   - Select template from assets/
   - Configure theme and colors

**Exit Criteria:**
- [ ] Implementation method selected
- [ ] Required tools/APIs available
- [ ] Output directory created
- [ ] Template or starting point ready

---

### Stage 3A: Nano Banana Pro PDF Generation (Default)

**Objective:** Generate visually stunning slides as images

**Steps:**
1. Define formatting goal for consistency:
   ```
   FORMATTING GOAL: [background color], [text color], [accent color],
   minimal professional design, sans-serif fonts, generous margins
   ```

2. Generate title slide (establishes style):
   ```bash
   python scripts/generate_slide_image.py "Title slide: '[Title]'.
   Subtitle: '[Conference]'. Speaker: [Name].
   FORMATTING GOAL: Dark blue (#1a237e), white text, gold accents (#ffc107),
   minimal design." -o slides/01_title.png
   ```

3. Generate subsequent slides (ALWAYS attach previous):
   ```bash
   python scripts/generate_slide_image.py "Slide titled '[Title]'.
   Key points: 1) [Point 1], 2) [Point 2].
   CITATIONS: Include at bottom: (Author et al., Year).
   FORMATTING GOAL: Match attached slide style exactly." \
   -o slides/02_intro.png --attach slides/01_title.png
   ```

4. For results slides, attach existing figures:
   ```bash
   # Check for existing figures first
   ls figures/ results/ plots/

   # Attach data figures to results slides
   python scripts/generate_slide_image.py "Slide titled 'Results'.
   Present attached chart with key findings.
   FORMATTING GOAL: Match attached slide style." \
   -o slides/05_results.png --attach slides/04_methods.png \
   --attach figures/accuracy_chart.png
   ```

5. Combine to PDF:
   ```bash
   python scripts/slides_to_pdf.py slides/*.png -o presentation.pdf
   ```

**Exit Criteria:**
- [ ] All slides generated with consistent formatting
- [ ] Each subsequent slide attached previous for style matching
- [ ] Citations included where appropriate
- [ ] Data figures attached to results slides
- [ ] PDF successfully generated

---

### Stage 3B: PowerPoint Generation (Alternative)

**Objective:** Create editable slides with AI-generated visuals

**Steps:**
1. Generate visuals with `--visual-only` flag:
   ```bash
   python scripts/generate_slide_image.py "[Visual description]" \
   -o figures/slide_visual.png --visual-only
   ```

2. Build PPTX using document-skills/pptx workflow

3. Add text content separately from visuals

**Exit Criteria:**
- [ ] All slide visuals generated
- [ ] PPTX created with proper layout
- [ ] Text and visuals integrated
- [ ] Animations/transitions added if needed

---

### Stage 3C: LaTeX Beamer Generation (Alternative)

**Objective:** Create reproducible slides with equations

**Steps:**
1. Copy template from assets/:
   ```bash
   cp assets/beamer_template_conference.tex presentation.tex
   ```

2. Customize theme and colors

3. Add content in LaTeX

4. Compile:
   ```bash
   pdflatex presentation.tex
   ```

**Exit Criteria:**
- [ ] LaTeX compiles without errors
- [ ] All equations render correctly
- [ ] Figures display properly
- [ ] Bibliography working (if applicable)

---

### Stage 4: Visual Validation

**Objective:** Systematic review and issue correction

**Steps:**
1. Convert to images for review (if not already):
   ```bash
   python scripts/pdf_to_images.py presentation.pdf review/slide --dpi 150
   ```

2. Check each slide against validation checklist:

   | Check | Pass Criteria |
   |-------|---------------|
   | Text overflow | No text extends beyond boundaries |
   | Element overlap | No overlapping text/figures |
   | Font size | Body ≥24pt, titles ≥36pt |
   | Contrast | Text-background ratio ≥4.5:1 |
   | Visual content | ≥1 strong visual per slide |
   | White space | 40-50% of slide empty |

3. Document issues:
   ```
   Slide # | Issue | Priority
   --------|-------|----------
   3       | Text overflow in bullet 4 | High
   7       | Figure overlaps caption | High
   ```

4. Fix issues and regenerate

**Exit Criteria:**
- [ ] All slides visually inspected
- [ ] No text overflow or overlap
- [ ] All fonts readable (≥24pt body)
- [ ] Adequate contrast throughout
- [ ] Visual content on every slide

---

### Stage 5: Practice & Refinement

**Objective:** Ensure timing and delivery quality

**Steps:**
1. Practice with timer (minimum 3 runs):
   - Run 1: Rough draft timing
   - Run 2: Smooth transitions
   - Run 3: Exact timing

2. Set timing checkpoints for 15-min talk:
   - 3-4 min: Finishing introduction
   - 7-8 min: Halfway through results
   - 12-13 min: Starting conclusions

3. Prepare backup slides for Q&A

4. Create multiple copies (laptop, cloud, USB)

**Exit Criteria:**
- [ ] Practiced ≥3 times with timer
- [ ] Timing within ±10% of target
- [ ] Timing checkpoints identified
- [ ] Backup slides prepared
- [ ] Multiple file copies saved

</workflow>

<success_criteria>
## Success Criteria

**Quantitative Thresholds:**

| Metric | Minimum | Target | Excellent |
|--------|---------|--------|-----------|
| Slides per minute | 0.7 | 1.0 | 1.2 |
| Body font size | 18pt | 24pt | 28pt |
| Title font size | 32pt | 36pt | 44pt |
| Visual content | 40% | 60% | 70% |
| White space | 30% | 45% | 50% |
| Practice runs | 2 | 3 | 5 |
| Contrast ratio | 4.5:1 | 7:1 | 10:1 |
| Citations (intro) | 2 | 3-5 | 6+ |
| Citations (discussion) | 2 | 3-5 | 6+ |

**Completion Checklist:**
- [ ] Slide count appropriate for duration (±20%)
- [ ] Every slide has strong visual element
- [ ] No text-only bullet point slides
- [ ] Font sizes ≥24pt body, ≥36pt titles
- [ ] High contrast (≥4.5:1) throughout
- [ ] Research context established with citations
- [ ] Clear narrative arc (hook → problem → approach → results → implications)
- [ ] Practiced ≥3 times with exact timing
- [ ] Backup slides prepared
- [ ] Multiple file copies saved

</success_criteria>

<scope>
## Scope

**In Scope:**
- Conference presentations (5-60 minutes)
- Academic seminars and colloquia
- Thesis and dissertation defenses
- Grant pitch presentations
- Journal club presentations
- Research talks
- PDF, PowerPoint, and LaTeX Beamer outputs
- AI-generated slide visuals
- Visual validation workflow

**Out of Scope** (use specialized resources):
- Research posters → use `latex-posters` or `pptx-posters`
- Scientific papers → use `scientific-writing`
- Data visualization code → use `plotting-libraries`
- Technical diagrams → use `scientific-schematics`
- Photorealistic images → use `generate-image`
- Document conversion → use `markitdown`

</scope>

<anti_patterns>
## Common Pitfalls

### 1. Text-Heavy Slides

**Anti-pattern:**
- Walls of text (>6 bullets per slide)
- Full sentences instead of keywords
- Reading slides verbatim
- No visual elements

**Solution:**
- Maximum 3-4 bullets, 4-6 words each
- Use visuals as primary content (60-70%)
- Text supports visuals, doesn't replace them
- If you need more text, split into multiple slides

---

### 2. Missing Research Context

**Anti-pattern:**
- No citations in introduction or discussion
- Claims without supporting evidence
- Unclear positioning relative to prior work

**Solution:**
- Use `research-lookup` to find 8-15 papers before creating slides
- Cite 3-5 papers in introduction establishing context
- Cite 3-5 papers in discussion for comparison
- Include citations directly in slide generation prompts

---

### 3. Inconsistent Formatting

**Anti-pattern:**
- Different colors/fonts across slides
- Varying layout styles
- No visual continuity

**Solution:**
- Define FORMATTING GOAL in first slide prompt
- ALWAYS attach previous slide when generating subsequent slides
- Use consistent color scheme (3-5 colors max)
- Maintain same font family throughout

---

### 4. Skipping Visual Validation

**Anti-pattern:**
- Not reviewing generated slides
- Missing text overflow or overlap issues
- Unreadable fonts in final presentation

**Solution:**
- Convert to images for systematic review
- Check every slide against validation checklist
- Fix issues and regenerate
- Test on presentation computer before event

---

### 5. Insufficient Practice

**Anti-pattern:**
- First run-through is during actual presentation
- No timing checkpoints
- No backup plan for running over

**Solution:**
- Practice minimum 3 times with timer
- Set 3-4 timing checkpoints
- Prepare skip-able slides for time pressure
- Never skip conclusions—cut earlier content instead

</anti_patterns>

<templates>
## Output Templates

### Slide Plan Template

```markdown
# Presentation: [TITLE]
Duration: [X] minutes | Slides: [N] | Type: [conference/seminar/defense/grant/journal]

## Formatting Goal
Background: [color]
Text: [color]
Accent: [color]
Style: [minimal/corporate/academic]

## Slide 1: Title
- Title: "[Full title]"
- Subtitle: "[Conference/Event]"
- Speaker: [Name, Affiliation]
- Visual: [Background description]

## Slide 2: Hook
- Title: "[Compelling question or statement]"
- Key point: [What grabs attention]
- Citations: (Author1 et al., Year; Author2 et al., Year)
- Visual: [Image/diagram description]

## Slide 3: Background
- Title: "[Context title]"
- Points: [2-3 bullet points]
- Citations: [Key papers]
- Visual: [Diagram/image]

... [continue for all slides]

## Backup Slides
- [Topic 1]: For Q&A about [specific topic]
- [Topic 2]: Additional data if asked
```

### Timing Checkpoint Template

```markdown
# Timing Checkpoints: [TITLE]
Total duration: [X] minutes

| Checkpoint | Target Time | Actual Run 1 | Actual Run 2 | Actual Run 3 |
|------------|-------------|--------------|--------------|--------------|
| End of intro | [X:XX] | | | |
| Start of results | [X:XX] | | | |
| Halfway through results | [X:XX] | | | |
| Start of conclusions | [X:XX] | | | |
| End | [X:XX] | | | |

Notes:
- Run 1: [observations]
- Run 2: [observations]
- Run 3: [observations]
```

</templates>

<cross_references>
## Related Skills

| Skill | Relationship |
|-------|-------------|
| `research-lookup` | Use BEFORE creating slides to gather 8-15 citations |
| `scientific-schematics` | Use for complex technical diagrams |
| `plotting-libraries` | Use for data visualizations to include in slides |
| `latex-posters` | Use for research posters (not presentations) |
| `pptx-posters` | Use for PowerPoint posters (not presentations) |
| `document-skills/pptx` | Use for editable PowerPoint creation |
| `venue-templates` | Consult for conference-specific styles |

**Skill Selection:**
See `SKILL_ROUTER.md` for decision trees when multiple skills may apply.

</cross_references>

<references>
## Reference Documents

| Document | Purpose |
|----------|---------|
| `references/presentation_structure.md` | Detailed structure for all talk types, timing allocation |
| `references/slide_design_principles.md` | Typography, color theory, layout, accessibility |
| `references/data_visualization_slides.md` | Simplifying figures, chart types, progressive disclosure |
| `references/talk_types_guide.md` | Specific guidance for conferences, seminars, defenses, grants |
| `references/beamer_guide.md` | Complete LaTeX Beamer documentation |
| `references/visual_review_workflow.md` | PDF to images, systematic inspection, iteration |
| `assets/powerpoint_design_guide.md` | PowerPoint design and implementation |
| `assets/timing_guidelines.md` | Timing, pacing, and practice strategies |

## Script Reference

| Script | Usage |
|--------|-------|
| `scripts/generate_slide_image.py` | Generate slides with Nano Banana Pro |
| `scripts/slides_to_pdf.py` | Combine slide images into PDF |
| `scripts/pdf_to_images.py` | Convert PDF to images for review |
| `scripts/validate_presentation.py` | Validate slide count, timing, file size |

</references>
