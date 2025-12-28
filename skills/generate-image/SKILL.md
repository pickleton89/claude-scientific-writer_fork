---
name: generate-image
version: 2.0.0
description: "Generate and edit images using AI models (FLUX, Gemini) via OpenRouter. For photos, illustrations, artwork, and visual assets. Use scientific-schematics for technical diagrams."
allowed-tools: [Read, Write, Edit, Bash]
---

# Generate Image

<overview>
Generate and edit high-quality images using OpenRouter's image generation models including FLUX.2 Pro and Gemini 3 Pro. Produces photos, illustrations, artwork, and visual assets for scientific documents, presentations, and publications. For technical diagrams (flowcharts, circuits, pathways), use `scientific-schematics` instead.
</overview>

<when_to_use>
## Trigger Conditions

Use this skill when:
- Generating photorealistic images for documents
- Creating artistic illustrations and concept art
- Producing visual assets for presentations or posters
- Editing existing images (background changes, modifications)
- Creating conceptual visualizations for scientific concepts

Do NOT use this skill when:
- Creating flowcharts or process diagrams → use `scientific-schematics`
- Creating circuit diagrams or electrical schematics → use `scientific-schematics`
- Creating biological pathways or signaling cascades → use `scientific-schematics`
- Creating system architecture diagrams → use `scientific-schematics`
- Creating data visualizations (charts, plots) → use `plotting-libraries`
</when_to_use>

<decision_framework>
## Image Type Selection

```
What kind of image is needed?
│
├─ Technical diagram?
│  │
│  ├─ Flowchart, process diagram → scientific-schematics
│  ├─ Circuit, electrical schematic → scientific-schematics
│  ├─ Biological pathway → scientific-schematics
│  └─ Architecture diagram → scientific-schematics
│
├─ Data visualization?
│  │
│  └─ Charts, plots, graphs → plotting-libraries
│
└─ Photo, illustration, or artwork?
   │
   ├─ Photorealistic image → generate-image (this skill)
   ├─ Scientific illustration → generate-image (this skill)
   ├─ Concept art → generate-image (this skill)
   └─ Image editing → generate-image (this skill)
```

## Model Selection Matrix

| Use Case | Recommended Model | Rationale |
|----------|------------------|-----------|
| High-quality generation | `google/gemini-3-pro-image-preview` | Best quality |
| Fast generation | `black-forest-labs/flux.2-pro` | Speed + quality |
| Image editing | `google/gemini-3-pro-image-preview` | Best editing |
| Budget generation | `black-forest-labs/flux.2-flex` | Lower cost |

## Model Capabilities

| Model | Generation | Editing | Quality | Cost |
|-------|------------|---------|---------|------|
| gemini-3-pro-image-preview | Yes | Yes | Excellent | $$ |
| flux.2-pro | Yes | Yes | Excellent | $$ |
| flux.2-flex | Yes | No | Good | $ |

</decision_framework>

<workflow>
## Image Generation Workflow

### Stage 1: Requirement Analysis

**Objective:** Define image requirements and select approach

**Steps:**
1. Clarify the image purpose (document, presentation, poster)
2. Determine image type using decision tree
3. Select appropriate model based on requirements
4. Define output specifications (size, format, style)

**Exit Criteria:**
- [ ] Image purpose defined
- [ ] Confirmed this skill is appropriate (not scientific-schematics)
- [ ] Model selected from matrix
- [ ] Output format specified (PNG, size)

---

### Stage 2: Prompt Engineering

**Objective:** Craft effective prompt for optimal generation

**Steps:**
1. Start with subject description
2. Add style and medium specifications
3. Include technical requirements
4. Add quality modifiers

**Prompt Structure:**

```
[SUBJECT] + [STYLE] + [MEDIUM] + [TECHNICAL] + [QUALITY]
```

**Component Guide:**

| Component | Examples |
|-----------|----------|
| Subject | "DNA double helix", "laboratory setting", "cancer cells" |
| Style | "scientific illustration", "photorealistic", "modern", "minimalist" |
| Medium | "digital art", "3D render", "photograph", "watercolor" |
| Technical | "high contrast", "neutral background", "well-lit" |
| Quality | "detailed", "professional", "publication-quality", "4K" |

**Scientific Image Prompts:**

| Image Type | Prompt Pattern |
|------------|----------------|
| Conceptual | "Conceptual illustration of [topic], scientific visualization style, clean background" |
| Microscopy | "Microscopic view of [subject], scientific imaging style, detailed, high contrast" |
| Lab setting | "Modern laboratory with [equipment], photorealistic, professional lighting" |
| Molecular | "[Molecule/structure] visualization, 3D render, scientific accuracy, clean background" |

**Exit Criteria:**
- [ ] Prompt includes all 5 components
- [ ] Prompt is specific (no vague terms)
- [ ] Style appropriate for intended use
- [ ] Technical requirements specified

---

### Stage 3: Generation Execution

**Objective:** Generate image using script

**Steps:**
1. Verify API key configuration
2. Execute generation command
3. Monitor for errors
4. Save output

**Commands:**

**Basic Generation:**
```bash
python scripts/generate_image.py "Your detailed prompt here"
```

**With Model Selection:**
```bash
python scripts/generate_image.py "Your prompt" --model "google/gemini-3-pro-image-preview"
```

**Custom Output Path:**
```bash
python scripts/generate_image.py "Your prompt" --output figures/my_image.png
```

**Image Editing:**
```bash
python scripts/generate_image.py "Make the background blue" --input photo.jpg
```

**API Key Configuration:**
```bash
# Check .env file exists
cat .env | grep OPENROUTER_API_KEY

# Or set environment variable
export OPENROUTER_API_KEY=your-api-key-here
```

**Exit Criteria:**
- [ ] Command executed without errors
- [ ] Output file created
- [ ] File size reasonable (not corrupted)
- [ ] Image opens correctly

---

### Stage 4: Quality Validation

**Objective:** Verify image meets requirements

**Steps:**
1. Open and inspect generated image
2. Check against quality criteria
3. Assess fitness for intended use
4. Decide: accept, regenerate, or edit

**Quality Checklist:**

| Criterion | Pass | Fail Action |
|-----------|------|-------------|
| Subject correct | Matches prompt | Refine prompt |
| Style appropriate | Fits intended use | Add style modifiers |
| Resolution sufficient | ≥1024px dimension | Specify size |
| No artifacts | Clean, professional | Regenerate |
| Color appropriate | Matches context | Specify colors |
| Composition balanced | Visually appealing | Add composition terms |

**Common Issues and Fixes:**

| Issue | Prompt Fix |
|-------|------------|
| Wrong style | Add explicit style: "in the style of scientific illustration" |
| Low detail | Add: "highly detailed", "intricate" |
| Wrong colors | Add: "blue tones", "monochromatic", specific colors |
| Busy composition | Add: "simple composition", "clean background", "minimalist" |
| Generic look | Add specific details about the subject |

**Exit Criteria:**
- [ ] Image passes quality checklist
- [ ] Appropriate for intended use
- [ ] Resolution sufficient for output format
- [ ] No visible artifacts or errors

---

### Stage 5: Integration

**Objective:** Incorporate image into project

**Steps:**
1. Move to appropriate directory
2. Rename following project conventions
3. Add to document/presentation
4. Update references if needed

**File Organization:**
```
project/
├── figures/
│   ├── fig1_concept_illustration.png
│   ├── fig2_lab_setting.png
│   └── fig3_molecular_view.png
└── slides/
    └── background_hero.png
```

**Exit Criteria:**
- [ ] File in correct project location
- [ ] Named following conventions
- [ ] Referenced in document/presentation
- [ ] Backup preserved

</workflow>

<success_criteria>
## Success Criteria

**Quantitative Thresholds:**

| Metric | Minimum | Target | Excellent |
|--------|---------|--------|-----------|
| Resolution | 512px | 1024px | 2048px |
| Generation time | <60s | <30s | <10s |
| Iteration count | ≤5 | ≤3 | 1 |
| Prompt length | 10 words | 20-50 words | Optimal coverage |

**Completion Checklist:**
- [ ] Image matches intended subject
- [ ] Style appropriate for use case
- [ ] Resolution sufficient for output
- [ ] No artifacts or distortions
- [ ] Colors and composition professional
- [ ] File saved in correct format and location

</success_criteria>

<scope>
## Scope

**In Scope:**
- Photorealistic image generation
- Artistic illustrations and concept art
- Scientific conceptual visualizations
- Image editing and modification
- Visual assets for documents and presentations
- Background images and hero shots

**Out of Scope** (use specialized resources):
- Flowcharts, process diagrams → use `scientific-schematics`
- Circuit/electrical diagrams → use `scientific-schematics`
- Biological pathways → use `scientific-schematics`
- Data charts and plots → use `plotting-libraries`
- CONSORT/methodology diagrams → use `scientific-schematics`
- Technical architecture diagrams → use `scientific-schematics`

</scope>

<anti_patterns>
## Common Pitfalls

### 1. Vague Prompts

**Anti-pattern:**
```bash
python scripts/generate_image.py "a cell"
# Result: Generic, unusable image
```

**Solution:**
```bash
python scripts/generate_image.py "Microscopic view of a human cancer cell being attacked by immunotherapy T-cells, scientific illustration style, high detail, colorful, professional medical visualization"
```

---

### 2. Missing API Key

**Anti-pattern:**
```bash
python scripts/generate_image.py "prompt"
# Error: API key not found
```

**Solution:**
```bash
# Check for API key first
if [ -z "$OPENROUTER_API_KEY" ] && ! grep -q OPENROUTER_API_KEY .env 2>/dev/null; then
    echo "API key required. Add to .env: OPENROUTER_API_KEY=your-key"
    exit 1
fi

python scripts/generate_image.py "prompt"
```

---

### 3. Wrong Skill Selection

**Anti-pattern:**
Using generate-image for flowcharts or diagrams (produces artistic interpretation, not accurate diagram)

**Solution:**
Use decision tree:
- Technical diagrams → `scientific-schematics`
- Data plots → `plotting-libraries`
- Photos/illustrations → `generate-image`

---

### 4. No Quality Validation

**Anti-pattern:**
```bash
python scripts/generate_image.py "prompt" --output final_figure.png
# Used directly without checking
```

**Solution:**
```bash
# Generate to temp location
python scripts/generate_image.py "prompt" --output temp_image.png

# Open and validate
open temp_image.png  # macOS
# Check: subject, style, resolution, artifacts

# If acceptable, move to final location
mv temp_image.png figures/fig1_final.png
```

---

### 5. Excessive Iterations

**Anti-pattern:**
Regenerating 10+ times with minor prompt tweaks

**Solution:**
Systematic refinement:
1. First generation: Test basic prompt
2. If wrong subject: Clarify subject terms
3. If wrong style: Add explicit style modifiers
4. If wrong composition: Add composition terms
5. Maximum 5 iterations; if still failing, reconsider approach

</anti_patterns>

<templates>
## Prompt Templates

### Template 1: Scientific Illustration

```
[SUBJECT DESCRIPTION], scientific illustration style,
clean white background, high detail,
professional medical/scientific visualization,
suitable for academic publication
```

**Example:**
```bash
python scripts/generate_image.py "DNA double helix structure with highlighted mutation site, scientific illustration style, clean white background, high detail, professional molecular visualization, suitable for academic publication" --output figures/dna_mutation.png
```

### Template 2: Conceptual Visualization

```
Conceptual illustration representing [CONCEPT],
modern scientific visualization,
abstract yet professional,
[COLOR SCHEME] color palette,
clean composition
```

**Example:**
```bash
python scripts/generate_image.py "Conceptual illustration representing machine learning in drug discovery, modern scientific visualization, abstract yet professional, blue and teal color palette, clean composition" --output figures/ml_drug_discovery.png
```

### Template 3: Laboratory/Equipment

```
Modern [EQUIPMENT/SETTING] in a contemporary research laboratory,
photorealistic, professional photography style,
well-lit, clean environment,
suitable for presentation or publication
```

**Example:**
```bash
python scripts/generate_image.py "Modern CRISPR gene editing equipment in a contemporary research laboratory, photorealistic, professional photography style, well-lit, clean environment, suitable for presentation" --output slides/lab_hero.png
```

### Template 4: Image Editing

```
[EDIT INSTRUCTION]. Keep the [ELEMENTS TO PRESERVE].
Make the change look natural and professional.
```

**Example:**
```bash
python scripts/generate_image.py "Change the background to a clean white color. Keep the main subject and its details. Make the change look natural and professional." --input original.png --output edited.png
```

</templates>

<cross_references>
## Related Skills

| Skill | Relationship |
|-------|-------------|
| `scientific-schematics` | Use for technical diagrams (flowcharts, pathways, circuits) |
| `plotting-libraries` | Use for data visualizations (charts, plots, graphs) |
| `scientific-slides` | Combine for visually rich presentations |
| `latex-posters` | Combine for poster visuals and hero images |
| `visual-design` | Consult for design principles and accessibility |

**Skill Selection:**
See `SKILL_ROUTER.md` for decision trees when multiple skills may apply.

</cross_references>

<references>
## Reference Documents

| Document | Purpose |
|----------|---------|
| `scripts/generate_image.py` | Main generation script |

## External Resources

- **OpenRouter**: https://openrouter.ai
- **API Keys**: https://openrouter.ai/keys
- **Model List**: https://openrouter.ai/models

</references>
