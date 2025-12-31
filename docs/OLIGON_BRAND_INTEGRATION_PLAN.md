# Oligon-Brand Skill Integration Plan

> Implementation plan for integrating the `oligon-brand` skill into the claude-scientific-writer skill library
> Created: 2025-12-30
> Status: **In Progress** (Phase 0 Complete)

---

## Executive Summary

The `oligon-brand` skill was migrated from `~/.claude/skills/` and needs restructuring to match the skill library's conventions. This plan transforms it from a standalone brand guidelines document into a properly integrated Claude skill that:

1. Follows the standard skill XML structure
2. Integrates with SKILL_ROUTER.md for proper routing
3. Cross-references visual-design and plotting-libraries
4. Supports future alternative brand implementations

**Architecture Decision:** Independent skill with bidirectional references (Option C from analysis).

---

## Current State

### Files in `skills/oligon-brand/`

```
skills/oligon-brand/
├── SKILL.md                          # ⚠️ Brand guidelines, not skill format
├── adapters/
│   ├── matplotlib_adapter.py         # ✅ Complete adapter
│   ├── ggplot2_adapter.R             # ✅ Complete adapter
│   ├── reportlab_adapter.py          # ✅ PDF generation
│   ├── docx_adapter.py               # ✅ Word document styling
│   ├── pptx_adapter.py               # ✅ PowerPoint styling
│   └── html_adapter.css              # ✅ CSS variables
├── tokens/
│   └── brand-tokens.json             # ✅ Single source of truth
├── references/
│   └── brand-colors-full.md          # ✅ V4 complete specification
├── scripts/
│   ├── matplotlib_brand_setup.py     # ✅ Setup automation
│   └── ggplot2_brand_setup.R         # ✅ R setup automation
└── assets/
    ├── oligon_color_brand.mplstyle   # ✅ Matplotlib style file
    ├── Oligon_Logo_BrandBlue.svg     # ✅ Logo assets
    ├── Oligon_Logo_White.svg
    ├── Oligon_Icon_BrandBlue.svg
    └── Oligon_Icon_White.svg
```

### Issues to Fix

| Issue | Current | Target |
|-------|---------|--------|
| SKILL.md structure | Markdown guidelines | XML skill sections |
| Frontmatter | `name`, `description` only | Full metadata |
| SKILL_ROUTER.md | Not included | Brand routing section |
| visual-design cross-ref | Mentions but undefined | Clear relationship |
| plotting-libraries | Generic brand import | Concrete oligon-brand example |

---

## Implementation Phases

### Phase 0: Preparation ✅ COMPLETE
- [x] Read and backup current SKILL.md content (→ `SKILL.md.backup`)
- [x] Verify all adapter files are functional (6 adapters verified)
- [x] Create implementation tracking in this document

**Phase 0 Verification Results:**
- `brand-tokens.json`: Valid JSON, v4.0, all color categories present
- `matplotlib_adapter.py`: Well-structured, requires matplotlib dependency
- `ggplot2_adapter.R`: Well-structured, requires ggplot2 dependency
- `reportlab_adapter.py`: Well-structured, requires reportlab dependency
- `docx_adapter.py`: Well-structured, requires python-docx dependency
- `pptx_adapter.py`: Well-structured, requires python-pptx dependency
- `html_adapter.css`: Complete, no dependencies

### Phase 1: SKILL.md Restructure
- [ ] Step 1.1: Create new SKILL.md with proper frontmatter
- [ ] Step 1.2: Add all required XML sections
- [ ] Step 1.3: Move detailed guidelines to quick-reference format
- [ ] Step 1.4: Verify references to brand-colors-full.md

### Phase 2: SKILL_ROUTER.md Integration
- [ ] Step 2.1: Add "Brand Application" decision tree
- [ ] Step 2.2: Update skill count
- [ ] Step 2.3: Add oligon-brand to Quick Reference Matrix

### Phase 3: Cross-Reference Updates
- [ ] Step 3.1: Update visual-design/SKILL.md cross-references
- [ ] Step 3.2: Update plotting-libraries/SKILL.md with concrete examples
- [ ] Step 3.3: Update any other skills that reference brand colors

### Phase 4: Future-Proofing
- [ ] Step 4.1: Document brand skill interface contract
- [ ] Step 4.2: Add brand selection to SKILL_ROUTER.md

### Phase 5: Validation & Documentation
- [ ] Step 5.1: Test adapter imports
- [ ] Step 5.2: Verify cross-reference paths
- [ ] Step 5.3: Update CHANGELOG.md

---

## Detailed Implementation Steps

### Phase 1: SKILL.md Restructure

#### Step 1.1: New Frontmatter

Replace current frontmatter with:

```yaml
---
name: oligon-brand
version: 4.0.0
description: "Oligon scientific brand implementation for figures, presentations, and documents. Provides color palettes, typography standards, and format-specific adapters for matplotlib, ggplot2, ReportLab, python-docx, python-pptx, and HTML/CSS."
allowed-tools: [Read, Glob]
brand-type: organizational
tokens-file: tokens/brand-tokens.json
---
```

#### Step 1.2: Required XML Sections

The new SKILL.md must include these sections (in order):

```markdown
# Oligon Scientific Brand

<overview>
[Brief description of purpose and core principle]
</overview>

<when_to_use>
## Trigger Conditions

Use this skill when:
- Creating figures for Oligon-branded publications
- Applying consistent brand colors to matplotlib/seaborn/ggplot2 plots
- Styling presentations or documents with Oligon identity
- Need to access brand color hex codes programmatically

Do NOT use this skill when:
- Working on non-Oligon projects → use visual-design defaults
- Need general design philosophy → use visual-design
- Writing plotting code → use plotting-libraries (which references this skill)
</when_to_use>

<prerequisites>
## Prerequisites

**For Python (matplotlib/seaborn):**
```bash
pip install matplotlib numpy
```

**For R (ggplot2):**
```r
install.packages("ggplot2")
```

**Setup (one-time):**
```bash
# Install mplstyle to matplotlib config
python skills/oligon-brand/scripts/matplotlib_brand_setup.py

# Or for R
Rscript skills/oligon-brand/scripts/ggplot2_brand_setup.R
```
</prerequisites>

<decision_framework>
## Color Cycle Selection

[Move existing decision tree here - the "How many data series?" flowchart]

### Quick Selection Table

| Context | Cycle | Colors |
|---------|-------|--------|
| Control vs Treatment | `treatment_control` | Dark Gray → Brand Blue → Medium Gray |
| All equal weight | `neutrals_only` | Dark Gray → Medium Gray → Muted Gray |
| Up vs Down regulation | `opposing` | Brand Blue → Contrast Orange → Dark Gray |
| 4 categories (+ shapes!) | `multi_category` | Dark Gray → Brand Blue → Medium Gray → Medium Blue |
</decision_framework>

<workflow>
## Workflow

### Stage 1: Load Brand Tokens

**Python:**
```python
from skills.oligon_brand.adapters.matplotlib_adapter import (
    set_brand_style,
    BRAND_COLORS,
    get_cycle
)
```

**Or load JSON directly:**
```python
import json
from pathlib import Path

tokens_path = Path("skills/oligon-brand/tokens/brand-tokens.json")
with open(tokens_path) as f:
    brand = json.load(f)
```

### Stage 2: Apply Style

```python
# Set matplotlib defaults
set_brand_style('treatment_control')

# Or apply mplstyle
import matplotlib.pyplot as plt
plt.style.use('oligon_color_brand')
```

### Stage 3: Create Figure

```python
# Your plotting code here
# Brand colors are now the default cycle
```

### Stage 4: Preflight Check

Before finalizing, verify:
- [ ] White background
- [ ] Black axes and text
- [ ] Top/right spines removed
- [ ] Brand Blue limited to 1-2 elements
- [ ] Tested in grayscale
</workflow>

<success_criteria>
## Success Criteria

**Visual Compliance Checklist:**
- [ ] White background (`#FFFFFF`)
- [ ] Black text and axes (`#000000`)
- [ ] Top and right spines removed
- [ ] Correct color cycle for scientific context
- [ ] Brand Blue (`#2DB2E8`) limited to 1-2 elements per panel
- [ ] No library default colors
- [ ] Redundant encoding (shapes/styles with color)
- [ ] Tested in grayscale
- [ ] Error bars labeled (SEM/SD/95% CI, n=?)
- [ ] Arial font, 7-10 pt range
- [ ] Resolution ≥300 DPI for print
- [ ] Logo excluded from publication figures
</success_criteria>

<scope>
## Scope

**In Scope:**
- Color palette definitions and usage guidance
- Format-specific adapters (matplotlib, ggplot2, ReportLab, etc.)
- Typography standards for figures and documents
- Logo usage guidelines
- Preflight verification checklist

**Out of Scope:**
- Design philosophy and principles → use `visual-design`
- Writing plotting code → use `plotting-libraries`
- Journal-specific requirements → use `venue-templates`
- Diagram creation → use `scientific-schematics`
</scope>

<anti_patterns>
## Common Pitfalls

### 1. Using Library Defaults
**Problem:** Matplotlib/seaborn default colors instead of brand palette.
**Solution:** Always call `set_brand_style()` or load mplstyle before plotting.

### 2. Overusing Brand Blue
**Problem:** Too many elements highlighted, diluting emphasis.
**Solution:** Limit Brand Blue to 1-2 focal elements per panel. Use neutrals for everything else.

### 3. Color-Only Encoding
**Problem:** Relying solely on color to distinguish groups.
**Solution:** Combine color with marker shapes, line styles, or direct labels.

### 4. Ignoring Grayscale Test
**Problem:** Figure unreadable when printed in black and white.
**Solution:** Always test in grayscale before finalizing.

### 5. Brand Blue vs Medium Blue Confusion
**Problem:** Using `#2DB2E8` and `#158BBB` as primary differentiators (too similar).
**Solution:** Pair Brand Blue with neutrals or Contrast Orange, not Medium Blue.
</anti_patterns>

<error_handling>
## Error Handling

### Font Not Found
**Symptom:** Matplotlib warning about Arial font.
**Solution:**
```python
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
```

### Token File Not Found
**Symptom:** `FileNotFoundError` when loading brand-tokens.json
**Solution:** The matplotlib_adapter has embedded fallback values. For direct JSON loading, verify path:
```python
from pathlib import Path
tokens_path = Path(__file__).parent / "tokens" / "brand-tokens.json"
assert tokens_path.exists(), f"Tokens not found at {tokens_path}"
```

### Color Cycle Not Applying
**Symptom:** Default matplotlib colors appearing despite `set_brand_style()`.
**Solution:** Call `set_brand_style()` BEFORE creating any figures:
```python
set_brand_style('treatment_control')  # Must be first!
fig, ax = plt.subplots()
```
</error_handling>

<adapters>
## Format Adapters

| Format | Adapter | Usage |
|--------|---------|-------|
| matplotlib/seaborn | `adapters/matplotlib_adapter.py` | `from adapters.matplotlib_adapter import set_brand_style` |
| ggplot2 | `adapters/ggplot2_adapter.R` | `source("adapters/ggplot2_adapter.R")` |
| ReportLab (PDF) | `adapters/reportlab_adapter.py` | `from adapters.reportlab_adapter import get_brand_colors` |
| python-docx | `adapters/docx_adapter.py` | `from adapters.docx_adapter import apply_brand_styles` |
| python-pptx | `adapters/pptx_adapter.py` | `from adapters.pptx_adapter import get_brand_theme` |
| HTML/CSS | `adapters/html_adapter.css` | `<link rel="stylesheet" href="html_adapter.css">` |

See individual adapter files for detailed API documentation.
</adapters>

<cross_references>
## Related Skills

| Skill | Relationship |
|-------|-------------|
| `visual-design` | **Design philosophy** — provides universal principles; oligon-brand provides Oligon-specific implementation |
| `plotting-libraries` | **Code implementation** — use plotting-libraries for matplotlib/ggplot2 code patterns; load oligon-brand first for brand colors |
| `scientific-slides` | **Presentations** — apply oligon-brand colors to slide designs |
| `latex-posters` / `pptx-posters` | **Posters** — apply oligon-brand colors and logo placement |
| `markdown-to-pdf` | **PDF reports** — uses brand colors for styled output |
| `document-skills/docx` | **Word documents** — use docx_adapter for branded styling |

**Skill Selection:** See `SKILL_ROUTER.md` for decision trees when multiple skills may apply.

**Typical Workflow:**
1. `oligon-brand` → load tokens and set style
2. `plotting-libraries` → implement figure code
3. `visual-design` → review for design quality
</cross_references>

<references>
## Reference Documents

| Document | Purpose |
|----------|---------|
| `references/brand-colors-full.md` | Complete V4 brand specification with all color codes, typography, and embedded logo SVGs |
| `tokens/brand-tokens.json` | Machine-readable brand tokens for programmatic access |
| `assets/oligon_color_brand.mplstyle` | Matplotlib style file for one-line brand application |

### Color Quick Reference

| Role | Color | HEX |
|------|-------|-----|
| Brand Blue (hero) | ![#2DB2E8](https://via.placeholder.com/12/2DB2E8/2DB2E8.png) | `#2DB2E8` |
| Contrast Orange | ![#E8622D](https://via.placeholder.com/12/E8622D/E8622D.png) | `#E8622D` |
| Dark Gray (control) | ![#222222](https://via.placeholder.com/12/222222/222222.png) | `#222222` |
| Medium Gray | ![#666666](https://via.placeholder.com/12/666666/666666.png) | `#666666` |
| Muted Gray | ![#999999](https://via.placeholder.com/12/999999/999999.png) | `#999999` |
| Light Gray (annotations) | ![#BDBDBD](https://via.placeholder.com/12/BDBDBD/BDBDBD.png) | `#BDBDBD` |

For complete color tables, see `references/brand-colors-full.md`.
</references>
```

#### Step 1.3: Content Migration

Move detailed content as follows:

| Current Location | Target Location |
|------------------|-----------------|
| Color palette tables | `<references>` quick reference + `references/brand-colors-full.md` |
| Decision tree | `<decision_framework>` |
| Figure standards | `<success_criteria>` checklist |
| Typography tables | `references/brand-colors-full.md` (already there) |
| Logo usage | `references/brand-colors-full.md` (already there) |
| Preflight checklist | `<success_criteria>` |

#### Step 1.4: Verify Reference Paths

Ensure these paths are correct:
- `references/brand-colors-full.md` ✓ (exists)
- `tokens/brand-tokens.json` ✓ (exists)
- `assets/oligon_color_brand.mplstyle` ✓ (exists)

---

### Phase 2: SKILL_ROUTER.md Integration

#### Step 2.1: Add Brand Application Section

Insert new section after "Decision Tree: Figure & Visual Creation":

```markdown
---

## Decision Tree: Brand Application

```
User needs to apply organizational brand standards
│
├─ Which brand?
│  │
│  ├─ Oligon brand?
│  │  └─ YES → oligon-brand
│  │            │
│  │            ├─ Load brand tokens
│  │            ├─ Select appropriate color cycle
│  │            └─ Then proceed to format-specific skill:
│  │               ├─ Figures → plotting-libraries
│  │               ├─ Slides → scientific-slides
│  │               ├─ Posters → latex-posters / pptx-posters
│  │               ├─ Documents → document-skills/docx
│  │               └─ PDFs → markdown-to-pdf
│  │
│  └─ Other/None?
│     └─ Use visual-design defaults
│
├─ What output format?
│  │
│  ├─ matplotlib/seaborn figure?
│  │  └─ oligon-brand (matplotlib_adapter) → plotting-libraries
│  │
│  ├─ ggplot2 figure?
│  │  └─ oligon-brand (ggplot2_adapter) → plotting-libraries
│  │
│  ├─ PowerPoint presentation?
│  │  └─ oligon-brand (pptx_adapter) → scientific-slides → document-skills/pptx
│  │
│  ├─ Word document?
│  │  └─ oligon-brand (docx_adapter) → document-skills/docx
│  │
│  └─ PDF report?
│     └─ oligon-brand (reportlab_adapter) → markdown-to-pdf
```

### Brand Application: Quick Selection

| Output | Brand Skill | Implementation Skill |
|--------|-------------|---------------------|
| Publication figure (Python) | oligon-brand | plotting-libraries |
| Publication figure (R) | oligon-brand | plotting-libraries |
| Conference poster | oligon-brand | latex-posters / pptx-posters |
| Presentation slides | oligon-brand | scientific-slides |
| Word document | oligon-brand | document-skills/docx |
| PDF report | oligon-brand | markdown-to-pdf |

---
```

#### Step 2.2: Update Skill Count

Change footer from:
```
*This router provides deterministic skill selection for the 23-skill scientific writing library (19 top-level + 4 document-skills sub-skills).*
```

To:
```
*This router provides deterministic skill selection for the 24-skill scientific writing library (20 top-level + 4 document-skills sub-skills).*
```

#### Step 2.3: Update Quick Reference Matrix

Add to "Task → Skill Mapping" table:

```markdown
| **Brand** | | |
| Apply Oligon brand colors | oligon-brand | + plotting-libraries, document-skills |
| Get brand color hex codes | oligon-brand | — |
| Style matplotlib for brand | oligon-brand | + plotting-libraries |
```

---

### Phase 3: Cross-Reference Updates

#### Step 3.1: Update visual-design/SKILL.md

Replace current cross_references section:

```markdown
<cross_references>
## Related Skills

| Skill | Relationship |
|-------|-------------|
| `plotting-libraries` | Use for matplotlib/seaborn/R implementation after design decisions made |
| `scientific-schematics` | Use for diagram generation; visual-design provides principles |
| `generate-image` | Use for photorealistic images; visual-design provides aesthetics guidance |
| `scientific-slides` | Visual-design principles apply to presentation design |
| `latex-posters` / `pptx-posters` | Visual-design principles apply to poster design |
| `oligon-brand` | **Brand implementation** — provides Oligon-specific colors, typography, and logos. For programmatic access, load `oligon-brand/tokens/brand-tokens.json`. For matplotlib, use `oligon-brand/adapters/matplotlib_adapter.py`. |

**Skill Selection:**
See `SKILL_ROUTER.md` for decision trees when multiple skills may apply.

**Brand Integration:**
When Oligon brand standards apply:
1. Load `oligon-brand` tokens first
2. Apply `visual-design` principles
3. Implement with `plotting-libraries` or format-specific skill

**Typical Workflow:**
1. `visual-design` → establish design direction and specifications
2. `oligon-brand` → load brand colors and apply style (if applicable)
3. `plotting-libraries` or `scientific-schematics` → implement visuals
4. `visual-design` → review quality before submission

</cross_references>
```

Update references section path:

Change:
```markdown
| `references/BRAND_COLORS_v4.md` | Oligon brand color palette definitions |
```

To:
```markdown
| `../oligon-brand/references/brand-colors-full.md` | Oligon brand color palette definitions (V4) |
```

#### Step 3.2: Update plotting-libraries/SKILL.md

Replace `<styling_integration>` section:

```markdown
<styling_integration>
## Brand Integration

For scientific figures requiring organizational brand compliance, load brand tokens before plotting.

### Oligon Brand (recommended)

**Python (matplotlib/seaborn):**
```python
# Option 1: Use the adapter (recommended)
import sys
sys.path.insert(0, 'skills/oligon-brand')
from adapters.matplotlib_adapter import set_brand_style, BRAND_COLORS, get_cycle

# Apply brand defaults
set_brand_style('treatment_control')

# Access specific colors
highlight_color = BRAND_COLORS['brand_blue']  # #2DB2E8

# Or use specific cycle
colors = get_cycle('opposing')  # ['#2DB2E8', '#E8622D', '#222222']
```

**Option 2: Load tokens directly:**
```python
import json
from pathlib import Path

tokens_path = Path('skills/oligon-brand/tokens/brand-tokens.json')
with open(tokens_path) as f:
    brand = json.load(f)

# Access colors
brand_blue = brand['colors']['primary']['brand_blue']['hex']  # #2DB2E8
cycle = brand['color_cycles']['treatment_control']['colors']
```

**Option 3: Use mplstyle file:**
```python
import matplotlib.pyplot as plt

# If installed to matplotlib config
plt.style.use('oligon_color_brand')

# Or load directly
plt.style.use('skills/oligon-brand/assets/oligon_color_brand.mplstyle')
```

**R (ggplot2):**
```r
source('skills/oligon-brand/adapters/ggplot2_adapter.R')

# Apply brand theme
ggplot(df, aes(x, y, color = group)) +
  geom_point() +
  scale_color_oligon('treatment_control') +
  theme_oligon()
```

### Generic Brand Integration

For non-Oligon projects or custom brands:

```python
# Define your own brand colors
BRAND_COLORS = {
    'primary': '#YOUR_COLOR',
    'secondary': '#YOUR_COLOR',
    'accent': '#YOUR_COLOR',
}

# Apply to seaborn
import seaborn as sns
sns.set_palette(list(BRAND_COLORS.values()))
```

See `oligon-brand` skill for complete brand implementation patterns.
</styling_integration>
```

#### Step 3.3: Check Other Skills

Search for references to brand colors in other skills and update as needed:

- `markdown-to-pdf` — may reference brand colors for PDF styling
- `scientific-slides` — may have brand color references
- `latex-posters` / `pptx-posters` — may have brand color references

---

### Phase 4: Future-Proofing

#### Step 4.1: Document Brand Skill Interface Contract

Create `skills/BRAND_SKILL_INTERFACE.md`:

```markdown
# Brand Skill Interface Contract

> Specification for creating alternative brand skills compatible with the skill library
> Version: 1.0.0

## Purpose

This document defines the interface that brand skills must implement to integrate with the skill library. Following this contract ensures that:

1. Skills like `plotting-libraries` can load any brand without code changes
2. SKILL_ROUTER.md can route to the correct brand skill
3. Future brands integrate seamlessly

## Required Structure

A brand skill must have this directory structure:

```
skills/<brand-name>/
├── SKILL.md                    # Skill definition (required)
├── tokens/
│   └── brand-tokens.json       # Machine-readable tokens (required)
├── adapters/                   # Format adapters (at least one required)
│   ├── matplotlib_adapter.py   # Python/matplotlib
│   ├── ggplot2_adapter.R       # R/ggplot2
│   └── ...
├── references/
│   └── brand-spec.md           # Human-readable specification
└── assets/                     # Logo files, style files
    └── ...
```

## Required: brand-tokens.json Schema

The tokens file must include at minimum:

```json
{
  "meta": {
    "name": "Brand Name",
    "version": "X.Y.Z"
  },
  "colors": {
    "primary": {
      "<primary_color_name>": {"hex": "#XXXXXX", "role": "description"}
    },
    "neutrals": {
      "<neutral_1>": {"hex": "#XXXXXX", "role": "description"},
      "<neutral_2>": {"hex": "#XXXXXX", "role": "description"}
    },
    "structure": {
      "white": {"hex": "#FFFFFF"},
      "black": {"hex": "#000000"}
    }
  },
  "color_cycles": {
    "default": {
      "colors": ["#XXXXXX", "#XXXXXX", "#XXXXXX"]
    }
  }
}
```

## Required: Matplotlib Adapter API

If providing a Python/matplotlib adapter, it must export:

```python
# Constants
BRAND_COLORS: dict[str, str]  # Color name → hex mapping

# Functions
def set_brand_style(cycle: str = 'default') -> None:
    """Apply brand rcParams to matplotlib."""
    ...

def get_cycle(name: str) -> list[str]:
    """Return color cycle by name."""
    ...

def get_color(name: str) -> str:
    """Return single color by name."""
    ...
```

## Required: SKILL.md Frontmatter

```yaml
---
name: <brand-name>
version: X.Y.Z
description: "..."
brand-type: organizational  # or 'personal', 'institutional'
tokens-file: tokens/brand-tokens.json
---
```

## Integration Checklist

Before submitting a new brand skill:

- [ ] `tokens/brand-tokens.json` follows schema above
- [ ] At least one adapter provided
- [ ] Adapter exports required API
- [ ] SKILL.md has proper frontmatter
- [ ] SKILL.md has `<cross_references>` to visual-design and plotting-libraries
- [ ] Colors tested for accessibility (grayscale, colorblind)
- [ ] Added to SKILL_ROUTER.md "Brand Application" section
```

#### Step 4.2: Update SKILL_ROUTER.md for Multiple Brands

The "Brand Application" decision tree already supports this with:

```
├─ Which brand?
│  │
│  ├─ Oligon brand?
│  │  └─ YES → oligon-brand
│  │
│  ├─ [Future brand]?
│  │  └─ YES → [future-brand-skill]
│  │
│  └─ Other/None?
│     └─ Use visual-design defaults
```

---

### Phase 5: Validation & Documentation

#### Step 5.1: Test Adapter Imports

Create test script `scripts/test_oligon_brand_integration.py`:

```python
#!/usr/bin/env python3
"""Test oligon-brand skill integration."""

import sys
from pathlib import Path

# Add skill to path
skill_path = Path(__file__).parent.parent / "skills" / "oligon-brand"
sys.path.insert(0, str(skill_path))

def test_matplotlib_adapter():
    """Test matplotlib adapter imports and functions."""
    from adapters.matplotlib_adapter import (
        set_brand_style,
        BRAND_COLORS,
        get_cycle,
        get_color,
        get_divergent_cmap,
        add_panel_label,
    )

    # Test color access
    assert BRAND_COLORS['brand_blue'] == '#2DB2E8'
    assert get_color('brand_blue') == '#2DB2E8'

    # Test cycles
    cycle = get_cycle('treatment_control')
    assert len(cycle) == 3
    assert cycle[1] == '#2DB2E8'  # Treatment is Brand Blue

    # Test style application
    set_brand_style('treatment_control')

    import matplotlib.pyplot as plt
    assert plt.rcParams['axes.spines.top'] == False
    assert plt.rcParams['axes.spines.right'] == False

    print("✅ matplotlib_adapter: All tests passed")

def test_tokens_loading():
    """Test direct token loading."""
    import json

    tokens_path = skill_path / "tokens" / "brand-tokens.json"
    assert tokens_path.exists(), f"Tokens file not found: {tokens_path}"

    with open(tokens_path) as f:
        tokens = json.load(f)

    assert 'colors' in tokens
    assert 'color_cycles' in tokens
    assert tokens['colors']['primary']['brand_blue']['hex'] == '#2DB2E8'

    print("✅ brand-tokens.json: All tests passed")

def test_mplstyle():
    """Test mplstyle file exists and is valid."""
    mplstyle_path = skill_path / "assets" / "oligon_color_brand.mplstyle"
    assert mplstyle_path.exists(), f"mplstyle not found: {mplstyle_path}"

    # Try to load it
    import matplotlib.pyplot as plt
    plt.style.use(str(mplstyle_path))

    print("✅ oligon_color_brand.mplstyle: Loaded successfully")

if __name__ == '__main__':
    test_tokens_loading()
    test_matplotlib_adapter()
    test_mplstyle()
    print("\n✅ All oligon-brand integration tests passed!")
```

#### Step 5.2: Verify Cross-Reference Paths

Check these paths resolve correctly:

| From | Reference | Target |
|------|-----------|--------|
| visual-design/SKILL.md | `../oligon-brand/references/brand-colors-full.md` | ✓ Verify exists |
| plotting-libraries/SKILL.md | `skills/oligon-brand/tokens/brand-tokens.json` | ✓ Verify exists |
| SKILL_ROUTER.md | `oligon-brand` skill | ✓ Verify SKILL.md exists |

#### Step 5.3: Update CHANGELOG.md

Add to Unreleased section:

```markdown
### Added
- `oligon-brand` skill integration into skill library
  - Restructured SKILL.md to standard Claude skill format with XML sections
  - Added to SKILL_ROUTER.md with "Brand Application" decision tree
  - Updated cross-references in visual-design and plotting-libraries
  - Created BRAND_SKILL_INTERFACE.md for future brand skill development
  - Added integration test script

### Changed
- Updated skill count to 24 (20 top-level + 4 document sub-skills)
- visual-design now explicitly references oligon-brand for brand implementation
- plotting-libraries now includes concrete oligon-brand usage examples
```

---

## Execution Order

| Step | Description | Dependencies | Est. Lines Changed |
|------|-------------|--------------|-------------------|
| 1.1 | Create new SKILL.md frontmatter | None | ~10 |
| 1.2 | Add XML sections to SKILL.md | 1.1 | ~300 |
| 1.3 | Migrate content to proper sections | 1.2 | ~-100 (removal) |
| 1.4 | Verify reference paths | 1.3 | ~0 |
| 2.1 | Add brand routing to SKILL_ROUTER.md | 1.4 | ~50 |
| 2.2 | Update skill count | 2.1 | ~1 |
| 2.3 | Add to Quick Reference Matrix | 2.2 | ~5 |
| 3.1 | Update visual-design cross-refs | 2.3 | ~30 |
| 3.2 | Update plotting-libraries examples | 3.1 | ~60 |
| 3.3 | Check other skill references | 3.2 | ~10 |
| 4.1 | Create BRAND_SKILL_INTERFACE.md | 3.3 | ~150 (new file) |
| 4.2 | Verify router supports future brands | 4.1 | ~0 |
| 5.1 | Create integration test script | 4.2 | ~60 (new file) |
| 5.2 | Verify all paths | 5.1 | ~0 |
| 5.3 | Update CHANGELOG.md | 5.2 | ~15 |

**Total estimated changes:** ~590 lines added/modified

---

## Rollback Plan

If integration causes issues:

1. Restore original `skills/oligon-brand/SKILL.md` from git
2. Revert SKILL_ROUTER.md changes
3. Revert cross-reference updates in visual-design and plotting-libraries
4. Delete BRAND_SKILL_INTERFACE.md

---

## Success Criteria

Integration is complete when:

- [ ] `skills/oligon-brand/SKILL.md` follows standard skill XML structure
- [ ] SKILL_ROUTER.md includes oligon-brand in "Brand Application" section
- [ ] visual-design references oligon-brand with correct path
- [ ] plotting-libraries shows concrete oligon-brand usage examples
- [ ] Integration test script passes
- [ ] CHANGELOG.md documents all changes
- [ ] All cross-reference paths resolve correctly

---

*Plan created: 2025-12-30*
*Status: Ready for implementation*
