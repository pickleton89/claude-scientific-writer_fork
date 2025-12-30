# Scientific Visual Identity & Figure Standards

**Version 5.0 (Scientific Minimalist)**

This document outlines the strict visual identity for **figures, charts, diagrams, and presentation slides** for scientific communication (manuscripts, grants, and conferences).

Core Principle: **"Ascetic Precision."** Avoid chartjunk. Use Slate for structure and Teal for primary data. Color is used only to denote significance or identity, never for decoration. Maximize data-ink ratio with a modern "lab notebook" aesthetic.

**Always follow the approved color palette, layout rules, and typography standards below.** These rules apply to **all visuals**, including when generating **plotting code** (Matplotlib, ggplot2, Plotly, etc.).

---

## 1. Color Palettes (HEX + Role)

### Primary Palette (Teal & Slate)

The foundation of all scientific figures. Teal commands attention; Slate provides structure.

| Color | HEX | Role |
|-------|-----|------|
| Teal Primary | #0F766E | Primary Data / Treatment / Focus — the main action color |
| Slate Dark | #0F172A | Text / Axes — replaces pure black for softer contrast |
| Slate Medium | #64748B | Secondary Data / Context — control groups or baseline |
| Slate Light | #94A3B8 | Tertiary / Nonsignificant — de-emphasized elements |

### Background Colors

| Element | Color | HEX |
|---------|-------|-----|
| Page Background | Off-white (Slate-50) | #F8FAFC |
| Plot/Card Background | Pure White | #FFFFFF |
| Gridlines (if essential) | Slate-100 | #F1F5F9 |
| Axis lines | Slate-200 | #E2E8F0 |

> **Note:** Use off-white (#F8FAFC) for page/document backgrounds and pure white (#FFFFFF) for plot areas to create subtle depth and focus attention on data.

### Semantic Status Colors (Significance & Contrast)

Use these purposefully to encode scientific meaning. Each color has a specific semantic role.

| Color | HEX | Role & When to Use |
|-------|-----|-------------------|
| Indigo | #6366F1 | Highly Significant / Novelty — (p < 0.001), key findings |
| Emerald | #10B981 | Validation / Concordance — replication success, positive outcomes |
| Rose | #E11D48 | Negative / Warning / Discordance — adverse events, down-regulated, failures |
| Amber | #F59E0B | Limitations / Provisional — preprint status, caveats, uncertain results |

### Extended Slate Palette (Reference)

For fine-grained control over neutrals:

| Name | HEX | Use |
|------|-----|-----|
| Slate-50 | #F8FAFC | Page background |
| Slate-100 | #F1F5F9 | Gridlines, subtle borders |
| Slate-200 | #E2E8F0 | Axis lines, dividers |
| Slate-300 | #CBD5E1 | Disabled states |
| Slate-400 | #94A3B8 | Tertiary data (Slate Light) |
| Slate-500 | #64748B | Secondary data (Slate Medium) |
| Slate-600 | #475569 | Body text alternative |
| Slate-700 | #334155 | Emphasized text |
| Slate-800 | #1E293B | Headers |
| Slate-900 | #0F172A | Primary text (Slate Dark) |

---

## 2. Scenario-Specific Color Cycles

**Do not use default library color cycles.** Select the cycle that matches your scientific context:

### Cycle 1: Control vs. Treatment (Most Common)

Use when comparing an experimental condition against a baseline.

| Order | Color | HEX | Assigned To |
|-------|-------|-----|-------------|
| 1 | Slate Medium | #64748B | Control / Baseline |
| 2 | Teal Primary | #0F766E | Treatment / Experimental |
| 3 | Slate Light | #94A3B8 | Secondary comparison (if needed) |

```python
CYCLE_TREATMENT_CONTROL = ["#64748B", "#0F766E", "#94A3B8"]
```

### Cycle 2: Neutrals Only (No Designated Highlight)

Use when all groups are peers and no single condition deserves emphasis.

| Order | Color | HEX | Assigned To |
|-------|-------|-----|-------------|
| 1 | Slate Dark | #0F172A | Group A |
| 2 | Slate Medium | #64748B | Group B |
| 3 | Slate Light | #94A3B8 | Group C |

```python
CYCLE_NEUTRALS_ONLY = ["#0F172A", "#64748B", "#94A3B8"]
```

### Cycle 3: Opposing Effects (Teal vs. Rose)

Use when showing bidirectional changes (up/down, positive/negative, gain/loss).

| Order | Color | HEX | Assigned To |
|-------|-------|-----|-------------|
| 1 | Teal Primary | #0F766E | Direction A (e.g., Down-regulated, Positive) |
| 2 | Rose | #E11D48 | Direction B (e.g., Up-regulated, Negative) |
| 3 | Slate Medium | #64748B | Neutral / Control / Unchanged |

```python
CYCLE_OPPOSING = ["#0F766E", "#E11D48", "#64748B"]
```

### Cycle 4: Significance Levels

Use when encoding statistical significance or confidence levels.

| Order | Color | HEX | Assigned To |
|-------|-------|-----|-------------|
| 1 | Slate Light | #94A3B8 | Not Significant (NS) |
| 2 | Amber | #F59E0B | Marginal (p < 0.05) |
| 3 | Teal Primary | #0F766E | Significant (p < 0.01) |
| 4 | Indigo | #6366F1 | Highly Significant (p < 0.001) |

```python
CYCLE_SIGNIFICANCE = ["#94A3B8", "#F59E0B", "#0F766E", "#6366F1"]
```

### Cycle 5: Multiple Categories (≤4 groups, use sparingly)

Use only when you must distinguish 4 categories and cannot restructure the figure.

| Order | Color | HEX | Marker Shape |
|-------|-------|-----|--------------|
| 1 | Slate Dark | #0F172A | Circle (o) |
| 2 | Teal Primary | #0F766E | Square (s) |
| 3 | Slate Medium | #64748B | Triangle (^) |
| 4 | Indigo | #6366F1 | Diamond (D) |

```python
CYCLE_MULTI_CATEGORY = ["#0F172A", "#0F766E", "#64748B", "#6366F1"]
```

> **Important:** When using Cycle 5, marker shapes or line styles are **mandatory** for accessibility. If you need more than 4 categories, redesign the figure using faceting or small multiples.

---

## 3. Accessibility & Colorblind Safety

### Core Requirements

1. **All figures must remain interpretable in grayscale.** Test by converting to grayscale before submission.

2. **Do not rely solely on color** to distinguish data series. Combine color with:
   - Line style (solid, dashed, dotted)
   - Marker shape (circle, square, triangle, diamond)
   - Direct labels or annotations

3. **The Teal–Rose pair is colorblind-safe.** The combination of #0F766E and #E11D48 plus Slate grays provides good discrimination for deuteranopia and protanopia (red-green color blindness, ~8% of males).

4. **Indigo provides excellent contrast** against both Teal and Rose, making it suitable for highlighting the most significant findings.

### Testing Tools

Before finalizing figures, verify discriminability:

- [Coblis Colorblind Simulator](https://www.color-blindness.com/coblis-color-blindness-simulator/)
- [Viz Palette](https://projects.susielu.com/viz-palette)
- Matplotlib: `colorspacious` package
- ImageJ/FIJI: Built-in colorblind simulation

### Contrast Ratios

| Element | Minimum Contrast vs. White |
|---------|---------------------------|
| Text | 4.5:1 |
| Data points/lines | 3:1 |

All approved colors meet these thresholds against white backgrounds:
- Teal #0F766E: 5.4:1 ✓
- Slate Dark #0F172A: 16.1:1 ✓
- Slate Medium #64748B: 4.6:1 ✓
- Rose #E11D48: 4.5:1 ✓
- Indigo #6366F1: 4.5:1 ✓

> **Caution:** Slate Light (#94A3B8) has a contrast ratio of ~3:1 — use only for de-emphasized elements, not primary data.

---

## 4. Global Scientific Principles

1. **White Plot Backgrounds:** All plot areas use pure white (#FFFFFF). Page backgrounds may use off-white (#F8FAFC) for document/slide context.

2. **Slate Text & Structure:** Use Slate Dark (#0F172A) for text and axes instead of pure black for a modern, softer appearance.

3. **Data-Ink Ratio:** Maximize data visibility. Remove unnecessary borders, heavy fills, and distinct gridlines. Every pixel should serve the data.

4. **Neutrals First:** Do not use Teal (#0F766E) as the default color for *all* series. Most baseline data should be Slate; Teal highlights the most important comparison or effect. **Teal should appear in no more than one or two elements per panel.**

5. **No Effects:** Avoid gradients, shadows, or 3D elements.

6. **Redundant Encoding:** Never rely on color alone to convey critical information (see Accessibility Standards).

7. **Semantic Color Use:** Reserve Indigo for significance, Emerald for validation, Rose for warnings/negatives, Amber for limitations.

---

## 5. Typography Standards

Typography is functional. Use specific families for specific roles.

### Font Families by Role

| Role | Font Family | Characteristics |
|------|-------------|-----------------|
| UI / Headers | Inter, Roboto, Helvetica | Clean sans-serif. High legibility. |
| Data / Code / Tables | JetBrains Mono, Fira Code | Monospace. `tabular-nums`. Crucial for aligned numbers. |
| Long Text / Captions | Merriweather, Lora | Serif. Optimized for extended reading. |
| Fallback | Arial, DejaVu Sans | System fonts when primaries unavailable. |

### Font Sizes (Print/Manuscript)

| Element | Size |
|---------|------|
| Panel Labels (A, B, C) | 9–10 pt, **bold** |
| Panel Titles/Headings | 8–9 pt |
| Axis Labels | 8 pt |
| Axis Tick Labels | 7–8 pt |
| Legend Text | 7–8 pt |
| Annotations | 7 pt |

### Panel Labels (Multi-Panel Figures)

- **Format:** Uppercase letters (**A**, **B**, **C**)
- **Style:** Bold, Slate Dark text, 9–10 pt
- **Placement:** Top-left corner of each panel, aligned with the left edge of the y-axis

---

## 6. Output Resolution Standards

| Output Context | Resolution | File Format |
|----------------|------------|-------------|
| Print/Publication | 300 DPI minimum (600 DPI for line art) | TIFF, EPS, PDF |
| Screen/Presentation | 150 DPI | PNG, PDF |
| Web/Preprint | 150 DPI | PNG, SVG |

> **Vector formats (PDF, EPS, SVG) are preferred** for all line-based figures as they scale without quality loss.

---

## 7. Chart & Figure Design Standards

### Structure & Layout ("Despined")

Figures must use an "open" or "despined" layout common in high-impact journals.

- **Axes & Spines:** Only keep the **left** and **bottom** spines. **Remove top and right spines.**
- **Axis Lines:** Slate-200 (#E2E8F0), 0.8 pt weight
- **Gridlines:** Off by default. If essential for readability, use Slate-100 (#F1F5F9), 0.5 pt weight

### Aspect Ratio

| Chart Type | Recommended Ratio |
|------------|-------------------|
| Time Series / Line Charts | 4:3 (wider than tall) |
| Scatter Plots / Correlations | 1:1 (square) |
| Bar Charts | Context-dependent; avoid excessive width |
| Heatmaps | Determined by data dimensions |

### Legends

- **Placement:** Inside the plotting area in an empty region whenever possible
- **Style:** Minimalist, no surrounding frame (`frameon=False`)
- **Exception:** Use a very thin Slate-200 frame only if legend overlaps data

---

## 8. Color Selection Decision Tree

```
START: How many data series/categories?
│
├─► 1 series (single condition)
│   └─► Use #64748B (Slate Medium) — no highlight needed
│
├─► 2 series
│   │
│   ├─► Treatment vs. Control?
│   │   └─► Cycle 1: #64748B (Control) + #0F766E (Treatment)
│   │
│   └─► Opposing effects (Up vs. Down)?
│       └─► Cycle 3: #0F766E + #E11D48
│
├─► 3 series
│   │
│   ├─► One highlighted condition?
│   │   └─► Cycle 1: #64748B, #0F766E, #94A3B8
│   │
│   ├─► All conditions equal weight?
│   │   └─► Cycle 2: #0F172A, #64748B, #94A3B8
│   │
│   └─► Two opposing + one neutral?
│       └─► Cycle 3: #0F766E, #E11D48, #64748B
│
├─► 4 series
│   │
│   ├─► Significance levels?
│   │   └─► Cycle 4: #94A3B8, #F59E0B, #0F766E, #6366F1
│   │
│   └─► Categorical groups?
│       └─► Cycle 5: #0F172A, #0F766E, #64748B, #6366F1
│       └─► MANDATORY: Add marker shapes for accessibility
│
└─► 5+ series
    └─► STOP: Redesign the figure
        • Split into multiple panels
        • Use faceting/small multiples
        • Group categories hierarchically
```

---

## 9. Specific Chart Guidance

### Line & Scatter Plots

- **Line width:** 1.0–1.5 pt for primary data; 0.75 pt for secondary
- **Lines:** Solid for primary data; dashed/dotted for secondary (accessibility redundancy)
- **Markers:** Clean, defined shapes; size 4–6 pt
- **Marker edges:** Avoid heavy outlines unless needed for overlap contrast

### Error Representation Standards

| Error Type | When to Use | Visual Style |
|------------|-------------|--------------|
| **SEM** (Standard Error of Mean) | Estimating precision of the mean | Default for most experimental comparisons |
| **SD** (Standard Deviation) | Showing data spread/variability | When variability itself is of interest |
| **95% CI** (Confidence Interval) | Statistical inference | Preferred for clinical/translational data |

**Visual encoding:**

- **Discrete points (bar charts, dot plots):** Thin error bars (0.75 pt), with or without caps (be consistent)
- **Continuous data (line plots):** Semi-transparent fill bands (20–30% opacity) without borders
- **Always state in figure legend:** "Error bars represent [SEM/SD/95% CI], n = X"

### Bar Charts

- **Use Slate for most bars**; reserve Teal (#0F766E) or semantic colors for emphasis
- Use solid fills based on the appropriate cycle
- **Do not** outline bars with dark lines (reduces data-ink ratio)
- Consider dot plots or box plots as alternatives for small n

### Box Plots & Violin Plots

- Box/violin fill: 50% opacity of assigned color
- Median line: solid, same color at full opacity
- Outliers: same color, smaller marker size

### Heatmaps (Divergent Scales)

For data centered around zero (z-scores, log2 fold change):

| Value | Color | HEX |
|-------|-------|-----|
| Extreme Negative | Teal Primary | #0F766E |
| Moderate Negative | Teal (lighter) | #14B8A6 |
| Zero/Neutral | White | #FFFFFF |
| Moderate Positive | Rose (lighter) | #FB7185 |
| Extreme Positive | Rose | #E11D48 |

**Implementation requirements:**

1. **Use perceptually uniform interpolation** — linear RGB blending creates artifacts
2. **Center the colormap at zero** — do not let automatic scaling shift the neutral point
3. **Consider symmetric limits** — if data ranges -3 to +5, use limits of -5 to +5
4. **Add a colorbar** — always include with tick marks at meaningful intervals

---

## 10. Canonical Plot Templates for Biology

These patterns cover common visualization needs in biological and biomedical research.

### Kaplan-Meier Survival Curves

**Context:** Comparing survival between treatment groups (common in preclinical oncology).

| Element | Color | Style |
|---------|-------|-------|
| Control group | #64748B | Solid line |
| Treatment group | #0F766E | Solid line |
| Additional treatment | #6366F1 | Solid line (if showing enhanced effect) |
| Adverse outcome group | #E11D48 | Solid line (if showing negative effect) |
| Confidence bands | Same as line | 20% opacity fill, no border |
| Censored tick marks | #0F172A | Short vertical ticks on curve |
| At-risk table | #0F172A | Below x-axis, matching group colors for labels |

**Additional guidance:**

- Line width: 1.25–1.5 pt
- Include median survival lines (thin dashed Slate) if medians are reached
- P-value and hazard ratio: place in upper right or lower left, 7 pt text
- X-axis: "Time (days)" or appropriate unit; Y-axis: "Survival probability" or "Percent survival"

### Volcano Plots

**Context:** Differential expression analysis showing fold change vs. significance.

| Element | Color | HEX |
|---------|-------|-----|
| Non-significant points | Slate Light | #94A3B8 |
| Significant up-regulated | Rose | #E11D48 |
| Significant down-regulated | Teal Primary | #0F766E |
| Highly significant | Indigo | #6366F1 (optional overlay) |
| Threshold lines | Slate Medium | #64748B (dashed, 0.5 pt) |
| Labeled genes | Slate Dark text | #0F172A (7 pt, with leader lines if needed) |

**Additional guidance:**

- Point size: 3–4 pt for non-significant, 5–6 pt for significant (optional size encoding)
- Alpha: 0.6–0.8 for non-significant points to reduce overplotting
- X-axis: "log₂(Fold Change)"; Y-axis: "-log₁₀(p-value)" or "-log₁₀(FDR)"
- Symmetric x-axis limits recommended

### Grouped Bar Charts (Treatment Comparisons)

**Context:** Comparing measurements across conditions and groups.

| Element | Color | When |
|---------|-------|------|
| Control bars | #64748B | Default/baseline condition |
| Treatment bars | #0F766E | Experimental condition of interest |
| Validated/replicated | #10B981 | Confirmed findings |
| Negative/adverse | #E11D48 | Adverse outcomes or failures |

**Additional guidance:**

- Group bars by biological condition (e.g., cell line), not by treatment
- Use Slate for most bars; reserve Teal for the key comparison
- Error bars: thin (0.75 pt), same color as bar or Slate Dark
- Consider replacing with dot plots if n < 10

### Dose-Response Curves

**Context:** IC50/EC50 determination, drug sensitivity assays.

| Element | Color | Style |
|---------|-------|-------|
| Drug A (focus) | #0F766E | Solid line, circle markers |
| Drug B (comparison) | #64748B | Dashed line, square markers |
| Vehicle control | #94A3B8 | Dotted line (if shown as curve) |
| IC50 reference line | #CBD5E1 | Thin horizontal dashed line |

**Additional guidance:**

- X-axis: log scale concentration; Y-axis: "% Viability" or "% Response"
- Include 95% CI bands (20% opacity) for fitted curves
- Mark IC50 values with vertical droplines or annotations
- Error bars at each concentration: SEM with n stated

### Flow Cytometry / Population Comparisons

**Context:** Comparing cell populations or marker expression.

| Element | Color | When |
|---------|-------|------|
| Isotype/FMO control | #94A3B8 | Negative control (de-emphasized) |
| Untreated/baseline | #64748B | Reference population |
| Treatment condition | #0F766E | Experimental condition |
| Secondary marker | #E11D48 | If showing opposing populations |

**Additional guidance:**

- Histogram overlays: 50% opacity fills, no borders
- Order: controls behind, treatment in front
- Include gate boundaries as thin Slate Dark lines
- MFI values: annotate directly on plot or in legend

---

## 11. Presentation Slide Rules (Exceptions)

While publication standards are preferred, presentation slides allow slight deviations for visibility.

| Element | Publication Standard | Slide Exception |
|---------|---------------------|-----------------|
| Titles | Slate Dark only | May use #0F766E (Teal) for emphasis |
| Font sizes | 7–9 pt | Scale up 2–3× for readability |
| Accent elements | Minimal | Teal (#0F766E) for callouts, icons |
| Body text | Slate Dark | Slate Dark (#0F172A) |
| Page background | Off-white (#F8FAFC) | Off-white or white |
| Plot background | White | White |

---

## 12. Code Implementation

### Matplotlib Complete Setup

```python
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# =============================================================================
# SCIENTIFIC MINIMALIST COLOR DEFINITIONS
# =============================================================================

BRAND_COLORS = {
    # Primary
    'teal': '#0F766E',

    # Slate palette (neutrals)
    'slate_dark': '#0F172A',
    'slate_medium': '#64748B',
    'slate_light': '#94A3B8',
    'slate_100': '#F1F5F9',
    'slate_200': '#E2E8F0',
    'slate_50': '#F8FAFC',

    # Semantic status colors
    'indigo': '#6366F1',
    'emerald': '#10B981',
    'rose': '#E11D48',
    'amber': '#F59E0B',

    # Structure
    'white': '#FFFFFF',
}

# =============================================================================
# SCENARIO-SPECIFIC COLOR CYCLES
# =============================================================================

# Cycle 1: Control vs Treatment (most common)
CYCLE_TREATMENT_CONTROL = ['#64748B', '#0F766E', '#94A3B8']

# Cycle 2: Neutrals only (no designated highlight)
CYCLE_NEUTRALS_ONLY = ['#0F172A', '#64748B', '#94A3B8']

# Cycle 3: Opposing effects (teal vs rose)
CYCLE_OPPOSING = ['#0F766E', '#E11D48', '#64748B']

# Cycle 4: Significance levels
CYCLE_SIGNIFICANCE = ['#94A3B8', '#F59E0B', '#0F766E', '#6366F1']

# Cycle 5: Multiple categories (use sparingly, add shapes!)
CYCLE_MULTI_CATEGORY = ['#0F172A', '#0F766E', '#64748B', '#6366F1']

# Volcano plot specific
VOLCANO_COLORS = {
    'nonsig': '#94A3B8',
    'up': '#E11D48',
    'down': '#0F766E',
}

# =============================================================================
# RCPARAMS CONFIGURATION
# =============================================================================

def set_scientific_style(cycle='treatment_control'):
    """
    Apply Scientific Minimalist visual standards to matplotlib.

    Parameters
    ----------
    cycle : str
        Which color cycle to use as default:
        - 'treatment_control' (default)
        - 'neutrals_only'
        - 'opposing'
        - 'significance'
        - 'multi_category'
    """

    cycle_map = {
        'treatment_control': CYCLE_TREATMENT_CONTROL,
        'neutrals_only': CYCLE_NEUTRALS_ONLY,
        'opposing': CYCLE_OPPOSING,
        'significance': CYCLE_SIGNIFICANCE,
        'multi_category': CYCLE_MULTI_CATEGORY,
    }

    selected_cycle = cycle_map.get(cycle, CYCLE_TREATMENT_CONTROL)

    plt.rcParams.update({
        # Figure
        'figure.facecolor': '#F8FAFC',  # Off-white page background
        'figure.edgecolor': '#F8FAFC',
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.facecolor': '#F8FAFC',
        'savefig.edgecolor': '#F8FAFC',
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,

        # Axes
        'axes.facecolor': '#FFFFFF',  # White plot background
        'axes.edgecolor': '#E2E8F0',  # Slate-200
        'axes.linewidth': 0.8,
        'axes.grid': False,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.labelcolor': '#64748B',  # Slate-500
        'axes.labelsize': 8,
        'axes.titlesize': 9,
        'axes.titleweight': 'normal',
        'axes.prop_cycle': plt.cycler(color=selected_cycle),

        # Ticks
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'xtick.color': '#64748B',
        'ytick.color': '#64748B',
        'xtick.direction': 'out',
        'ytick.direction': 'out',

        # Grid (off by default)
        'grid.color': '#F1F5F9',  # Slate-100
        'grid.linewidth': 0.5,
        'grid.alpha': 1.0,

        # Legend
        'legend.fontsize': 7,
        'legend.frameon': False,
        'legend.loc': 'best',

        # Font
        'font.family': 'sans-serif',
        'font.sans-serif': ['Inter', 'Roboto', 'Arial', 'Helvetica', 'DejaVu Sans'],
        'font.monospace': ['JetBrains Mono', 'Fira Code', 'Courier New'],
        'font.size': 8,

        # Text
        'text.color': '#0F172A',  # Slate-900

        # Lines
        'lines.linewidth': 1.25,
        'lines.markersize': 5,

        # Patches (bars, etc.)
        'patch.edgecolor': 'none',
    })

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_divergent_cmap():
    """
    Create a perceptually uniform divergent colormap for heatmaps.
    Centered at white: teal (neg) -> white -> rose (pos)
    """
    from matplotlib.colors import LinearSegmentedColormap

    colors = [
        BRAND_COLORS['teal'],
        '#14B8A6',  # Teal-500 (lighter)
        BRAND_COLORS['white'],
        '#FB7185',  # Rose-400 (lighter)
        BRAND_COLORS['rose'],
    ]

    positions = [0.0, 0.25, 0.5, 0.75, 1.0]

    return LinearSegmentedColormap.from_list(
        'scientific_divergent',
        list(zip(positions, colors))
    )


def add_panel_label(ax, label, offset=(-0.12, 1.05)):
    """Add a panel label (A, B, C, etc.) to an axes."""
    ax.text(
        offset[0], offset[1],
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight='bold',
        va='top',
        ha='left',
        color=BRAND_COLORS['slate_dark']
    )


def add_error_band(ax, x, y, yerr, color, alpha=0.25, label=None):
    """Add a semi-transparent error band around a line."""
    if isinstance(yerr, tuple):
        lower, upper = y - yerr[0], y + yerr[1]
    else:
        lower, upper = y - yerr, y + yerr

    ax.fill_between(x, lower, upper, color=color, alpha=alpha, linewidth=0)
    line, = ax.plot(x, y, color=color, label=label)
    return line


def plot_volcano(ax, log2fc, pval, fc_thresh=1.0, p_thresh=0.05,
                 gene_labels=None, n_labels=10):
    """
    Create a volcano plot with Scientific Minimalist colors.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    log2fc : array-like
        Log2 fold change values
    pval : array-like
        P-values (will be -log10 transformed)
    fc_thresh : float
        Fold change threshold (absolute log2)
    p_thresh : float
        P-value threshold
    gene_labels : array-like, optional
        Gene names for labeling
    n_labels : int
        Number of top genes to label
    """
    neg_log_p = -np.log10(pval)

    # Classify points
    nonsig = (np.abs(log2fc) < fc_thresh) | (pval >= p_thresh)
    sig_up = (log2fc >= fc_thresh) & (pval < p_thresh)
    sig_down = (log2fc <= -fc_thresh) & (pval < p_thresh)

    # Plot
    ax.scatter(log2fc[nonsig], neg_log_p[nonsig],
               c=VOLCANO_COLORS['nonsig'], s=12, alpha=0.6, label='NS')
    ax.scatter(log2fc[sig_up], neg_log_p[sig_up],
               c=VOLCANO_COLORS['up'], s=20, alpha=0.8, label='Up')
    ax.scatter(log2fc[sig_down], neg_log_p[sig_down],
               c=VOLCANO_COLORS['down'], s=20, alpha=0.8, label='Down')

    # Threshold lines
    ax.axhline(-np.log10(p_thresh), color='#64748B', linestyle='--', linewidth=0.5)
    ax.axvline(-fc_thresh, color='#64748B', linestyle='--', linewidth=0.5)
    ax.axvline(fc_thresh, color='#64748B', linestyle='--', linewidth=0.5)

    ax.set_xlabel('log₂(Fold Change)')
    ax.set_ylabel('-log₁₀(p-value)')

    return ax

# =============================================================================
# EXAMPLE: KAPLAN-MEIER STYLE SURVIVAL PLOT
# =============================================================================

def plot_survival_example():
    """Example survival curve following Scientific Minimalist standards."""
    set_scientific_style('treatment_control')

    fig, ax = plt.subplots(figsize=(5, 3.75))

    # Simulated survival data
    time_ctrl = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40])
    surv_ctrl = np.array([1.0, 0.95, 0.85, 0.70, 0.55, 0.40, 0.30, 0.22, 0.18])

    time_tx = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40])
    surv_tx = np.array([1.0, 0.98, 0.94, 0.88, 0.80, 0.72, 0.65, 0.58, 0.52])

    # Plot with confidence bands
    ax.fill_between(time_ctrl, surv_ctrl - 0.08, surv_ctrl + 0.08,
                    color=BRAND_COLORS['slate_medium'], alpha=0.2)
    ax.step(time_ctrl, surv_ctrl, where='post',
            color=BRAND_COLORS['slate_medium'], linewidth=1.5, label='Control')

    ax.fill_between(time_tx, surv_tx - 0.06, surv_tx + 0.06,
                    color=BRAND_COLORS['teal'], alpha=0.2)
    ax.step(time_tx, surv_tx, where='post',
            color=BRAND_COLORS['teal'], linewidth=1.5, label='Treatment')

    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Survival probability')
    ax.set_ylim(0, 1.05)
    ax.legend(loc='lower left')
    add_panel_label(ax, 'A')

    # Add p-value annotation
    ax.text(0.95, 0.95, 'p = 0.003\nHR = 0.45', transform=ax.transAxes,
            ha='right', va='top', fontsize=7, color=BRAND_COLORS['slate_dark'])

    plt.tight_layout()
    return fig

if __name__ == "__main__":
    fig = plot_survival_example()
    plt.savefig('survival_example.png', dpi=300)
    plt.show()
```

### ggplot2 Complete Setup (R)

```r
library(ggplot2)

# =============================================================================
# SCIENTIFIC MINIMALIST COLOR DEFINITIONS
# =============================================================================

brand_colors <- list(
  # Primary
  teal = "#0F766E",

  # Slate palette
  slate_dark = "#0F172A",
  slate_medium = "#64748B",
  slate_light = "#94A3B8",
  slate_100 = "#F1F5F9",
  slate_200 = "#E2E8F0",
  slate_50 = "#F8FAFC",

  # Semantic status colors
  indigo = "#6366F1",
  emerald = "#10B981",
  rose = "#E11D48",
  amber = "#F59E0B",

  # Structure
  white = "#FFFFFF"
)

# =============================================================================
# SCENARIO-SPECIFIC SCALES
# =============================================================================

# Cycle 1: Treatment vs Control
scale_treatment_control <- c(
  "Control" = "#64748B",
  "Treatment" = "#0F766E"
)

# Cycle 2: Neutrals only
scale_neutrals <- c("#0F172A", "#64748B", "#94A3B8")

# Cycle 3: Opposing effects
scale_opposing <- c(
  "Down" = "#0F766E",
  "Up" = "#E11D48",
  "NS" = "#64748B"
)

# Cycle 4: Significance levels
scale_significance <- c(
  "NS" = "#94A3B8",
  "p<0.05" = "#F59E0B",
  "p<0.01" = "#0F766E",
  "p<0.001" = "#6366F1"
)

# Cycle 5: Multiple categories
scale_multi <- c("#0F172A", "#0F766E", "#64748B", "#6366F1")

# Volcano plot
scale_volcano <- c(
  "NS" = "#94A3B8",
  "Up" = "#E11D48",
  "Down" = "#0F766E"
)

# =============================================================================
# SCIENTIFIC MINIMALIST THEME
# =============================================================================

theme_scientific <- function(base_size = 8, base_family = "Inter") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Background
      plot.background = element_rect(fill = "#F8FAFC", color = NA),
      panel.background = element_rect(fill = "#FFFFFF", color = NA),

      # Grid
      panel.grid.major = element_line(color = "#F1F5F9", linewidth = 0.5),
      panel.grid.minor = element_blank(),

      # Axes
      axis.line = element_line(color = "#E2E8F0", linewidth = 0.8 / .pt),
      axis.ticks = element_line(color = "#E2E8F0", linewidth = 0.4 / .pt),
      axis.text = element_text(color = "#64748B", size = 7),
      axis.title = element_text(color = "#0F172A", size = 8),

      # Text
      text = element_text(color = "#0F172A"),
      plot.title = element_text(color = "#0F172A", size = 9, hjust = 0),
      plot.subtitle = element_text(color = "#64748B", size = 8, hjust = 0),

      # Legend
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 7, color = "#64748B"),
      legend.title = element_text(size = 8, color = "#0F172A"),

      # Facets
      strip.background = element_rect(fill = "#F1F5F9", color = NA),
      strip.text = element_text(color = "#0F172A", size = 8, face = "bold")
    )
}

# =============================================================================
# COLOR SCALES
# =============================================================================

scale_fill_scientific <- function(...) {
  scale_fill_manual(values = c("#0F766E", "#64748B", "#6366F1", "#E11D48"), ...)
}

scale_color_scientific <- function(...) {
  scale_color_manual(values = c("#0F766E", "#64748B", "#6366F1", "#E11D48"), ...)
}

# Divergent scale for heatmaps
scale_fill_scientific_divergent <- function(limits = NULL, ...) {
  scale_fill_gradientn(
    colors = c(
      "#0F766E",  # Teal (negative extreme)
      "#14B8A6",  # Teal-500
      "#FFFFFF",  # White (center)
      "#FB7185",  # Rose-400
      "#E11D48"   # Rose (positive extreme)
    ),
    values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
    limits = limits,
    ...
  )
}

# =============================================================================
# EXAMPLE: VOLCANO PLOT
# =============================================================================

# Simulated data
set.seed(42)
volcano_df <- data.frame(
  log2FC = rnorm(1000, 0, 1.5),
  pval = 10^(-runif(1000, 0, 5))
)
volcano_df$category <- with(volcano_df, ifelse(
  abs(log2FC) < 1 | pval > 0.05, "NS",
  ifelse(log2FC > 0, "Up", "Down")
))

ggplot(volcano_df, aes(x = log2FC, y = -log10(pval), color = category)) +
  geom_point(aes(size = category, alpha = category)) +
  scale_color_manual(values = scale_volcano) +
  scale_size_manual(values = c("NS" = 1, "Up" = 2, "Down" = 2)) +
  scale_alpha_manual(values = c("NS" = 0.5, "Up" = 0.8, "Down" = 0.8)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "#64748B", linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "#64748B", linewidth = 0.3) +
  labs(x = expression(log[2](Fold~Change)),
       y = expression(-log[10](p-value))) +
  theme_scientific() +
  theme(legend.position = "none")

ggsave("volcano_example.png", width = 4, height = 4, dpi = 300, bg = "#F8FAFC")
```

---

## 13. Quick Reference Card

### Essential Colors

```
PRIMARY:       #0F766E  (Teal)
SLATE DARK:    #0F172A  (Text/Primary)
SLATE MEDIUM:  #64748B  (Secondary)
SLATE LIGHT:   #94A3B8  (Tertiary/NS)
BACKGROUND:    #F8FAFC  (Page) / #FFFFFF (Plot)

SEMANTIC:
  INDIGO:      #6366F1  (Highly Significant)
  EMERALD:     #10B981  (Validation)
  ROSE:        #E11D48  (Warning/Negative)
  AMBER:       #F59E0B  (Limitations)
```

### Cycle Quick Reference

| Scenario | Colors (in order) |
|----------|-------------------|
| Control vs Treatment | #64748B, #0F766E |
| Neutrals only | #0F172A, #64748B, #94A3B8 |
| Opposing effects | #0F766E, #E11D48, #64748B |
| Significance levels | #94A3B8, #F59E0B, #0F766E, #6366F1 |
| 4 categories | #0F172A, #0F766E, #64748B, #6366F1 + shapes |

### Pre-Flight Checklist

Before submitting any figure:

- [ ] Off-white page background (#F8FAFC) or white
- [ ] White plot background (#FFFFFF)
- [ ] Slate Dark text (#0F172A)
- [ ] Top and right spines removed
- [ ] Correct color cycle selected for scientific context
- [ ] Teal limited to 1–2 elements per panel
- [ ] Colors match palette (no library defaults)
- [ ] Redundant encoding present (shapes/styles with color)
- [ ] Tested in grayscale
- [ ] Colorblind simulation checked
- [ ] Error bars labeled in legend (SEM/SD/95% CI, n = ?)
- [ ] Font is Inter/Arial, 7–9 pt range
- [ ] Resolution ≥300 DPI for print
- [ ] Panel labels (A, B, C) bold, 9–10 pt, top-left aligned

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | — | Initial release |
| 2.0 | — | Shifted to neutrals-first philosophy |
| 3.0 | — | Added accessibility standards, expanded code templates, error bar guidance, decision tree |
| 4.0 | — | Scenario-specific color cycles, canonical bio-plot templates (K-M, volcano, dose-response, flow), clarified supporting blue roles, grayscale testing requirement, enhanced bar chart guidance |
| 5.0 | — | Scientific Minimalist redesign: Teal/Slate palette, semantic status colors (Indigo, Emerald, Rose, Amber), off-white backgrounds, expanded typography (Inter, JetBrains Mono, Merriweather), significance-based color cycle |

---

*Document maintained for scientific figure standardization. For questions or proposed amendments, contact the visual standards committee.*
