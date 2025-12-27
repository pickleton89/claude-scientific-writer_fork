---
name: plotting-libraries
description: "Data visualization for scientific research. Decision frameworks for Python (matplotlib, seaborn) and R (ggplot2, Bioconductor), common patterns, and bioinformatics-specific plots (volcano, heatmaps, survival curves, genome tracks). Use when creating plots, charts, figures, or visualizations."
---

<objective>
Guide effective use of data visualization libraries for creating publication-quality scientific figures. This skill provides decision frameworks for choosing between Python and R libraries, common patterns for bioinformatics visualizations, and links to detailed reference documentation.
</objective>

<quick_start>
**Choose your ecosystem:**

| Preference | Ecosystem | Primary Libraries |
|------------|-----------|-------------------|
| Python-centric workflow | Python | matplotlib, seaborn |
| Statistical graphics, publications | R | ggplot2, Bioconductor |
| Bioinformatics (omics data) | Either | See specialized guides below |

**Quick patterns:**

Python (matplotlib/seaborn):
```python
import seaborn as sns
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10, 6))
sns.boxplot(data=df, x='group', y='value', ax=ax)
plt.savefig('figure.pdf', dpi=300, bbox_inches='tight')
```

R (ggplot2):
```r
library(ggplot2)

ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  theme_classic()
ggsave("figure.pdf", width = 7, height = 5)
```
</quick_start>

<decision_framework>
## Python vs R Decision

| Factor | Python | R |
|--------|--------|---|
| Existing codebase | Already using Python | Already using R |
| Statistical plots | seaborn adequate | ggplot2 preferred |
| Publication quality | Good with effort | Excellent defaults |
| Bioinformatics packages | pyGenomeTracks, scanpy | Bioconductor (ComplexHeatmap, etc.) |
| Learning curve | Lower for programmers | Lower for statisticians |

## Python Library Selection

**Start with seaborn when:**
- Working with pandas DataFrames
- Need statistical estimation (means, confidence intervals, regression)
- Creating distribution plots (histplot, kdeplot, violinplot)
- Want faceted/small-multiple layouts quickly
- Exploring data (pairplot, jointplot)

**Start with matplotlib when:**
- Need pixel-level control over layout
- Creating non-standard plot types
- Building complex multi-panel figures with different plot types
- Animating visualizations
- Embedding in GUI applications

**Use both together when:**
- seaborn creates the base plot, matplotlib fine-tunes it
- Need statistical analysis AND custom annotations
- Creating publication-quality figures

## R Library Selection

**ggplot2 for:**
- Most statistical visualizations
- Publication-quality figures
- Consistent grammar of graphics approach
- Extensive customization via themes

**Bioconductor packages for:**
- ComplexHeatmap: Advanced expression heatmaps
- EnhancedVolcano: Differential expression volcano plots
- survminer: Kaplan-Meier survival curves
- Gviz: Genome track visualization
- clusterProfiler/enrichplot: Enrichment analysis plots

</decision_framework>

<bioinformatics_plots>
## Specialized Bioinformatics Visualizations

| Plot Type | Python | R | Reference |
|-----------|--------|---|-----------|
| Volcano plot | matplotlib + adjustText | EnhancedVolcano, ggplot2 | [volcano_plots.md](references/volcano_plots.md) |
| Expression heatmap | seaborn clustermap | ComplexHeatmap, pheatmap | [heatmaps.md](references/heatmaps.md) |
| Survival curves | lifelines | survminer | [survival_curves.md](references/survival_curves.md) |
| Genome tracks | pyGenomeTracks | Gviz | [genome_tracks.md](references/genome_tracks.md) |
| PCA/UMAP | matplotlib, scanpy | ggplot2 | See respective guides |

**When to use which:**

| Visualization Need | Recommended Tool |
|-------------------|------------------|
| Differential expression results | EnhancedVolcano (R) or volcano_plots.md patterns (Python) |
| Gene expression clustering | ComplexHeatmap (R) for publication, seaborn clustermap for exploration |
| Clinical survival analysis | survminer (R) for risk tables |
| Multi-omics integration | ComplexHeatmap (R) for multiple linked heatmaps |
| Genome browser views | pyGenomeTracks (Python) or Gviz (R) |

</bioinformatics_plots>

<common_patterns>
**Combined Python workflow (publications):**
```python
import seaborn as sns
import matplotlib.pyplot as plt

# seaborn for statistical plot
fig, ax = plt.subplots(figsize=(8, 6))
sns.boxplot(data=df, x='group', y='value', ax=ax)

# matplotlib for customization
ax.set_title('Distribution by Group', fontsize=14, fontweight='bold')
ax.set_xlabel('Experimental Group', fontsize=12)
ax.set_ylabel('Measurement (units)', fontsize=12)

# Add custom elements
ax.axhline(y=threshold, color='red', linestyle='--', label='Threshold')
ax.legend()

plt.tight_layout()
plt.savefig('figure.pdf', dpi=300, bbox_inches='tight')
```

**R publication workflow:**
```r
library(ggplot2)
library(ggpubr)

ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "anova") +
  labs(x = "Group", y = "Value") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("figure.pdf", width = 7, height = 5)
```

**Export for publications:**
```python
# Vector format (scalable, required by most journals)
plt.savefig('figure.pdf', bbox_inches='tight')
plt.savefig('figure.svg', bbox_inches='tight')

# Raster format (presentations, web)
plt.savefig('figure.png', dpi=300, bbox_inches='tight')
```

```r
# Vector format
ggsave("figure.pdf", width = 180, height = 120, units = "mm")

# Raster format
ggsave("figure.png", width = 7, height = 5, dpi = 300)
```
</common_patterns>

<styling_integration>
For scientific figures, coordinate with:

- **[visual-design](../visual-design/SKILL.md)**: Design principles, publication requirements, brand colors

**Apply consistent styling:**

Python:
```python
from your_brand_module import COLORS  # See visual-design skill

# matplotlib
ax.plot(x, y, color=COLORS['primary'])

# seaborn
sns.set_palette([COLORS['primary'], COLORS['secondary'], COLORS['accent']])
```

R:
```r
# ggsci for journal palettes
library(ggsci)
ggplot(...) + scale_color_npg()  # Nature Publishing Group
ggplot(...) + scale_color_lancet()  # Lancet

# Manual colors
scale_fill_manual(values = c("Control" = "#377EB8", "Treatment" = "#E41A1C"))
```
</styling_integration>

<reference_guides>
## Detailed Library Documentation

### Python
- **[references/matplotlib.md](references/matplotlib.md)**: Complete matplotlib reference—OO API, plot types, subplots, 3D, styling, saving
- **[references/seaborn.md](references/seaborn.md)**: Complete seaborn reference—statistical plots, figure-level vs axes-level, palettes, theming

### R
- **[references/ggplot2.md](references/ggplot2.md)**: Complete ggplot2 reference—grammar of graphics, geoms, faceting, themes, publication export
- **[references/bioconductor_viz.md](references/bioconductor_viz.md)**: Bioconductor packages—ComplexHeatmap, EnhancedVolcano, ggpubr, clusterProfiler

### Bioinformatics-Specific
- **[references/volcano_plots.md](references/volcano_plots.md)**: Volcano plots for differential expression (Python + R)
- **[references/heatmaps.md](references/heatmaps.md)**: Expression heatmaps with clustering (Python + R)
- **[references/survival_curves.md](references/survival_curves.md)**: Kaplan-Meier and survival analysis (Python lifelines + R survminer)
- **[references/genome_tracks.md](references/genome_tracks.md)**: Genome visualization (pyGenomeTracks + Gviz)
</reference_guides>

<cross_references>
- **[visual-design](../visual-design/SKILL.md)**: Design principles and publication specifications
- **[statistical-analysis](../statistical-analysis/SKILL.md)**: Statistics underlying visualizations
- **[scientific-writing](../scientific-writing/SKILL.md)**: Figure legends and methods descriptions
- **[reproducible-research](../reproducible-research/SKILL.md)**: Environment management for reproducible figures
</cross_references>

<success_criteria>
- Appropriate library/ecosystem chosen for the task
- Figure renders correctly with all elements visible
- Exported at appropriate resolution (300 DPI for print, 150 for web)
- Labels, legends, and titles are readable
- Color palette is accessible (colorblind-friendly when possible)
- Vector format used for publications (PDF, SVG)
- Figure size matches journal requirements
</success_criteria>
