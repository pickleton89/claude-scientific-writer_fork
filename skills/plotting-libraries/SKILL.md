---
name: plotting-libraries
version: 1.1.0
description: "Data visualization for scientific research. Decision frameworks for Python (matplotlib, seaborn) and R (ggplot2, Bioconductor), common patterns, and bioinformatics-specific plots (volcano, heatmaps, survival curves, genome tracks). Use when creating plots, charts, figures, or visualizations."
extends: visual-design
---

> **Design Foundation:** This skill extends `visual-design` for design philosophy, typography principles, color theory, and publication specifications. See [visual-design](../visual-design/SKILL.md) for the foundational design layer.

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

<statistical_visualization>
## Visualizing Statistical Results

Coordinate with **statistical-analysis** skill for test selection. This section provides implementation patterns for common statistical visualizations.

### Group Comparisons (t-test, ANOVA)

**Python:**
```python
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

fig, ax = plt.subplots(figsize=(8, 6))

# Box plot with individual data points (preferred over bar plots)
sns.boxplot(data=df, x='group', y='value', ax=ax, width=0.5)
sns.stripplot(data=df, x='group', y='value', ax=ax, color='black', alpha=0.5, size=4)

# Add significance annotation
# For automated significance bars, use statannotations package
ax.set_ylabel('Measurement (units)')
ax.set_title('Group Comparison')
```

**R:**
```r
library(ggplot2)
library(ggpubr)

ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  stat_compare_means(method = "t.test") +  # or "anova", "wilcox.test"
  theme_classic() +
  theme(legend.position = "none")
```

### Normality Assessment (Q-Q Plots)

**Python:**
```python
from scipy import stats
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Histogram with normal curve overlay
axes[0].hist(data, bins=30, density=True, alpha=0.7)
x = np.linspace(data.min(), data.max(), 100)
axes[0].plot(x, stats.norm.pdf(x, data.mean(), data.std()), 'r-', lw=2)
axes[0].set_title('Distribution')

# Q-Q plot
stats.probplot(data, dist="norm", plot=axes[1])
axes[1].set_title('Q-Q Plot')
```

**R:**
```r
library(ggplot2)
library(patchwork)

p1 <- ggplot(df, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue") +
  stat_function(fun = dnorm, args = list(mean = mean(df$value), sd = sd(df$value)),
                color = "red", linewidth = 1)

p2 <- ggplot(df, aes(sample = value)) +
  stat_qq() + stat_qq_line(color = "red")

p1 + p2
```

### Regression Diagnostics

**Python:**
```python
import statsmodels.api as sm
import matplotlib.pyplot as plt

# Fit model
model = sm.OLS(y, sm.add_constant(X)).fit()

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Residuals vs Fitted
axes[0, 0].scatter(model.fittedvalues, model.resid, alpha=0.5)
axes[0, 0].axhline(y=0, color='r', linestyle='--')
axes[0, 0].set_xlabel('Fitted values')
axes[0, 0].set_ylabel('Residuals')

# Q-Q plot of residuals
sm.qqplot(model.resid, line='45', ax=axes[0, 1])

# Scale-Location
axes[1, 0].scatter(model.fittedvalues, np.sqrt(np.abs(model.resid)), alpha=0.5)
axes[1, 0].set_xlabel('Fitted values')
axes[1, 0].set_ylabel('√|Residuals|')

# Residuals vs Leverage
sm.graphics.influence_plot(model, ax=axes[1, 1])

plt.tight_layout()
```

**R:**
```r
# Base R diagnostic plots (4 classic plots)
par(mfrow = c(2, 2))
plot(model)

# ggplot2 version with ggfortify
library(ggfortify)
autoplot(model, which = 1:4)
```

### Correlation Matrices

**Python:**
```python
import seaborn as sns

# Compute correlation matrix
corr = df.corr()

# Create mask for upper triangle (avoid redundancy)
mask = np.triu(np.ones_like(corr, dtype=bool))

fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(corr, mask=mask, annot=True, fmt='.2f', cmap='RdBu_r',
            center=0, vmin=-1, vmax=1, ax=ax)
ax.set_title('Correlation Matrix')
```

**R:**
```r
library(corrplot)

cor_matrix <- cor(df, use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "color", type = "lower",
         addCoef.col = "black", number.cex = 0.7,
         tl.col = "black", tl.srt = 45)
```

### Forest Plots (Effect Sizes)

**Python:**
```python
import matplotlib.pyplot as plt

# Data: estimates, lower CI, upper CI, labels
fig, ax = plt.subplots(figsize=(8, 6))

y_pos = range(len(labels))
ax.errorbar(estimates, y_pos, xerr=[estimates - lower_ci, upper_ci - estimates],
            fmt='o', capsize=5, color='steelblue')
ax.axvline(x=0, color='gray', linestyle='--')  # or x=1 for ratios
ax.set_yticks(y_pos)
ax.set_yticklabels(labels)
ax.set_xlabel('Effect Size (95% CI)')
```

**R:**
```r
library(ggplot2)

ggplot(results, aes(x = estimate, y = reorder(variable, estimate))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Effect Size (95% CI)", y = NULL) +
  theme_minimal()
```

### Recommended Packages

| Task | Python | R |
|------|--------|---|
| Significance annotations | `statannotations` | `ggpubr::stat_compare_means` |
| Regression diagnostics | `statsmodels.graphics` | `ggfortify::autoplot` |
| Correlation plots | `seaborn.heatmap` | `corrplot`, `ggcorrplot` |
| Forest plots | `forestplot` | `forestplot`, `ggforestplot` |
| Effect size visualization | Custom matplotlib | `effectsize` + ggplot2 |

</statistical_visualization>

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
## Related Skills

**Skill Selection:** See `SKILL_ROUTER.md` for decision trees when multiple skills may apply.

- **[visual-design](../visual-design/SKILL.md)**: Design principles and publication specifications
- **[statistical-analysis](../statistical-analysis/SKILL.md)**: Statistics underlying visualizations—see `<statistical_visualization>` above for implementation patterns
- **[scientific-writing](../scientific-writing/SKILL.md)**: Figure legends and methods descriptions
- **[reproducible-research](../reproducible-research/SKILL.md)**: Environment management for reproducible figures
- **[code-documentation](../code-documentation/SKILL.md)**: Docstring patterns for plotting functions—see `visualization_documentation` section for templates
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
