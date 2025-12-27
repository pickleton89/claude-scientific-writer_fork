---
name: code-documentation
description: "Documentation for scientific code. Docstrings (Python Google/NumPy, R roxygen2), README generation, notebook best practices, and inline commenting. Use when documenting analysis scripts, creating project READMEs, or writing reproducible notebooks."
---

<objective>
Generate clear, consistent documentation for scientific code that serves both human readers and automated tools. This skill covers function-level documentation (docstrings), project-level documentation (READMEs), and narrative documentation (notebooks).
</objective>

<scope>
**In Scope:**
- Docstring generation (Python: Google, NumPy, Sphinx styles; R: roxygen2)
- README templates for analysis projects, software tools, and data repositories
- Notebook best practices (Jupyter, R Markdown, Quarto)
- Inline commenting guidelines (when and how)
- API documentation for reusable packages

**Out of Scope:**
- Full documentation site generation (Sphinx/pkgdown setup and hosting)
- CI/CD integration for documentation builds
- Documentation translation or localization
</scope>

<documentation_hierarchy>
## Levels of Documentation

### 1. Inline Comments
- **When**: Complex logic that isn't self-evident
- **Style**: Brief, explains "why" not "what"
- **Anti-pattern**: Commenting obvious code (`i += 1  # increment i`)

### 2. Docstrings
- **When**: Every public function, class, module
- **Style**: Structured format (Google, NumPy, or roxygen2)
- **Content**: Purpose, parameters, returns, examples

### 3. README
- **When**: Every project/repository
- **Content**: Purpose, installation, usage, examples, citation

### 4. Notebooks
- **When**: Analysis workflows, tutorials, exploratory work
- **Style**: Narrative flow with code cells

### 5. API Documentation
- **When**: Reusable packages/libraries
- **Tools**: Sphinx (Python), pkgdown (R)
</documentation_hierarchy>

<style_selection>
## Which Docstring Style?

| Context | Recommended Style | Why |
|---------|-------------------|-----|
| Scientific Python | NumPy style | Standard in numpy, scipy, pandas ecosystem |
| General Python | Google style | Clean, readable, less verbose |
| Sphinx documentation | NumPy or Google | Both have Sphinx extensions |
| R packages | roxygen2 | Required for CRAN submission |
| Quick scripts | Minimal | One-line description sufficient |

## Quick Style Decision

```
Is this code for a package/library?
├── Yes → Is it R?
│         ├── Yes → roxygen2
│         └── No → NumPy style (scientific) or Google style (general)
└── No → Is it a one-off script?
          ├── Yes → Minimal docstrings (one-liners)
          └── No → Match the existing codebase style
```
</style_selection>

<readme_structure>
## README Essential Sections

### For Analysis Projects
1. **Title & Description**: What analysis this contains
2. **Data Requirements**: Input data sources and formats
3. **Environment Setup**: How to recreate the environment
4. **Usage**: How to run the analysis
5. **Output**: What files/figures are produced
6. **Citation**: How to cite this work

### For Software Tools
1. **Title & Badges**: Name, CI status, version, license
2. **Description**: What the tool does (one paragraph)
3. **Installation**: pip/conda/source instructions
4. **Quick Start**: Minimal working example
5. **Documentation**: Link to full docs
6. **Contributing**: How to contribute
7. **License & Citation**: Legal and academic credit

### For Data Repositories
1. **Title & DOI**: Dataset name and persistent identifier
2. **Description**: What data this contains
3. **Data Dictionary**: Column/field descriptions
4. **Collection Methods**: How data was generated
5. **Usage Terms**: License and citation requirements
6. **Related Publications**: Papers using this data
</readme_structure>

<notebook_guidelines>
## Notebook Best Practices

### Structure
1. **Title cell**: Clear, descriptive title as H1
2. **Purpose cell**: What this notebook accomplishes
3. **Setup cell**: All imports at the top
4. **Sections**: Logical divisions with markdown headers
5. **Conclusion cell**: Summary of findings

### Code Cells
- One concept per cell
- Cells should be runnable in order
- Use markdown cells to explain complex logic
- Keep cells under 20-30 lines

### Reproducibility
- Set random seeds at the start
- Print package versions (`%load_ext watermark`)
- Use relative paths or document data locations
- Clear outputs before committing (optional, depends on use case)

### Output Management
- Hide code for reports: `#| echo: false` (Quarto)
- Show code for methods: Full visibility
- Limit large outputs with `.head()` or sampling
</notebook_guidelines>

<usage_patterns>
## When to Apply This Skill

1. **Starting a new project**: "Create a README for my RNA-seq analysis project"
2. **Documenting functions**: "Write NumPy-style docstrings for these analysis functions"
3. **Code review preparation**: "Add appropriate documentation to this script"
4. **Publishing code**: "Prepare this code for a GitHub release with proper docs"
5. **Creating tutorials**: "Structure this notebook as a teaching resource"

## Integration with Other Skills

- After writing analysis code with `statistical-analysis`
- Before publishing with `reproducible-research`
- When preparing Methods sections with `scientific-writing`
</usage_patterns>

<common_pitfalls>
## Documentation Anti-Patterns

1. **Over-commenting**: Explaining what the code does line-by-line
2. **Under-documenting**: No docstrings on public functions
3. **Stale documentation**: Docs that don't match current behavior
4. **Missing examples**: Docstrings without usage examples
5. **Assumed knowledge**: Not defining domain-specific terms
6. **Inconsistent style**: Mixing Google and NumPy styles in one project

## Quality Signals

Good documentation:
- Can be understood without reading the implementation
- Includes examples that actually run
- Explains *why*, not just *what*
- Is updated when code changes
- Uses consistent formatting throughout
</common_pitfalls>

<cross_references>
- **reproducible-research**: Project structure and environment documentation
- **scientific-writing**: Methods section describing analysis
- **statistical-analysis**: Documenting statistical choices and parameters
- **plotting-libraries**: Documenting figure generation code
</cross_references>
