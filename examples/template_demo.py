#!/usr/bin/env python3
"""
Template Demo - Demonstration of the markdown-to-pdf template system.

Shows how to:
1. List available templates
2. Parse markdown documents with frontmatter
3. Validate against schemas
4. Generate branded PDFs

Run with: uv run python examples/template_demo.py
"""

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Spacer, Paragraph, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle

from oligon_reports import (
    BRAND_COLORS,
    TemplateParser,
    FindingCard,
    StatusTable,
    CalloutBox,
    MetadataHeader,
    SectionDivider,
)


def demo_list_templates():
    """Demonstrate listing available templates."""
    print("\n" + "=" * 60)
    print("AVAILABLE TEMPLATES")
    print("=" * 60)

    templates = TemplateParser.list_templates()

    # Group by category
    by_category: dict[str, list[dict]] = {}
    for t in templates:
        cat = t.get("category", "unknown")
        by_category.setdefault(cat, []).append(t)

    for category, items in sorted(by_category.items()):
        print(f"\n{category.upper()}:")
        for item in items:
            print(f"  - {item['type']}: {item['description'][:60]}...")

    print(f"\nTotal: {len(templates)} templates across {len(by_category)} categories")
    return templates


def demo_parse_document():
    """Demonstrate parsing a markdown document."""
    print("\n" + "=" * 60)
    print("PARSING SAMPLE DOCUMENT")
    print("=" * 60)

    # Sample analysis report markdown
    sample_md = """---
type: analysis-report
title: RNA-seq Differential Expression Analysis
date: 2025-12-28
author: Research Team
project: Gene Expression Study
assessment: pass-fail
---

# Executive Summary

This analysis identifies differentially expressed genes between treatment and control groups using DESeq2.

## Key Metrics

| Metric | Value |
|--------|-------|
| Total samples | 24 |
| Genes tested | 18,432 |
| DE genes (FDR<0.05) | 342 |

## Quality Control

| Check | Criterion | Result |
|-------|-----------|--------|
| Sample quality | RIN > 7 | ✓ |
| Sequencing depth | >20M reads | ✓ |
| Mapping rate | >80% | ✓ |

# Methods

## Data Processing

Raw reads were trimmed using Trimmomatic and aligned with STAR.

```bash
STAR --runThreadN 8 --genomeDir $GENOME_DIR --readFilesIn $READS
```

# Results

## Differential Expression

#### Finding 1: Significant upregulation in treatment group

We identified 342 significantly upregulated genes in the treatment group.

#### Finding 2: Pathway enrichment confirms mechanism

KEGG analysis showed enrichment in oxidative phosphorylation.

> **Note:** All p-values are FDR-corrected using Benjamini-Hochberg.

# Conclusions

The treatment induces widespread transcriptional changes consistent with metabolic adaptation.
"""

    # Parse the document
    parser = TemplateParser()
    tree = parser.parse(sample_md)

    print(f"\nDocument Type: {tree.frontmatter.type}")
    print(f"Title: {tree.frontmatter.title}")
    print(f"Date: {tree.frontmatter.date}")

    print(f"\nSections ({len(tree.sections)} top-level):")
    for section in tree.sections:
        print(f"  - [{section.level}] {section.heading}")
        for sub in section.subsections:
            print(f"      - [{sub.level}] {sub.heading}")
            # Show detected elements
            if sub.elements:
                for elem in sub.elements:
                    print(f"          → {elem.type}")

    return tree, sample_md


def demo_validate_document(tree):
    """Demonstrate document validation."""
    print("\n" + "=" * 60)
    print("VALIDATING DOCUMENT")
    print("=" * 60)

    parser = TemplateParser("analysis-report")
    errors = parser.validate(tree)

    if not errors:
        print("\n✓ Document is valid!")
    else:
        print(f"\nValidation issues ({len(errors)}):")
        for err in errors:
            icon = "✗" if err.level == "error" else "⚠"
            print(f"  {icon} [{err.level}] {err.message}")
            if err.location:
                print(f"      at: {err.location}")

    return errors


def demo_generate_pdf(tree, output_path: str = "template_demo.pdf"):
    """Demonstrate generating a branded PDF from parsed content."""
    print("\n" + "=" * 60)
    print("GENERATING PDF")
    print("=" * 60)

    doc = SimpleDocTemplate(
        output_path,
        pagesize=letter,
        leftMargin=0.75 * inch,
        rightMargin=0.75 * inch,
        topMargin=0.75 * inch,
        bottomMargin=0.75 * inch,
    )

    # Styles
    styles = getSampleStyleSheet()
    h2_style = ParagraphStyle(
        "H2",
        parent=styles["Heading2"],
        fontSize=14,
        spaceBefore=16,
        spaceAfter=8,
        textColor=BRAND_COLORS.BRAND_BLUE,
    )
    body_style = ParagraphStyle(
        "Body",
        parent=styles["Normal"],
        fontSize=10,
        spaceAfter=8,
        textColor=BRAND_COLORS.DARK_GRAY,
    )

    story = []

    # Header from frontmatter
    story.append(MetadataHeader(
        doc_type=tree.frontmatter.type.replace("-", " ").title(),
        metadata={
            "Title": tree.frontmatter.title,
            "Date": tree.frontmatter.date,
            **{k: str(v) for k, v in tree.frontmatter.raw.items()
               if k not in ("type", "title", "date")}
        },
        columns=3,
    ))
    story.append(Spacer(1, 20))

    # Process each section
    def process_section(section, depth=0):
        """Recursively process sections into PDF elements."""
        elements = []

        # Section heading
        if depth == 0:
            elements.append(SectionDivider(section.heading, style="line"))
        else:
            elements.append(Paragraph(section.heading, h2_style))

        # Section content (simplified - just first paragraph)
        if section.content:
            # Skip element markup, get plain text
            lines = section.content.split("\n")
            text_lines = [l for l in lines if l.strip() and not l.startswith("|")
                         and not l.startswith("```") and not l.startswith(">")
                         and not l.startswith("####")]
            if text_lines:
                elements.append(Paragraph(text_lines[0], body_style))

        # Process detected elements
        for elem in section.elements:
            if elem.type == "table":
                # Parse and render table
                table_data = parse_markdown_table(elem.content)
                if table_data:
                    # Check if it's a status table
                    has_status = any("✓" in str(cell) or "✗" in str(cell)
                                    for row in table_data[1:] for cell in row)
                    if has_status:
                        elements.append(StatusTable(
                            headers=table_data[0],
                            rows=table_data[1:],
                        ))
                    else:
                        elements.append(make_simple_table(table_data))

            elif elem.type == "finding_card":
                num = elem.attributes.get("number", 1)
                # Extract title from content
                title = elem.content.split("\n")[0] if elem.content else f"Finding {num}"
                desc = "\n".join(elem.content.split("\n")[1:]).strip() if elem.content else ""
                elements.append(FindingCard(number=num, title=title, description=desc))

            elif elem.type == "callout":
                callout_type = elem.attributes.get("callout_type", "Note")
                elements.append(CalloutBox(
                    title=callout_type,
                    text=elem.content,
                    box_type="info" if callout_type.lower() == "note" else "warning",
                ))

        elements.append(Spacer(1, 8))

        # Recursively process subsections
        for sub in section.subsections:
            elements.extend(process_section(sub, depth + 1))

        return elements

    for section in tree.sections:
        story.extend(process_section(section))

    # Build PDF
    doc.build(story)
    print(f"\n✓ Generated: {output_path}")
    return output_path


def parse_markdown_table(table_str: str) -> list[list[str]]:
    """Parse a markdown table into rows."""
    lines = [l.strip() for l in table_str.strip().split("\n") if l.strip()]
    if len(lines) < 2:
        return []

    result = []
    for i, line in enumerate(lines):
        # Skip separator line
        if i == 1 and set(line.replace("|", "").replace("-", "").replace(":", "").strip()) == set():
            continue

        cells = [c.strip() for c in line.split("|")]
        cells = [c for c in cells if c]  # Remove empty strings from edges
        if cells:
            result.append(cells)

    return result


def make_simple_table(data: list[list[str]]):
    """Create a simple styled table."""
    if not data:
        return Spacer(1, 0)

    table = Table(data)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), BRAND_COLORS.LIGHT_GRAY),
        ("TEXTCOLOR", (0, 0), (-1, 0), BRAND_COLORS.DARK_GRAY),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("GRID", (0, 0), (-1, -1), 0.5, BRAND_COLORS.MEDIUM_GRAY),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
    ]))
    return table


def main():
    """Run all demonstrations."""
    print("\n" + "=" * 60)
    print("MARKDOWN-TO-PDF TEMPLATE SYSTEM DEMO")
    print("=" * 60)

    # 1. List templates
    templates = demo_list_templates()

    # 2. Parse document
    tree, _ = demo_parse_document()

    # 3. Validate
    errors = demo_validate_document(tree)

    # 4. Generate PDF
    output = demo_generate_pdf(tree)

    print("\n" + "=" * 60)
    print("DEMO COMPLETE")
    print("=" * 60)
    print(f"\nGenerated files:")
    print(f"  - {output}")
    print("\nFor component showcase, run:")
    print("  uv run python examples/component_demo.py")


if __name__ == "__main__":
    main()
