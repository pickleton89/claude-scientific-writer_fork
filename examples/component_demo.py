#!/usr/bin/env python3
"""
Component Demo - Visual showcase of all oligon_reports PDF components.

Generates a PDF demonstrating each component with realistic example data.
Run with: uv run python examples/component_demo.py
"""

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Spacer, Paragraph
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle

from oligon_reports import (
    BRAND_COLORS,
    # Core components
    MetricCard,
    MetricCardRow,
    CalloutBox,
    Timeline,
    SectionDivider,
    FigurePlaceholder,
    # Phase 2 components
    FindingCard,
    StatusTable,
    GradedTable,
    MethodBlock,
    MetadataHeader,
)


def build_demo_pdf(output_path: str = "component_demo.pdf"):
    """Generate a PDF showcasing all components."""

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
    title_style = ParagraphStyle(
        "DemoTitle",
        parent=styles["Heading1"],
        fontSize=24,
        spaceAfter=20,
        textColor=BRAND_COLORS.DARK_GRAY,
    )
    section_style = ParagraphStyle(
        "SectionTitle",
        parent=styles["Heading2"],
        fontSize=14,
        spaceBefore=20,
        spaceAfter=10,
        textColor=BRAND_COLORS.BRAND_BLUE,
    )
    desc_style = ParagraphStyle(
        "Description",
        parent=styles["Normal"],
        fontSize=9,
        textColor=BRAND_COLORS.MEDIUM_GRAY,
        spaceAfter=8,
    )

    story = []

    # Title
    story.append(Paragraph("Oligon Reports Component Demo", title_style))
    story.append(Paragraph(
        "Visual showcase of all PDF components available in the oligon_reports package. "
        "These components provide brand-consistent styling for scientific documents.",
        desc_style
    ))
    story.append(Spacer(1, 20))

    # --- MetadataHeader ---
    story.append(SectionDivider("MetadataHeader", style="line"))
    story.append(Paragraph(
        "Document header with type badge and key-value metadata grid. "
        "Maps to YAML frontmatter in templates.",
        desc_style
    ))
    story.append(MetadataHeader(
        doc_type="Analysis Report",
        metadata={
            "Project": "RNA-seq Analysis",
            "Date": "2025-12-28",
            "Author": "Research Team",
            "Version": "1.0",
            "Status": "Draft",
            "Classification": "Internal",
        },
        columns=3,
    ))
    story.append(Spacer(1, 20))

    # --- MetricCardRow ---
    story.append(SectionDivider("MetricCard & MetricCardRow", style="line"))
    story.append(Paragraph(
        "Highlighted metric cards for key statistics. Use for executive summaries.",
        desc_style
    ))
    story.append(MetricCardRow(
        metrics=[
            ("2,847", "Samples"),
            ("98.2%", "Coverage"),
            ("< 0.001", "P-value"),
            ("1.8x", "Fold Change"),
        ],
        accent_colors=[
            BRAND_COLORS.BRAND_BLUE,
            BRAND_COLORS.MEDIUM_BLUE,
            BRAND_COLORS.DARK_TEAL,
            BRAND_COLORS.CONTRAST_ORANGE,
        ]
    ))
    story.append(Spacer(1, 20))

    # --- FindingCard ---
    story.append(SectionDivider("FindingCard", style="line"))
    story.append(Paragraph(
        "Numbered finding cards with badge, title, and description. "
        "Maps to '#### Finding N:' patterns in markdown.",
        desc_style
    ))
    story.append(FindingCard(
        number=1,
        title="Significant upregulation in treatment group",
        description="Gene expression analysis revealed 342 significantly upregulated genes "
                   "(FDR < 0.05, log2FC > 1) in the treatment group compared to controls.",
    ))
    story.append(Spacer(1, 8))
    story.append(FindingCard(
        number=2,
        title="Pathway enrichment confirms mechanism",
        description="KEGG pathway analysis showed significant enrichment in oxidative "
                   "phosphorylation (p = 3.2e-12) and ribosome biogenesis (p = 1.8e-9).",
        accent_color=BRAND_COLORS.CONTRAST_ORANGE,
    ))
    story.append(Spacer(1, 20))

    # --- StatusTable ---
    story.append(SectionDivider("StatusTable", style="line"))
    story.append(Paragraph(
        "Table with automatic pass/fail color coding. "
        "Detects symbols like checkmarks, X marks, PASS/FAIL, Yes/No.",
        desc_style
    ))
    story.append(StatusTable(
        headers=["Check", "Criterion", "Result", "Notes"],
        rows=[
            ["1", "Sample quality (RIN > 7)", "✓", "Mean RIN: 8.4"],
            ["2", "Sequence depth (>20M reads)", "✓", "Avg: 32M reads"],
            ["3", "Mapping rate (>80%)", "✓", "Mean: 92.1%"],
            ["4", "Duplicate rate (<30%)", "✗", "Mean: 34.2%"],
            ["5", "rRNA contamination (<5%)", "PASS", "Mean: 2.1%"],
        ],
        col_widths=[0.5 * inch, 2 * inch, 0.8 * inch, 2 * inch],
    ))
    story.append(Spacer(1, 20))

    # --- GradedTable ---
    story.append(SectionDivider("GradedTable", style="line"))
    story.append(Paragraph(
        "Table with tier-based row coloring. Supports Tier 1/2/3, A/B/C, High/Medium/Low.",
        desc_style
    ))
    story.append(GradedTable(
        headers=["Rank", "Gene", "Log2FC", "Adj. P-value", "Tier"],
        rows=[
            ["1", "BRCA1", "3.42", "1.2e-15", "Tier 1"],
            ["2", "TP53", "2.98", "4.5e-12", "Tier 1"],
            ["3", "EGFR", "2.11", "2.3e-08", "Tier 2"],
            ["4", "MYC", "1.87", "8.1e-06", "Tier 2"],
            ["5", "KRAS", "1.45", "3.4e-04", "Tier 3"],
        ],
        grade_column="Tier",
        col_widths=[0.6 * inch, 1 * inch, 1 * inch, 1.2 * inch, 1 * inch],
    ))
    story.append(Spacer(1, 20))

    # --- MethodBlock ---
    story.append(SectionDivider("MethodBlock", style="line"))
    story.append(Paragraph(
        "What/Why/How structured block for method documentation. "
        "Ideal for protocols and technical specifications.",
        desc_style
    ))
    story.append(MethodBlock(
        title="Differential Expression Analysis",
        what="Compare gene expression between treatment and control groups using DESeq2.",
        why="DESeq2 provides robust normalization for RNA-seq count data and controls "
            "for library size differences while modeling gene-wise dispersion.",
        how="""library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, colData, ~condition)
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)""",
    ))
    story.append(Spacer(1, 20))

    # --- CalloutBox ---
    story.append(SectionDivider("CalloutBox", style="line"))
    story.append(Paragraph(
        "Callout boxes for alerts, notes, and warnings. "
        "Types: info (blue), warning (orange), neutral (teal), success (medium blue).",
        desc_style
    ))
    story.append(CalloutBox(
        title="Key Insight",
        text="The observed upregulation pattern is consistent with cellular stress response "
             "pathways, suggesting the treatment induces adaptive metabolic changes.",
        box_type="info",
    ))
    story.append(Spacer(1, 8))
    story.append(CalloutBox(
        title="Limitation",
        text="Batch effects were detected between sequencing runs. Results should be "
             "validated with an independent cohort before drawing definitive conclusions.",
        box_type="warning",
    ))
    story.append(Spacer(1, 20))

    # --- Timeline ---
    story.append(SectionDivider("Timeline", style="line"))
    story.append(Paragraph(
        "Horizontal timeline for project milestones. Supports completed, current, and pending states.",
        desc_style
    ))
    story.append(Timeline(
        milestones=[
            {"label": "Design", "date": "Q1 2025", "completed": True},
            {"label": "Data Collection", "date": "Q2 2025", "completed": True},
            {"label": "Analysis", "date": "Q3 2025", "current": True},
            {"label": "Validation", "date": "Q4 2025"},
            {"label": "Publication", "date": "Q1 2026"},
        ]
    ))
    story.append(Spacer(1, 20))

    # --- FigurePlaceholder ---
    story.append(SectionDivider("FigurePlaceholder", style="line"))
    story.append(Paragraph(
        "Placeholder for figures with caption area. Use for layout planning.",
        desc_style
    ))
    story.append(FigurePlaceholder(
        figure_id="1",
        caption="Volcano plot showing differential gene expression between treatment and control groups.",
        width=4 * inch,
        height=2.5 * inch,
    ))
    story.append(Spacer(1, 20))

    # --- SectionDivider styles ---
    story.append(Paragraph("SectionDivider Styles", section_style))
    story.append(Paragraph(
        "Three divider styles: line, gradient, and dots.",
        desc_style
    ))
    story.append(SectionDivider("Line Style", style="line"))
    story.append(Spacer(1, 8))
    story.append(SectionDivider("Gradient Style", style="gradient"))
    story.append(Spacer(1, 8))
    story.append(SectionDivider("Dots Style", style="dots"))

    # Build PDF
    doc.build(story)
    print(f"Generated: {output_path}")
    return output_path


if __name__ == "__main__":
    import sys
    output = sys.argv[1] if len(sys.argv) > 1 else "component_demo.pdf"
    build_demo_pdf(output)
