"""Oligon Reports - Professional PDF report generation with Oligon brand standards."""

from .brand_colors import BRAND_COLORS, ColorCycles
from .components import (
    CalloutBox,
    FindingCard,
    FigurePlaceholder,
    GradedTable,
    MetadataHeader,
    MethodBlock,
    MetricCard,
    MetricCardRow,
    SectionDivider,
    StatusTable,
    Timeline,
)
from .report_generator import ReportGenerator
from .templates import (
    DocumentTree,
    Element,
    FrontmatterData,
    Section,
    TemplateParser,
    ValidationError,
)

__version__ = "0.1.0"
__all__ = [
    # Brand
    "BRAND_COLORS",
    "ColorCycles",
    # Components - Core
    "ReportGenerator",
    "MetricCard",
    "MetricCardRow",
    "CalloutBox",
    "Timeline",
    "SectionDivider",
    "FigurePlaceholder",
    # Components - Phase 2 (Template Integration)
    "FindingCard",
    "StatusTable",
    "GradedTable",
    "MethodBlock",
    "MetadataHeader",
    # Templates
    "TemplateParser",
    "DocumentTree",
    "FrontmatterData",
    "Section",
    "Element",
    "ValidationError",
]
