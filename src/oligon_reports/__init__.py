"""Oligon Reports - Professional PDF report generation with Oligon brand standards."""

from .brand_colors import BRAND_COLORS, ColorCycles
from .components import CalloutBox, MetricCard, SectionDivider, Timeline
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
    # Components
    "ReportGenerator",
    "MetricCard",
    "CalloutBox",
    "Timeline",
    "SectionDivider",
    # Templates
    "TemplateParser",
    "DocumentTree",
    "FrontmatterData",
    "Section",
    "Element",
    "ValidationError",
]
