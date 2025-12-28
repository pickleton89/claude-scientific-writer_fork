"""Template infrastructure for markdown-to-PDF conversion.

This module provides:
- TemplateParser: Parse markdown with YAML frontmatter into structured DocumentTree
- Document type schemas for validation
- Markdown templates for each document type

Example:
    from oligon_reports.templates import TemplateParser, DocumentTree

    parser = TemplateParser("analysis-report")
    tree = parser.parse(markdown_content)
    errors = parser.validate(tree)
"""

from dataclasses import dataclass, field
from typing import Any

__all__ = [
    "FrontmatterData",
    "Element",
    "Section",
    "DocumentTree",
    "ValidationError",
    "TemplateParser",
]


@dataclass
class FrontmatterData:
    """Parsed YAML frontmatter from a document.

    Attributes:
        type: Document type identifier (e.g., "analysis-report")
        title: Document title
        date: Document date as string (YYYY-MM-DD format)
        raw: Complete frontmatter dictionary for additional fields
    """

    type: str
    title: str
    date: str
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass
class Element:
    """A detected element within a section (table, callout, code block, etc.).

    Attributes:
        type: Element type ("table", "callout", "code_block", "finding_card", etc.)
        content: Raw content of the element
        attributes: Additional parsed attributes (e.g., language for code blocks)
    """

    type: str
    content: str
    attributes: dict[str, Any] = field(default_factory=dict)


@dataclass
class Section:
    """A document section with heading, content, and nested structure.

    Attributes:
        id: Section identifier (derived from heading, e.g., "objective")
        heading: Full heading text
        level: Heading level (1-6)
        content: Raw markdown content of the section
        subsections: Nested child sections
        elements: Detected elements within this section
    """

    id: str
    heading: str
    level: int
    content: str
    subsections: list["Section"] = field(default_factory=list)
    elements: list[Element] = field(default_factory=list)


@dataclass
class DocumentTree:
    """Complete parsed representation of a markdown document.

    Attributes:
        frontmatter: Parsed YAML frontmatter
        sections: Top-level sections in document order
        raw_content: Original markdown content
    """

    frontmatter: FrontmatterData
    sections: list[Section]
    raw_content: str


@dataclass
class ValidationError:
    """A validation error or warning from schema validation.

    Attributes:
        level: Severity ("error" or "warning")
        message: Human-readable description of the issue
        location: Optional location hint (e.g., "frontmatter.date", "section:methods")
    """

    level: str  # "error" | "warning"
    message: str
    location: str | None = None


# Import TemplateParser after dataclasses are defined to avoid circular imports
from oligon_reports.templates.parser import TemplateParser
