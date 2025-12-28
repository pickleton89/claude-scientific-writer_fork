"""Template parser for markdown documents with YAML frontmatter.

Parses markdown documents into structured DocumentTree representations,
validates against document type schemas, and provides template access.
"""

from __future__ import annotations

import importlib.resources
import re
from pathlib import Path
from typing import Any

import frontmatter
import yaml

from oligon_reports.templates import (
    DocumentTree,
    Element,
    FrontmatterData,
    Section,
    ValidationError,
)

# Valid document types (schemas must exist for each)
DOCUMENT_TYPES = [
    "analysis-report",
    "phase-plan",
    "data-report",
    "literature-review",
    "meeting-notes",
    "project-status",
    "technical-spec",
    "task-list",
    "standards-guide",
    "agent-definition",
    "readme",
    "method-guide",
]


class TemplateParser:
    """Parser for markdown documents with YAML frontmatter.

    Parses markdown content into a structured DocumentTree, validates
    against document type schemas, and provides access to templates.

    Args:
        template_type: Optional document type to use for validation.
            If not provided, will be detected from frontmatter.

    Example:
        parser = TemplateParser("analysis-report")
        tree = parser.parse(markdown_content)
        errors = parser.validate(tree)

        # Or auto-detect type:
        parser = TemplateParser()
        doc_type = parser.detect_type(content)
    """

    def __init__(self, template_type: str | None = None) -> None:
        self.template_type = template_type
        self._schema: dict[str, Any] | None = None

        if template_type is not None:
            self._schema = self._load_schema(template_type)

    def parse(self, content: str) -> DocumentTree:
        """Parse markdown content into a DocumentTree.

        Args:
            content: Markdown string with YAML frontmatter

        Returns:
            DocumentTree with parsed frontmatter, sections, and elements

        Raises:
            ValueError: If frontmatter is missing or invalid
        """
        # Parse frontmatter
        fm_data = self._parse_frontmatter(content)

        # If no template type set, detect from frontmatter
        if self.template_type is None and fm_data.type:
            self.template_type = fm_data.type
            self._schema = self._load_schema(fm_data.type)

        # Parse sections (skeleton - will be completed in Step 7)
        sections = self._parse_sections(content)

        return DocumentTree(
            frontmatter=fm_data,
            sections=sections,
            raw_content=content,
        )

    def validate(self, tree: DocumentTree) -> list[ValidationError]:
        """Validate a DocumentTree against the schema.

        Args:
            tree: Parsed document tree to validate

        Returns:
            List of validation errors/warnings (empty if valid)
        """
        errors: list[ValidationError] = []

        if self._schema is None:
            errors.append(
                ValidationError(
                    level="warning",
                    message="No schema loaded, skipping validation",
                )
            )
            return errors

        # Validate frontmatter required fields
        errors.extend(self._validate_frontmatter(tree.frontmatter))

        # Validate sections (skeleton - will be completed in Step 7)
        errors.extend(self._validate_sections(tree.sections))

        return errors

    def detect_type(self, content: str) -> str | None:
        """Detect document type from content frontmatter.

        Args:
            content: Markdown string with YAML frontmatter

        Returns:
            Document type string if found, None otherwise
        """
        try:
            post = frontmatter.loads(content)
            doc_type = post.metadata.get("type")
            if doc_type in DOCUMENT_TYPES:
                return doc_type
            return None
        except Exception:
            return None

    @staticmethod
    def list_templates() -> list[dict[str, str]]:
        """List all available document templates.

        Returns:
            List of dicts with 'type', 'category', and 'description' keys
        """
        templates = []
        schemas_path = Path(__file__).parent / "schemas"

        if not schemas_path.exists():
            return templates

        for schema_file in sorted(schemas_path.glob("*.yaml")):
            try:
                with open(schema_file) as f:
                    schema = yaml.safe_load(f)
                    templates.append(
                        {
                            "type": schema.get("document_type", schema_file.stem),
                            "category": schema.get("category", "unknown"),
                            "description": schema.get("description", "").split("\n")[0],
                        }
                    )
            except Exception:
                continue

        return templates

    @staticmethod
    def get_template(template_type: str) -> str:
        """Get the markdown template for a document type.

        Args:
            template_type: Document type identifier

        Returns:
            Template markdown content

        Raises:
            FileNotFoundError: If template doesn't exist
        """
        template_path = Path(__file__).parent / "markdown" / f"{template_type}.md"

        if not template_path.exists():
            raise FileNotFoundError(f"Template not found: {template_type}")

        return template_path.read_text()

    # -------------------------------------------------------------------------
    # Private methods
    # -------------------------------------------------------------------------

    def _parse_frontmatter(self, content: str) -> FrontmatterData:
        """Parse YAML frontmatter from markdown content.

        Args:
            content: Markdown string with YAML frontmatter

        Returns:
            FrontmatterData with parsed fields

        Raises:
            ValueError: If frontmatter is missing required fields
        """
        try:
            post = frontmatter.loads(content)
        except Exception as e:
            raise ValueError(f"Failed to parse frontmatter: {e}") from e

        metadata = post.metadata

        # Extract required fields with defaults
        doc_type = metadata.get("type", "")
        title = metadata.get("title", "Untitled")
        date = metadata.get("date", "")

        # Convert date to string if it's a date object
        if hasattr(date, "isoformat"):
            date = date.isoformat()
        else:
            date = str(date) if date else ""

        return FrontmatterData(
            type=doc_type,
            title=title,
            date=date,
            raw=dict(metadata),
        )

    def _load_schema(self, template_type: str) -> dict[str, Any] | None:
        """Load a document type schema from the schemas directory.

        Args:
            template_type: Document type identifier

        Returns:
            Parsed schema dict, or None if not found
        """
        schema_path = Path(__file__).parent / "schemas" / f"{template_type}.yaml"

        if not schema_path.exists():
            return None

        try:
            with open(schema_path) as f:
                return yaml.safe_load(f)
        except Exception:
            return None

    def _parse_sections(self, content: str) -> list[Section]:
        """Parse markdown content into sections.

        Skeleton implementation - will be completed in Step 7.

        Args:
            content: Raw markdown content

        Returns:
            List of top-level sections
        """
        # Remove frontmatter to get body content
        try:
            post = frontmatter.loads(content)
            body = post.content
        except Exception:
            body = content

        sections: list[Section] = []
        # Regex to match markdown headings
        heading_pattern = re.compile(r"^(#{1,6})\s+(.+)$", re.MULTILINE)

        matches = list(heading_pattern.finditer(body))
        if not matches:
            return sections

        for i, match in enumerate(matches):
            level = len(match.group(1))
            heading = match.group(2).strip()

            # Get content until next heading or end
            start = match.end()
            end = matches[i + 1].start() if i + 1 < len(matches) else len(body)
            section_content = body[start:end].strip()

            # Generate section ID from heading
            section_id = self._heading_to_id(heading)

            # Detect elements within section (skeleton)
            elements = self._detect_elements(section_content)

            sections.append(
                Section(
                    id=section_id,
                    heading=heading,
                    level=level,
                    content=section_content,
                    subsections=[],  # Will be populated with proper nesting in Step 7
                    elements=elements,
                )
            )

        return sections

    def _heading_to_id(self, heading: str) -> str:
        """Convert a heading to a section ID.

        Examples:
            "1. Objective" -> "objective"
            "Key Questions" -> "key_questions"
            "### Results" -> "results"
        """
        # Remove leading numbers and punctuation
        cleaned = re.sub(r"^\d+\.\s*", "", heading)
        # Convert to lowercase and replace spaces/special chars with underscores
        section_id = re.sub(r"[^a-z0-9]+", "_", cleaned.lower())
        # Remove leading/trailing underscores
        return section_id.strip("_")

    def _detect_elements(self, content: str) -> list[Element]:
        """Detect special elements within section content.

        Skeleton implementation - will be completed in Step 7.

        Args:
            content: Section content

        Returns:
            List of detected elements
        """
        elements: list[Element] = []

        # Detect tables
        if "|" in content and "---" in content:
            table_pattern = re.compile(
                r"(\|.+\|\n\|[-:| ]+\|\n(?:\|.+\|\n?)+)", re.MULTILINE
            )
            for match in table_pattern.finditer(content):
                elements.append(
                    Element(type="table", content=match.group(1).strip())
                )

        # Detect code blocks
        code_pattern = re.compile(r"```(\w*)\n(.*?)```", re.DOTALL)
        for match in code_pattern.finditer(content):
            elements.append(
                Element(
                    type="code_block",
                    content=match.group(2).strip(),
                    attributes={"language": match.group(1) or "text"},
                )
            )

        # Detect callouts (blockquotes with bold type)
        callout_pattern = re.compile(r"^>\s+\*\*(\w+):\*\*\s*(.+?)(?=\n[^>]|\Z)", re.MULTILINE | re.DOTALL)
        for match in callout_pattern.finditer(content):
            elements.append(
                Element(
                    type="callout",
                    content=match.group(2).strip(),
                    attributes={"callout_type": match.group(1)},
                )
            )

        return elements

    def _validate_frontmatter(self, fm: FrontmatterData) -> list[ValidationError]:
        """Validate frontmatter against schema requirements.

        Args:
            fm: Parsed frontmatter data

        Returns:
            List of validation errors
        """
        errors: list[ValidationError] = []

        if self._schema is None:
            return errors

        fm_schema = self._schema.get("frontmatter", {})
        required_fields = fm_schema.get("required", [])

        # Check required fields
        for field_spec in required_fields:
            if isinstance(field_spec, dict):
                field_name = list(field_spec.keys())[0]
            else:
                field_name = field_spec

            if field_name not in fm.raw or not fm.raw[field_name]:
                errors.append(
                    ValidationError(
                        level="error",
                        message=f"Missing required frontmatter field: {field_name}",
                        location=f"frontmatter.{field_name}",
                    )
                )

        # Check type matches if schema specifies expected type
        for field_spec in required_fields:
            if isinstance(field_spec, dict):
                field_name = list(field_spec.keys())[0]
                expected = field_spec[field_name].get("expected")
                if expected and fm.raw.get(field_name) != expected:
                    errors.append(
                        ValidationError(
                            level="error",
                            message=f"Field '{field_name}' expected '{expected}', got '{fm.raw.get(field_name)}'",
                            location=f"frontmatter.{field_name}",
                        )
                    )

        return errors

    def _validate_sections(self, sections: list[Section]) -> list[ValidationError]:
        """Validate sections against schema requirements.

        Skeleton implementation - will be completed in Step 7.

        Args:
            sections: List of parsed sections

        Returns:
            List of validation errors
        """
        errors: list[ValidationError] = []

        if self._schema is None:
            return errors

        validation_rules = self._schema.get("validation", {})
        min_sections = validation_rules.get("min_sections", 0)

        if len(sections) < min_sections:
            errors.append(
                ValidationError(
                    level="error",
                    message=f"Document has {len(sections)} sections, minimum required: {min_sections}",
                    location="sections",
                )
            )

        # More detailed section validation will be added in Step 7

        return errors
