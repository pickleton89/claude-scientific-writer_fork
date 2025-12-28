"""Template parser for markdown documents with YAML frontmatter.

Parses markdown documents into structured DocumentTree representations,
validates against document type schemas, and provides template access.
"""

from __future__ import annotations

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

        # Parse sections into hierarchical structure
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

        # Validate sections against schema requirements
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
        """Parse markdown content into hierarchical sections.

        Builds a tree structure where deeper heading levels become subsections
        of their parent sections (e.g., h3 under h2, h4 under h3).

        Args:
            content: Raw markdown content

        Returns:
            List of top-level sections with nested subsections
        """
        # Remove frontmatter to get body content
        try:
            post = frontmatter.loads(content)
            body = post.content
        except Exception:
            body = content

        # Regex to match markdown headings
        heading_pattern = re.compile(r"^(#{1,6})\s+(.+)$", re.MULTILINE)

        matches = list(heading_pattern.finditer(body))
        if not matches:
            return []

        # Build flat list of sections with content
        flat_sections: list[Section] = []
        for i, match in enumerate(matches):
            level = len(match.group(1))
            heading = match.group(2).strip()

            # Get content until next heading or end
            start = match.end()
            end = matches[i + 1].start() if i + 1 < len(matches) else len(body)
            section_content = body[start:end].strip()

            # Generate section ID from heading
            section_id = self._heading_to_id(heading)

            # Detect elements within section
            elements = self._detect_elements(section_content)

            flat_sections.append(
                Section(
                    id=section_id,
                    heading=heading,
                    level=level,
                    content=section_content,
                    subsections=[],
                    elements=elements,
                )
            )

        # Build hierarchical structure
        return self._build_section_hierarchy(flat_sections)

    def _build_section_hierarchy(self, flat_sections: list[Section]) -> list[Section]:
        """Convert flat section list into nested hierarchy.

        Uses a stack-based approach to properly nest sections by heading level.

        Args:
            flat_sections: Flat list of all sections

        Returns:
            List of top-level sections with nested subsections
        """
        if not flat_sections:
            return []

        root_sections: list[Section] = []
        # Stack of (section, level) tuples for tracking parent context
        stack: list[tuple[Section, int]] = []

        for section in flat_sections:
            level = section.level

            # Pop stack until we find a parent with lower level
            while stack and stack[-1][1] >= level:
                stack.pop()

            if stack:
                # Add as subsection of the current parent
                parent = stack[-1][0]
                parent.subsections.append(section)
            else:
                # No parent - this is a root-level section
                root_sections.append(section)

            # Push current section onto stack
            stack.append((section, level))

        return root_sections

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

        Identifies tables, code blocks, callouts, lists, checklists,
        and finding cards for component mapping.

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

        # Detect code blocks (including directory trees)
        code_pattern = re.compile(r"```(\w*)\n(.*?)```", re.DOTALL)
        for match in code_pattern.finditer(content):
            code_content = match.group(2).strip()
            lang = match.group(1) or "text"

            # Check if it's a directory tree
            if self._is_directory_tree(code_content):
                elements.append(
                    Element(
                        type="directory_tree",
                        content=code_content,
                        attributes={"language": lang},
                    )
                )
            else:
                elements.append(
                    Element(
                        type="code_block",
                        content=code_content,
                        attributes={"language": lang},
                    )
                )

        # Detect callouts (blockquotes with bold type)
        callout_pattern = re.compile(
            r"^>\s+\*\*(\w+):\*\*\s*(.+?)(?=\n[^>]|\Z)", re.MULTILINE | re.DOTALL
        )
        for match in callout_pattern.finditer(content):
            elements.append(
                Element(
                    type="callout",
                    content=match.group(2).strip(),
                    attributes={"callout_type": match.group(1)},
                )
            )

        # Detect checklists (lines with [ ] or [x])
        checklist_pattern = re.compile(r"^[-*]\s+\[([ xX])\]\s+(.+)$", re.MULTILINE)
        checklist_matches = list(checklist_pattern.finditer(content))
        if checklist_matches:
            items = [
                {
                    "checked": m.group(1).lower() == "x",
                    "text": m.group(2).strip(),
                }
                for m in checklist_matches
            ]
            elements.append(
                Element(
                    type="checklist",
                    content="\n".join(m.group(0) for m in checklist_matches),
                    attributes={"items": items},
                )
            )

        # Detect numbered lists (only if no checklist found to avoid overlap)
        if not checklist_matches:
            numbered_pattern = re.compile(r"^\d+\.\s+(.+)$", re.MULTILINE)
            numbered_matches = list(numbered_pattern.finditer(content))
            if numbered_matches:
                items = [m.group(1).strip() for m in numbered_matches]
                elements.append(
                    Element(
                        type="numbered_list",
                        content="\n".join(m.group(0) for m in numbered_matches),
                        attributes={"items": items},
                    )
                )

        # Detect bullet lists (only if no checklist found)
        if not checklist_matches:
            bullet_pattern = re.compile(r"^[-*]\s+(?!\[[ xX]\])(.+)$", re.MULTILINE)
            bullet_matches = list(bullet_pattern.finditer(content))
            if bullet_matches:
                items = [m.group(1).strip() for m in bullet_matches]
                elements.append(
                    Element(
                        type="bullet_list",
                        content="\n".join(m.group(0) for m in bullet_matches),
                        attributes={"items": items},
                    )
                )

        # Detect finding cards (#### Finding N: pattern)
        finding_pattern = re.compile(
            r"####\s+Finding\s+(\d+):\s*(.+?)(?=####\s+Finding|\Z)",
            re.DOTALL | re.IGNORECASE,
        )
        for match in finding_pattern.finditer(content):
            elements.append(
                Element(
                    type="finding_card",
                    content=match.group(2).strip(),
                    attributes={"number": int(match.group(1))},
                )
            )

        return elements

    def _is_directory_tree(self, content: str) -> bool:
        """Check if content appears to be a directory tree structure.

        Looks for common tree characters (├, └, │) or multiple directory-like lines.
        Avoids false positives from URLs, imports, or code with slashes.

        Args:
            content: Code block content

        Returns:
            True if content looks like a directory tree
        """
        # Check for ASCII tree drawing characters (high confidence)
        tree_chars = ["├", "└", "│", "─"]
        if any(char in content for char in tree_chars):
            return True

        # Check for multiple lines ending with / (directory indicators)
        # Requires at least 2 such lines to avoid false positives
        dir_lines = re.findall(r"^\s*[\w.-]+/$", content, re.MULTILINE)
        if len(dir_lines) >= 2:
            return True

        # Check for indented path structure (lines starting with spaces + path)
        # Pattern: indented lines with folder/file patterns, no URLs
        indented_paths = re.findall(
            r"^[ ]{2,}[\w.-]+(?:/[\w.-]+)+/?$", content, re.MULTILINE
        )
        if len(indented_paths) >= 2:
            return True

        return False

    def _validate_frontmatter(self, fm: FrontmatterData) -> list[ValidationError]:
        """Validate frontmatter against schema requirements.

        Schema format expects:
          required:
            - name: fieldname
              expected: "value"  # optional
              type: string       # optional

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
                # New format: {name: "fieldname", expected: "value", ...}
                field_name = field_spec.get("name", "")
            else:
                # Simple format: just the field name as string
                field_name = field_spec

            if not field_name:
                continue

            # Check field exists and has a value
            if field_name not in fm.raw or not fm.raw[field_name]:
                errors.append(
                    ValidationError(
                        level="error",
                        message=f"Missing required frontmatter field: {field_name}",
                        location=f"frontmatter.{field_name}",
                    )
                )
                continue

            # Check expected value if specified
            if isinstance(field_spec, dict):
                expected = field_spec.get("expected")
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

        Checks required sections exist and have expected elements.

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

        # Flatten sections for easier lookup
        all_sections = self._flatten_sections(sections)

        if len(all_sections) < min_sections:
            errors.append(
                ValidationError(
                    level="error",
                    message=f"Document has {len(all_sections)} sections, minimum required: {min_sections}",
                    location="sections",
                )
            )

        # Check required sections from schema
        schema_sections = self._schema.get("sections", [])
        section_ids = {s.id for s in all_sections}
        section_headings = {self._heading_to_id(s.heading) for s in all_sections}

        for schema_section in schema_sections:
            section_id = schema_section.get("id", "")
            required = schema_section.get("required", False)
            heading = schema_section.get("heading", "")

            if required:
                # Match by id or normalized heading
                heading_id = self._heading_to_id(heading)
                if section_id not in section_ids and heading_id not in section_headings:
                    errors.append(
                        ValidationError(
                            level="error",
                            message=f"Missing required section: {heading or section_id}",
                            location=f"section.{section_id}",
                        )
                    )

            # Check expected elements if section exists
            expected_elements = schema_section.get("expected_elements", [])
            if expected_elements:
                matching_section = self._find_section(all_sections, section_id, heading)
                if matching_section:
                    errors.extend(
                        self._validate_section_elements(
                            matching_section, expected_elements, section_id
                        )
                    )

            # Recursively check subsections
            schema_subsections = schema_section.get("subsections", [])
            for sub in schema_subsections:
                sub_id = sub.get("id", "")
                sub_required = sub.get("required", False)
                sub_heading = sub.get("heading", "")

                if sub_required:
                    sub_heading_id = self._heading_to_id(sub_heading)
                    if sub_id not in section_ids and sub_heading_id not in section_headings:
                        errors.append(
                            ValidationError(
                                level="warning",
                                message=f"Missing expected subsection: {sub_heading or sub_id}",
                                location=f"section.{section_id}.{sub_id}",
                            )
                        )

        return errors

    def _flatten_sections(self, sections: list[Section]) -> list[Section]:
        """Flatten hierarchical sections into a single list.

        Args:
            sections: Hierarchical section list

        Returns:
            Flat list of all sections including subsections
        """
        result: list[Section] = []
        for section in sections:
            result.append(section)
            if section.subsections:
                result.extend(self._flatten_sections(section.subsections))
        return result

    def _find_section(
        self, sections: list[Section], section_id: str, heading: str
    ) -> Section | None:
        """Find a section by id or heading.

        Args:
            sections: List of sections to search
            section_id: Section identifier
            heading: Section heading text

        Returns:
            Matching section or None
        """
        heading_id = self._heading_to_id(heading)
        for section in sections:
            if section.id == section_id or section.id == heading_id:
                return section
        return None

    def _validate_section_elements(
        self,
        section: Section,
        expected: list[str],
        section_id: str,
    ) -> list[ValidationError]:
        """Validate that a section contains expected element types.

        Args:
            section: Section to validate
            expected: List of expected element types (e.g., ["tables", "code_blocks"])
            section_id: Section identifier for error reporting

        Returns:
            List of validation errors/warnings
        """
        errors: list[ValidationError] = []
        element_types = {e.type for e in section.elements}

        # Map schema element names to detected element types
        type_mapping = {
            "tables": "table",
            "code_blocks": "code_block",
            "checklists": "checklist",
            "method_blocks": "code_block",  # Method blocks are code blocks
            "finding_cards": "finding_card",
            "callouts": "callout",
            "numbered_lists": "numbered_list",
            "bullet_lists": "bullet_list",
        }

        for expected_type in expected:
            mapped_type = type_mapping.get(expected_type, expected_type)
            if mapped_type not in element_types:
                errors.append(
                    ValidationError(
                        level="warning",
                        message=f"Section '{section_id}' expects {expected_type} but none found",
                        location=f"section.{section_id}.elements",
                    )
                )

        return errors
