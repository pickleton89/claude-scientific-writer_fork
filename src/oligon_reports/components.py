"""
Reusable PDF Report Components

Unique design elements for professional scientific reports:
- MetricCard: Highlighted key numbers/statistics
- MetricCardRow: Row of metric cards with automatic spacing
- CalloutBox: Critical findings or alerts
- Timeline: Visual milestone progression
- SectionDivider: Brand-colored section separators
- FindingCard: Numbered finding with title and description
- StatusTable: Table with automatic pass/fail color coding
- GradedTable: Table with tier-based color bands for rankings
- MethodBlock: What/Why/How structured block for methods
- MetadataHeader: Document header with key-value metadata pairs
- FigurePlaceholder: Placeholder for figures with caption area
"""

from reportlab.lib.units import inch
from reportlab.platypus import Flowable

from .brand_colors import BRAND_COLORS, Typography


class MetricCard(Flowable):
    """
    A highlighted metric card for displaying key numbers.

    Creates a rounded box with brand accent color containing
    a large number and descriptive label.
    """

    def __init__(
        self,
        value: str,
        label: str,
        width: float = 1.5 * inch,
        height: float = 0.9 * inch,
        accent_color=None,
    ):
        super().__init__()
        self.value = value
        self.label = label
        self.box_width = width
        self.box_height = height
        self.accent_color = accent_color or BRAND_COLORS.BRAND_BLUE
        self.width = width
        self.height = height

    def draw(self):
        canvas = self.canv

        # Draw rounded rectangle background
        canvas.setFillColor(self.accent_color)
        canvas.roundRect(
            0,
            0,
            self.box_width,
            self.box_height,
            radius=6,
            fill=1,
            stroke=0,
        )

        # Draw value (large, white, centered)
        canvas.setFillColor(BRAND_COLORS.WHITE)
        canvas.setFont(Typography.BOLD_FONT, 22)
        value_y = self.box_height * 0.55
        canvas.drawCentredString(self.box_width / 2, value_y, self.value)

        # Draw label (smaller, white, centered)
        canvas.setFont(Typography.PRIMARY_FONT, 9)
        label_y = self.box_height * 0.2
        canvas.drawCentredString(self.box_width / 2, label_y, self.label)


class MetricCardRow(Flowable):
    """A row of metric cards with automatic spacing."""

    def __init__(self, metrics: list[tuple[str, str]], accent_colors: list = None):
        """
        Args:
            metrics: List of (value, label) tuples
            accent_colors: Optional list of colors for each card
        """
        super().__init__()
        self.metrics = metrics
        self.accent_colors = accent_colors or [BRAND_COLORS.BRAND_BLUE] * len(metrics)
        self.card_width = 1.5 * inch
        self.card_height = 0.9 * inch
        self.spacing = 0.25 * inch
        self.width = len(metrics) * self.card_width + (len(metrics) - 1) * self.spacing
        self.height = self.card_height

    def draw(self):
        x_offset = 0
        for i, (value, label) in enumerate(self.metrics):
            if i < len(self.accent_colors):
                color = self.accent_colors[i]
            else:
                color = BRAND_COLORS.BRAND_BLUE
            card = MetricCard(value, label, self.card_width, self.card_height, color)
            card.canv = self.canv
            self.canv.saveState()
            self.canv.translate(x_offset, 0)
            card.draw()
            self.canv.restoreState()
            x_offset += self.card_width + self.spacing


class CalloutBox(Flowable):
    """
    A callout box for highlighting critical findings or alerts.

    Features a colored left border and subtle background.
    """

    def __init__(
        self,
        text: str,
        title: str = None,
        width: float = 5.5 * inch,
        box_type: str = "info",
    ):
        """
        Args:
            text: Main callout text
            title: Optional title for the callout
            width: Box width
            box_type: 'info' (brand blue), 'warning' (orange), 'neutral' (gray)
        """
        super().__init__()
        self.text = text
        self.title = title
        self.box_width = width
        self.box_type = box_type
        self.width = width

        # Set colors based on type
        color_map = {
            "info": BRAND_COLORS.BRAND_BLUE,
            "warning": BRAND_COLORS.CONTRAST_ORANGE,
            "neutral": BRAND_COLORS.DARK_TEAL,
            "success": BRAND_COLORS.MEDIUM_BLUE,
        }
        self.accent_color = color_map.get(box_type, BRAND_COLORS.BRAND_BLUE)

        # Calculate height based on text
        self._calculate_height()

    def _calculate_height(self):
        """Estimate box height based on text length."""
        # Rough estimate: 15 chars per line, 14pt line height
        text_width = self.box_width - 0.5 * inch
        chars_per_line = int(text_width / 6)  # ~6 points per char at 10pt
        num_lines = max(1, len(self.text) // chars_per_line + 1)
        title_height = 18 if self.title else 0
        self.height = num_lines * 14 + title_height + 24  # padding

    def draw(self):
        canvas = self.canv

        # Draw background
        canvas.setFillColor(BRAND_COLORS.GRIDLINE)
        canvas.rect(0, 0, self.box_width, self.height, fill=1, stroke=0)

        # Draw left accent bar
        bar_width = 4
        canvas.setFillColor(self.accent_color)
        canvas.rect(0, 0, bar_width, self.height, fill=1, stroke=0)

        # Draw title if present
        y_pos = self.height - 16
        if self.title:
            canvas.setFillColor(self.accent_color)
            canvas.setFont(Typography.BOLD_FONT, 11)
            canvas.drawString(bar_width + 12, y_pos, self.title)
            y_pos -= 16

        # Draw text
        canvas.setFillColor(BRAND_COLORS.DARK_GRAY)
        canvas.setFont(Typography.PRIMARY_FONT, 10)

        # Simple text wrapping
        words = self.text.split()
        line = ""
        x_start = bar_width + 12
        max_width = self.box_width - x_start - 12

        for word in words:
            test_line = f"{line} {word}".strip()
            if canvas.stringWidth(test_line, Typography.PRIMARY_FONT, 10) < max_width:
                line = test_line
            else:
                canvas.drawString(x_start, y_pos, line)
                y_pos -= 14
                line = word
        if line:
            canvas.drawString(x_start, y_pos, line)


class Timeline(Flowable):
    """
    A horizontal timeline showing project milestones.

    Creates a visual progression with labeled milestones.
    """

    def __init__(
        self,
        milestones: list[dict],
        width: float = 6 * inch,
        height: float = 1.2 * inch,
    ):
        """
        Args:
            milestones: List of dicts with 'label', 'date', and optional 'completed' bool
            width: Timeline width
            height: Timeline height
        """
        super().__init__()
        self.milestones = milestones
        self.timeline_width = width
        self.timeline_height = height
        self.width = width
        self.height = height

    def draw(self):
        canvas = self.canv
        n = len(self.milestones)
        if n == 0:
            return

        # Timeline bar position
        bar_y = self.timeline_height * 0.5
        bar_height = 4

        # Draw background bar (full timeline)
        canvas.setFillColor(BRAND_COLORS.GRIDLINE)
        canvas.roundRect(
            0, bar_y - bar_height / 2,
            self.timeline_width, bar_height,
            radius=2, fill=1, stroke=0
        )

        # Calculate milestone positions
        padding = 0.5 * inch
        usable_width = self.timeline_width - 2 * padding
        spacing = usable_width / (n - 1) if n > 1 else 0

        for i, milestone in enumerate(self.milestones):
            x = padding + i * spacing
            completed = milestone.get("completed", False)
            is_current = milestone.get("current", False)

            # Draw progress bar up to completed milestones
            if completed and i > 0:
                prev_x = padding + (i - 1) * spacing
                canvas.setFillColor(BRAND_COLORS.BRAND_BLUE)
                canvas.roundRect(
                    prev_x, bar_y - bar_height / 2,
                    spacing, bar_height,
                    radius=2, fill=1, stroke=0
                )

            # Draw milestone marker
            marker_radius = 8 if is_current else 6
            if completed:
                canvas.setFillColor(BRAND_COLORS.BRAND_BLUE)
            elif is_current:
                canvas.setFillColor(BRAND_COLORS.CONTRAST_ORANGE)
            else:
                canvas.setFillColor(BRAND_COLORS.MUTED_GRAY)

            canvas.circle(x, bar_y, marker_radius, fill=1, stroke=0)

            # Inner circle for current milestone
            if is_current:
                canvas.setFillColor(BRAND_COLORS.WHITE)
                canvas.circle(x, bar_y, marker_radius - 3, fill=1, stroke=0)

            # Draw label above
            canvas.setFillColor(BRAND_COLORS.DARK_GRAY)
            font = Typography.BOLD_FONT if (completed or is_current) else Typography.PRIMARY_FONT
            canvas.setFont(font, 8)
            label = milestone.get("label", "")
            canvas.drawCentredString(x, bar_y + 18, label)

            # Draw date below
            canvas.setFillColor(BRAND_COLORS.MEDIUM_GRAY)
            canvas.setFont(Typography.PRIMARY_FONT, 7)
            date = milestone.get("date", "")
            canvas.drawCentredString(x, bar_y - 22, date)


class SectionDivider(Flowable):
    """
    A branded section divider with optional title.
    """

    def __init__(
        self,
        title: str = None,
        width: float = 6.5 * inch,
        style: str = "line",
    ):
        """
        Args:
            title: Optional section title
            width: Divider width
            style: 'line', 'gradient', or 'dots'
        """
        super().__init__()
        self.title = title
        self.divider_width = width
        self.style = style
        self.width = width
        self.height = 30 if title else 12

    def draw(self):
        canvas = self.canv

        if self.style == "line":
            # Simple line with brand accent
            y = self.height / 2
            canvas.setStrokeColor(BRAND_COLORS.BRAND_BLUE)
            canvas.setLineWidth(2)
            canvas.line(0, y, self.divider_width * 0.15, y)

            if self.title:
                canvas.setFillColor(BRAND_COLORS.DARK_GRAY)
                canvas.setFont(Typography.BOLD_FONT, 12)
                canvas.drawString(self.divider_width * 0.15 + 10, y - 4, self.title)

        elif self.style == "gradient":
            # Gradient fade line
            y = self.height / 2
            steps = 20
            step_width = self.divider_width / steps
            for i in range(steps):
                alpha = 1 - (i / steps)
                canvas.setStrokeColor(BRAND_COLORS.BRAND_BLUE)
                canvas.setLineWidth(2)
                canvas.setStrokeAlpha(alpha)
                canvas.line(i * step_width, y, (i + 1) * step_width, y)
            canvas.setStrokeAlpha(1)  # Reset

        elif self.style == "dots":
            # Dotted pattern
            y = self.height / 2
            dot_spacing = 8
            num_dots = int(self.divider_width * 0.3 / dot_spacing)
            canvas.setFillColor(BRAND_COLORS.BRAND_BLUE)
            for i in range(num_dots):
                canvas.circle(i * dot_spacing + 2, y, 2, fill=1, stroke=0)

            if self.title:
                canvas.setFillColor(BRAND_COLORS.DARK_GRAY)
                canvas.setFont(Typography.BOLD_FONT, 12)
                canvas.drawString(num_dots * dot_spacing + 15, y - 4, self.title)


class FindingCard(Flowable):
    """
    A highlighted finding card with numbered badge.

    Used to display key research findings with visual emphasis.
    Maps to `#### Finding N:` patterns in markdown templates.
    """

    def __init__(
        self,
        number: int | str,
        title: str,
        description: str = "",
        width: float = 5.5 * inch,
        accent_color=None,
    ):
        """
        Args:
            number: Finding number (displayed in badge)
            title: Finding title (bold, prominent)
            description: Finding description (wrapped text)
            width: Card width
            accent_color: Badge and accent color (default: brand blue)
        """
        super().__init__()
        self.number = str(number)
        self.title = title
        self.description = description
        self.card_width = width
        self.accent_color = accent_color or BRAND_COLORS.BRAND_BLUE
        self.width = width

        # Calculate height based on content
        self._calculate_height()

    def _calculate_height(self):
        """Estimate card height based on content."""
        badge_size = 24  # Must match badge_size in draw()
        title_height = 20
        padding = 24

        # Estimate description lines
        if self.description:
            text_width = self.card_width - badge_size - 36  # margins
            chars_per_line = int(text_width / 5.5)  # ~5.5 points per char at 10pt
            num_lines = max(1, len(self.description) // chars_per_line + 1)
            desc_height = num_lines * 14
        else:
            desc_height = 0

        self.height = max(badge_size + 16, title_height + desc_height + padding)

    def draw(self):
        canvas = self.canv

        # Draw background
        canvas.setFillColor(BRAND_COLORS.GRIDLINE)
        canvas.roundRect(0, 0, self.card_width, self.height, radius=4, fill=1, stroke=0)

        # Draw left accent bar
        bar_width = 4
        canvas.setFillColor(self.accent_color)
        canvas.roundRect(0, 0, bar_width, self.height, radius=2, fill=1, stroke=0)

        # Draw number badge
        badge_size = 24
        badge_x = 16
        badge_y = self.height - badge_size - 8
        canvas.setFillColor(self.accent_color)
        canvas.circle(
            badge_x + badge_size / 2, badge_y + badge_size / 2, badge_size / 2, fill=1, stroke=0
        )

        # Draw number in badge
        canvas.setFillColor(BRAND_COLORS.WHITE)
        canvas.setFont(Typography.BOLD_FONT, 12)
        canvas.drawCentredString(
            badge_x + badge_size / 2, badge_y + badge_size / 2 - 4, self.number
        )

        # Draw title
        title_x = badge_x + badge_size + 12
        title_y = badge_y + badge_size / 2 - 5
        canvas.setFillColor(BRAND_COLORS.DARK_GRAY)
        canvas.setFont(Typography.BOLD_FONT, 11)
        canvas.drawString(title_x, title_y, self.title)

        # Draw description with wrapping
        if self.description:
            canvas.setFillColor(BRAND_COLORS.MEDIUM_GRAY)
            canvas.setFont(Typography.PRIMARY_FONT, 10)

            desc_x = title_x
            desc_y = title_y - 18
            max_width = self.card_width - desc_x - 12

            words = self.description.split()
            line = ""
            for word in words:
                test_line = f"{line} {word}".strip()
                if canvas.stringWidth(test_line, Typography.PRIMARY_FONT, 10) < max_width:
                    line = test_line
                else:
                    canvas.drawString(desc_x, desc_y, line)
                    desc_y -= 14
                    line = word
            if line:
                canvas.drawString(desc_x, desc_y, line)


class StatusTable(Flowable):
    """
    A table with automatic pass/fail color coding.

    Detects status symbols (✓, ✗, ✅, ❌, PASS, FAIL, Yes, No) and
    applies appropriate background colors to those cells.
    """

    # Symbols that indicate pass/success
    PASS_SYMBOLS = {"✓", "✅", "PASS", "Pass", "pass", "YES", "Yes", "yes", "True", "true"}
    # Symbols that indicate fail/failure
    FAIL_SYMBOLS = {"✗", "✘", "❌", "FAIL", "Fail", "fail", "NO", "No", "no", "False", "false"}

    def __init__(
        self,
        headers: list[str],
        rows: list[list[str]],
        col_widths: list[float] = None,
        width: float = 5.5 * inch,
    ):
        """
        Args:
            headers: Column header strings
            rows: List of row data (each row is a list of cell strings)
            col_widths: Optional column widths (defaults to equal distribution)
            width: Total table width
        """
        super().__init__()
        self.headers = headers
        self.rows = rows
        self.table_width = width
        self.width = width

        # Calculate column widths
        num_cols = len(headers)
        if col_widths:
            self.col_widths = col_widths
        else:
            self.col_widths = [width / num_cols] * num_cols

        # Build the table
        self._build_table()

    def _detect_status(self, cell_value: str) -> str | None:
        """Detect if cell contains a pass/fail indicator."""
        cell_str = str(cell_value).strip()
        if cell_str in self.PASS_SYMBOLS:
            return "pass"
        elif cell_str in self.FAIL_SYMBOLS:
            return "fail"
        return None

    def _build_table(self):
        """Build the ReportLab table with styling."""
        from reportlab.lib.colors import HexColor
        from reportlab.platypus import Table, TableStyle

        # Combine headers and rows
        data = [self.headers] + self.rows

        # Create base table
        self.table = Table(data, colWidths=self.col_widths)

        # Base style
        style_commands = [
            # Header row
            ("BACKGROUND", (0, 0), (-1, 0), BRAND_COLORS.BRAND_BLUE),
            ("TEXTCOLOR", (0, 0), (-1, 0), BRAND_COLORS.WHITE),
            ("FONTNAME", (0, 0), (-1, 0), Typography.BOLD_FONT),
            ("FONTSIZE", (0, 0), (-1, 0), 10),
            ("ALIGN", (0, 0), (-1, 0), "CENTER"),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            # Body rows
            ("FONTNAME", (0, 1), (-1, -1), Typography.PRIMARY_FONT),
            ("FONTSIZE", (0, 1), (-1, -1), 9),
            ("TEXTCOLOR", (0, 1), (-1, -1), BRAND_COLORS.DARK_GRAY),
            # Grid
            ("GRID", (0, 0), (-1, -1), 0.5, BRAND_COLORS.GRIDLINE),
            ("LINEBELOW", (0, 0), (-1, 0), 1, BRAND_COLORS.BRAND_BLUE),
            # Padding
            ("TOPPADDING", (0, 0), (-1, -1), 6),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ("LEFTPADDING", (0, 0), (-1, -1), 8),
            ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ]

        # Add pass/fail cell coloring
        pass_bg = HexColor("#E8F5E9")  # Light green
        fail_bg = HexColor("#FFEBEE")  # Light red

        for row_idx, row in enumerate(self.rows):
            for col_idx, cell in enumerate(row):
                status = self._detect_status(cell)
                table_row = row_idx + 1  # Account for header row
                if status == "pass":
                    style_commands.append(
                        ("BACKGROUND", (col_idx, table_row), (col_idx, table_row), pass_bg)
                    )
                elif status == "fail":
                    style_commands.append(
                        ("BACKGROUND", (col_idx, table_row), (col_idx, table_row), fail_bg)
                    )

        self.table.setStyle(TableStyle(style_commands))

        # Calculate height
        self.table.wrapOn(None, self.table_width, 1000)
        self.height = self.table._height

    def draw(self):
        self.table.drawOn(self.canv, 0, 0)


class GradedTable(Flowable):
    """
    A table with tier-based color bands for rankings.

    Applies row background colors based on tier/grade values in a specified column.
    Supports tier patterns like "Tier 1", "A", "High", or numeric ranks.
    """

    def __init__(
        self,
        headers: list[str],
        rows: list[list[str]],
        grade_column: int | str,
        col_widths: list[float] = None,
        width: float = 5.5 * inch,
        tier_colors: dict = None,
    ):
        """
        Args:
            headers: Column header strings
            rows: List of row data
            grade_column: Column index (int) or header name (str) containing tier/grade
            col_widths: Optional column widths
            width: Total table width
            tier_colors: Optional custom tier-to-color mapping
        """
        super().__init__()
        self.headers = headers
        self.rows = rows
        self.table_width = width
        self.width = width

        # Default tier color mapping
        from reportlab.lib.colors import HexColor

        default_tier_colors = {
            # Tier 1 / A / High - Gold tint
            1: HexColor("#FFF8E1"),
            "1": HexColor("#FFF8E1"),
            "Tier 1": HexColor("#FFF8E1"),
            "A": HexColor("#FFF8E1"),
            "High": HexColor("#FFF8E1"),
            # Tier 2 / B / Medium - Silver/blue tint
            2: HexColor("#E3F2FD"),
            "2": HexColor("#E3F2FD"),
            "Tier 2": HexColor("#E3F2FD"),
            "B": HexColor("#E3F2FD"),
            "Medium": HexColor("#E3F2FD"),
            # Tier 3 / C / Low - Bronze/tan tint
            3: HexColor("#EFEBE9"),
            "3": HexColor("#EFEBE9"),
            "Tier 3": HexColor("#EFEBE9"),
            "C": HexColor("#EFEBE9"),
            "Low": HexColor("#EFEBE9"),
        }
        self.tier_colors = tier_colors or default_tier_colors

        # Resolve grade column to index
        if isinstance(grade_column, str):
            try:
                self.grade_col_idx = headers.index(grade_column)
            except ValueError:
                raise ValueError(f"Column '{grade_column}' not found in headers: {headers}")
        else:
            self.grade_col_idx = grade_column

        # Calculate column widths
        num_cols = len(headers)
        if col_widths:
            self.col_widths = col_widths
        else:
            self.col_widths = [width / num_cols] * num_cols

        self._build_table()

    def _get_tier_color(self, tier_value: str):
        """Get background color for a tier value."""
        tier_str = str(tier_value).strip()

        # Direct lookup
        if tier_str in self.tier_colors:
            return self.tier_colors[tier_str]

        # Try numeric conversion
        try:
            tier_num = int(tier_str)
            if tier_num in self.tier_colors:
                return self.tier_colors[tier_num]
        except ValueError:
            pass

        return None  # No special coloring

    def _build_table(self):
        """Build the ReportLab table with tier-based styling."""
        from reportlab.platypus import Table, TableStyle

        data = [self.headers] + self.rows
        self.table = Table(data, colWidths=self.col_widths)

        style_commands = [
            # Header row
            ("BACKGROUND", (0, 0), (-1, 0), BRAND_COLORS.BRAND_BLUE),
            ("TEXTCOLOR", (0, 0), (-1, 0), BRAND_COLORS.WHITE),
            ("FONTNAME", (0, 0), (-1, 0), Typography.BOLD_FONT),
            ("FONTSIZE", (0, 0), (-1, 0), 10),
            ("ALIGN", (0, 0), (-1, 0), "CENTER"),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            # Body rows
            ("FONTNAME", (0, 1), (-1, -1), Typography.PRIMARY_FONT),
            ("FONTSIZE", (0, 1), (-1, -1), 9),
            ("TEXTCOLOR", (0, 1), (-1, -1), BRAND_COLORS.DARK_GRAY),
            # Grid
            ("GRID", (0, 0), (-1, -1), 0.5, BRAND_COLORS.GRIDLINE),
            ("LINEBELOW", (0, 0), (-1, 0), 1, BRAND_COLORS.BRAND_BLUE),
            # Padding
            ("TOPPADDING", (0, 0), (-1, -1), 6),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ("LEFTPADDING", (0, 0), (-1, -1), 8),
            ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ]

        # Apply tier-based row coloring
        for row_idx, row in enumerate(self.rows):
            tier_value = row[self.grade_col_idx]
            color = self._get_tier_color(tier_value)
            if color:
                table_row = row_idx + 1
                style_commands.append(
                    ("BACKGROUND", (0, table_row), (-1, table_row), color)
                )

        self.table.setStyle(TableStyle(style_commands))
        self.table.wrapOn(None, self.table_width, 1000)
        self.height = self.table._height

    def draw(self):
        self.table.drawOn(self.canv, 0, 0)


class MethodBlock(Flowable):
    """
    A What/Why/How structured block for method documentation.

    Displays three distinct sections with colored accents:
    - What: Brief description (brand blue)
    - Why: Rationale (contrast orange)
    - How: Code/commands (dark teal, monospace)
    """

    def __init__(
        self,
        what: str,
        why: str,
        how: str,
        title: str = None,
        width: float = 5.5 * inch,
    ):
        """
        Args:
            what: Brief description of the procedure
            why: Rationale for this approach
            how: Code or commands (displayed in monospace)
            title: Optional title for the method block
            width: Block width
        """
        super().__init__()
        self.what = what
        self.why = why
        self.how = how
        self.title = title
        self.block_width = width
        self.width = width

        self._calculate_height()

    def _calculate_height(self):
        """Estimate block height based on content."""
        section_padding = 12
        line_height = 14
        code_line_height = 12

        # Title height
        title_height = 24 if self.title else 0

        # What section
        what_lines = max(1, len(self.what) // 60 + 1)
        what_height = 20 + what_lines * line_height + section_padding

        # Why section
        why_lines = max(1, len(self.why) // 60 + 1)
        why_height = 20 + why_lines * line_height + section_padding

        # How section (code block)
        how_lines = max(1, self.how.count("\n") + 1)
        how_height = 20 + how_lines * code_line_height + section_padding + 16

        self.height = title_height + what_height + why_height + how_height + 8

    def _draw_section(self, canvas, y_pos, label, text, accent_color, is_code=False):
        """Draw a single section with accent bar."""
        section_height = 20
        line_height = 14 if not is_code else 12

        # Calculate text lines
        if is_code:
            lines = text.split("\n")
        else:
            # Simple word wrap
            words = text.split()
            lines = []
            current_line = ""
            for word in words:
                test = f"{current_line} {word}".strip()
                if len(test) < 65:
                    current_line = test
                else:
                    if current_line:
                        lines.append(current_line)
                    current_line = word
            if current_line:
                lines.append(current_line)

        text_height = len(lines) * line_height
        total_height = section_height + text_height + 8

        # Draw accent bar
        bar_width = 3
        canvas.setFillColor(accent_color)
        canvas.rect(0, y_pos - total_height, bar_width, total_height, fill=1, stroke=0)

        # Draw label
        canvas.setFillColor(accent_color)
        canvas.setFont(Typography.BOLD_FONT, 10)
        canvas.drawString(12, y_pos - 14, label)

        # Draw text/code
        text_y = y_pos - 30
        if is_code:
            # Code block background
            code_padding = 8
            code_width = self.block_width - 24
            code_height = text_height + code_padding * 2
            canvas.setFillColor(BRAND_COLORS.GRIDLINE)
            canvas.roundRect(
                12, text_y - code_height + line_height, code_width, code_height, radius=3, fill=1, stroke=0
            )

            canvas.setFillColor(BRAND_COLORS.DARK_GRAY)
            canvas.setFont("Courier", 9)
            for line in lines:
                canvas.drawString(20, text_y, line)
                text_y -= line_height
        else:
            canvas.setFillColor(BRAND_COLORS.MEDIUM_GRAY)
            canvas.setFont(Typography.PRIMARY_FONT, 10)
            for line in lines:
                canvas.drawString(12, text_y, line)
                text_y -= line_height

        return total_height

    def draw(self):
        canvas = self.canv

        # Draw background
        canvas.setFillColor(BRAND_COLORS.WHITE)
        canvas.setStrokeColor(BRAND_COLORS.GRIDLINE)
        canvas.roundRect(0, 0, self.block_width, self.height, radius=4, fill=1, stroke=1)

        y_pos = self.height - 8

        # Draw title if present
        if self.title:
            canvas.setFillColor(BRAND_COLORS.DARK_GRAY)
            canvas.setFont(Typography.BOLD_FONT, 12)
            canvas.drawString(12, y_pos - 12, self.title)
            y_pos -= 28

        # Draw What section
        h = self._draw_section(canvas, y_pos, "What", self.what, BRAND_COLORS.BRAND_BLUE)
        y_pos -= h

        # Draw Why section
        h = self._draw_section(canvas, y_pos, "Why", self.why, BRAND_COLORS.CONTRAST_ORANGE)
        y_pos -= h

        # Draw How section
        self._draw_section(canvas, y_pos, "How", self.how, BRAND_COLORS.DARK_TEAL, is_code=True)


class MetadataHeader(Flowable):
    """
    A document header displaying key-value metadata pairs.

    Used to display YAML frontmatter data at the top of documents.
    Renders as a compact grid with document type badge.
    """

    def __init__(
        self,
        doc_type: str,
        metadata: dict[str, str],
        width: float = 5.5 * inch,
        columns: int = 2,
    ):
        """
        Args:
            doc_type: Document type label (e.g., "Analysis Report")
            metadata: Dictionary of key-value pairs to display
            width: Header width
            columns: Number of columns for metadata grid (default: 2)
        """
        super().__init__()
        self.doc_type = doc_type.upper()
        self.metadata = metadata
        self.header_width = width
        self.width = width
        self.columns = columns

        self._calculate_height()

    def _calculate_height(self):
        """Calculate header height based on metadata count."""
        type_badge_height = 28
        row_height = 20
        padding = 16

        # Calculate rows needed for metadata
        num_items = len(self.metadata)
        num_rows = (num_items + self.columns - 1) // self.columns

        self.height = type_badge_height + (num_rows * row_height) + padding

    def draw(self):
        canvas = self.canv

        # Draw outer border
        canvas.setStrokeColor(BRAND_COLORS.GRIDLINE)
        canvas.setLineWidth(1)
        canvas.roundRect(0, 0, self.header_width, self.height, radius=4, fill=0, stroke=1)

        # Draw document type badge
        badge_height = 24
        badge_y = self.height - badge_height - 2
        canvas.setFillColor(BRAND_COLORS.BRAND_BLUE)
        canvas.roundRect(0, badge_y, self.header_width, badge_height + 2, radius=4, fill=1, stroke=0)
        # Cover bottom corners of badge
        canvas.rect(0, badge_y, self.header_width, 4, fill=1, stroke=0)

        # Draw document type text
        canvas.setFillColor(BRAND_COLORS.WHITE)
        canvas.setFont(Typography.BOLD_FONT, 11)
        canvas.drawCentredString(self.header_width / 2, badge_y + 7, self.doc_type)

        # Draw metadata grid
        col_width = self.header_width / self.columns
        row_height = 20
        start_y = badge_y - 8

        items = list(self.metadata.items())
        for i, (key, value) in enumerate(items):
            col = i % self.columns
            row = i // self.columns

            x = col * col_width + 12
            y = start_y - (row * row_height)

            # Draw key
            canvas.setFillColor(BRAND_COLORS.MUTED_GRAY)
            canvas.setFont(Typography.PRIMARY_FONT, 8)
            canvas.drawString(x, y, f"{key}:")

            # Draw value
            canvas.setFillColor(BRAND_COLORS.DARK_GRAY)
            canvas.setFont(Typography.BOLD_FONT, 9)
            key_width = canvas.stringWidth(f"{key}: ", Typography.PRIMARY_FONT, 8)
            canvas.drawString(x + key_width, y, str(value))


class FigurePlaceholder(Flowable):
    """
    A placeholder for figures with brand-styled border and caption area.
    """

    def __init__(
        self,
        width: float = 5 * inch,
        height: float = 3 * inch,
        caption: str = None,
        figure_id: str = None,
    ):
        super().__init__()
        self.fig_width = width
        self.fig_height = height
        self.caption = caption
        self.figure_id = figure_id
        self.width = width
        caption_height = 20 if caption else 0
        self.height = height + caption_height

    def draw(self):
        canvas = self.canv

        # Draw placeholder box
        canvas.setStrokeColor(BRAND_COLORS.LIGHT_GRAY)
        canvas.setLineWidth(1)
        canvas.setDash([4, 4])
        caption_offset = 20 if self.caption else 0
        canvas.rect(0, caption_offset, self.fig_width, self.fig_height, fill=0, stroke=1)
        canvas.setDash([])

        # Draw "Figure" text in center
        canvas.setFillColor(BRAND_COLORS.MUTED_GRAY)
        canvas.setFont(Typography.PRIMARY_FONT, 12)
        label = f"Figure {self.figure_id}" if self.figure_id else "Figure"
        canvas.drawCentredString(
            self.fig_width / 2,
            caption_offset + self.fig_height / 2,
            label
        )

        # Draw caption below
        if self.caption:
            canvas.setFillColor(BRAND_COLORS.DARK_GRAY)
            canvas.setFont(Typography.PRIMARY_FONT, 9)
            if self.figure_id:
                caption_text = f"Figure {self.figure_id}. {self.caption}"
            else:
                caption_text = self.caption
            canvas.drawString(0, 4, caption_text)
