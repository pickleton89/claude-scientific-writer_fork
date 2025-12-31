#!/usr/bin/env python3
"""Test oligon-brand skill integration.

This script validates that the oligon-brand skill is properly integrated
and all adapters/tokens are functional.

Usage:
    python scripts/test_oligon_brand_integration.py
    # or
    uv run python scripts/test_oligon_brand_integration.py
"""

import sys
from pathlib import Path

# Add skill to path
skill_path = Path(__file__).parent.parent / "skills" / "oligon-brand"
sys.path.insert(0, str(skill_path))


def test_tokens_loading():
    """Test direct token loading."""
    import json

    tokens_path = skill_path / "tokens" / "brand-tokens.json"
    assert tokens_path.exists(), f"Tokens file not found: {tokens_path}"

    with open(tokens_path) as f:
        tokens = json.load(f)

    # Verify required top-level keys
    assert "meta" in tokens, "Missing 'meta' key in tokens"
    assert "colors" in tokens, "Missing 'colors' key in tokens"
    assert "color_cycles" in tokens, "Missing 'color_cycles' key in tokens"

    # Verify meta
    assert tokens["meta"]["name"] == "Oligon Scientific Brand"
    assert "version" in tokens["meta"]

    # Verify primary color
    assert tokens["colors"]["primary"]["brand_blue"]["hex"] == "#2DB2E8"

    # Verify color cycles exist
    assert "treatment_control" in tokens["color_cycles"]
    assert "neutrals_only" in tokens["color_cycles"]
    assert "opposing" in tokens["color_cycles"]

    print("✅ brand-tokens.json: All tests passed")


def test_matplotlib_adapter():
    """Test matplotlib adapter imports and functions."""
    from adapters.matplotlib_adapter import (
        BRAND_COLORS,
        add_panel_label,
        get_color,
        get_cycle,
        get_divergent_cmap,
        set_brand_style,
    )

    # Test color access
    assert BRAND_COLORS["brand_blue"] == "#2DB2E8", "BRAND_COLORS['brand_blue'] incorrect"
    assert get_color("brand_blue") == "#2DB2E8", "get_color('brand_blue') incorrect"

    # Test cycles
    cycle = get_cycle("treatment_control")
    assert len(cycle) == 3, f"treatment_control cycle should have 3 colors, got {len(cycle)}"
    assert cycle[1] == "#2DB2E8", "Treatment (index 1) should be Brand Blue"

    # Test opposing cycle
    opposing = get_cycle("opposing")
    assert "#2DB2E8" in opposing, "opposing cycle should contain Brand Blue"
    assert "#E8622D" in opposing, "opposing cycle should contain Contrast Orange"

    # Test style application
    set_brand_style("treatment_control")

    import matplotlib.pyplot as plt

    assert plt.rcParams["axes.spines.top"] == False, "Top spine should be disabled"
    assert plt.rcParams["axes.spines.right"] == False, "Right spine should be disabled"
    assert plt.rcParams["figure.facecolor"] == "white", "Figure background should be white"

    # Test divergent colormap
    cmap = get_divergent_cmap()
    assert cmap is not None, "Divergent colormap should not be None"

    # Test panel label function exists and is callable
    assert callable(add_panel_label), "add_panel_label should be callable"

    print("✅ matplotlib_adapter: All tests passed")


def test_mplstyle():
    """Test mplstyle file exists and is valid."""
    mplstyle_path = skill_path / "assets" / "oligon_color_brand.mplstyle"
    assert mplstyle_path.exists(), f"mplstyle not found: {mplstyle_path}"

    # Try to load it
    import matplotlib.pyplot as plt

    plt.style.use(str(mplstyle_path))

    # Verify some key settings from the style
    assert plt.rcParams["axes.spines.top"] == False, "mplstyle should disable top spine"
    assert plt.rcParams["axes.spines.right"] == False, "mplstyle should disable right spine"

    print("✅ oligon_color_brand.mplstyle: Loaded successfully")


def test_reference_files():
    """Test that reference documentation exists."""
    brand_colors_full = skill_path / "references" / "brand-colors-full.md"
    assert brand_colors_full.exists(), f"Reference file not found: {brand_colors_full}"

    # Verify it has content
    content = brand_colors_full.read_text()
    assert len(content) > 1000, "brand-colors-full.md seems too short"
    assert "#2DB2E8" in content, "brand-colors-full.md should contain Brand Blue hex"

    print("✅ Reference files: All exist and have content")


def test_skill_md():
    """Test that SKILL.md exists and has required sections."""
    skill_md = skill_path / "SKILL.md"
    assert skill_md.exists(), f"SKILL.md not found: {skill_md}"

    content = skill_md.read_text()

    # Check for required XML sections
    required_sections = [
        "<overview>",
        "<when_to_use>",
        "<decision_framework>",
        "<workflow>",
        "<success_criteria>",
        "<cross_references>",
        "<references>",
    ]

    for section in required_sections:
        assert section in content, f"SKILL.md missing required section: {section}"

    # Check frontmatter
    assert "name: oligon-brand" in content, "SKILL.md missing name in frontmatter"
    assert "version:" in content, "SKILL.md missing version in frontmatter"
    assert "brand-type:" in content, "SKILL.md missing brand-type in frontmatter"

    print("✅ SKILL.md: All required sections present")


if __name__ == "__main__":
    print("=" * 60)
    print("Oligon-Brand Integration Tests")
    print("=" * 60)
    print()

    try:
        test_tokens_loading()
        test_matplotlib_adapter()
        test_mplstyle()
        test_reference_files()
        test_skill_md()

        print()
        print("=" * 60)
        print("✅ All oligon-brand integration tests passed!")
        print("=" * 60)
        sys.exit(0)

    except AssertionError as e:
        print()
        print("=" * 60)
        print(f"❌ Test failed: {e}")
        print("=" * 60)
        sys.exit(1)

    except ImportError as e:
        print()
        print("=" * 60)
        print(f"❌ Import error: {e}")
        print("   Make sure matplotlib is installed: uv add matplotlib")
        print("=" * 60)
        sys.exit(1)
