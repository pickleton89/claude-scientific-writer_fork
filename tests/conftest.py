"""Pytest configuration and fixtures for oligon_reports integration tests."""

from pathlib import Path

import pytest


@pytest.fixture
def fixtures_dir() -> Path:
    """Return path to test fixtures directory."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def sample_documents_dir(fixtures_dir: Path) -> Path:
    """Return path to sample documents directory."""
    return fixtures_dir / "sample_documents"


@pytest.fixture
def output_dir(tmp_path: Path) -> Path:
    """Return temporary output directory for generated PDFs."""
    output = tmp_path / "output"
    output.mkdir(exist_ok=True)
    return output


@pytest.fixture
def sample_analysis_report(sample_documents_dir: Path) -> str:
    """Load sample analysis report markdown."""
    return (sample_documents_dir / "analysis_report.md").read_text()


@pytest.fixture
def sample_meeting_notes(sample_documents_dir: Path) -> str:
    """Load sample meeting notes markdown."""
    return (sample_documents_dir / "meeting_notes.md").read_text()
