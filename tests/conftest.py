"""
Shared pytest fixtures for MitoEdit tests.
"""

import pytest
import tempfile
import os
from pathlib import Path


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_fasta_content():
    """Sample FASTA content for testing."""
    return """>test_sequence_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>test_sequence_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
"""


@pytest.fixture
def sample_fasta_file(temp_dir, sample_fasta_content):
    """Create a temporary FASTA file with sample content."""
    fasta_file = temp_dir / "test.fasta"
    fasta_file.write_text(sample_fasta_content)
    return str(fasta_file)


@pytest.fixture
def empty_fasta_file(temp_dir):
    """Create an empty FASTA file."""
    fasta_file = temp_dir / "empty.fasta"
    fasta_file.write_text(">empty_sequence\n\n")
    return str(fasta_file)
