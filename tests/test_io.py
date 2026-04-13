"""
Tests for mitoedit.io — sequence file reading and parsing utilities.

Covers:
  - .txt plain-text files (no header)
  - .fasta / .fa files (single and multi-line sequences, header ignored)
  - In-memory parsing via parse_sequence_content
  - Error cases: empty content, unsupported extension
"""

import os
import textwrap

import pytest

from mitoedit.io import parse_sequence_content, read_sequence_file

SEQ = "ATCGATCGATCG"
SEQ_UPPER = SEQ.upper()

# A minimal but realistic multi-line FASTA (two header lines would be unusual
# but we only have one here; the sequence is split across two lines)
FASTA_SINGLE = f">MT sample sequence\n{SEQ[:6]}\n{SEQ[6:]}\n"
FASTA_LOWER = f">MT lower\n{SEQ.lower()}\n"
TXT_CONTENT = SEQ + "\n"
TXT_CONTENT_LOWER = SEQ.lower() + "\n"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write(tmp_path, filename, content):
    p = tmp_path / filename
    p.write_text(content)
    return str(p)


# ---------------------------------------------------------------------------
# read_sequence_file — .txt
# ---------------------------------------------------------------------------

def test_txt_plain(tmp_path):
    path = _write(tmp_path, "seq.txt", TXT_CONTENT)
    assert read_sequence_file(path) == SEQ_UPPER


def test_txt_lower_is_uppercased(tmp_path):
    path = _write(tmp_path, "seq.txt", TXT_CONTENT_LOWER)
    assert read_sequence_file(path) == SEQ_UPPER


def test_txt_multiline_joined(tmp_path):
    content = f"{SEQ[:6]}\n{SEQ[6:]}\n"
    path = _write(tmp_path, "seq.txt", content)
    assert read_sequence_file(path) == SEQ_UPPER


# ---------------------------------------------------------------------------
# read_sequence_file — .fasta / .fa
# ---------------------------------------------------------------------------

def test_fasta_header_ignored(tmp_path):
    path = _write(tmp_path, "seq.fasta", FASTA_SINGLE)
    assert read_sequence_file(path) == SEQ_UPPER


def test_fa_extension(tmp_path):
    path = _write(tmp_path, "seq.fa", FASTA_SINGLE)
    assert read_sequence_file(path) == SEQ_UPPER


def test_fasta_lower_uppercased(tmp_path):
    path = _write(tmp_path, "seq.fasta", FASTA_LOWER)
    assert read_sequence_file(path) == SEQ_UPPER


def test_fasta_multiline_sequence(tmp_path):
    content = f">header\n{SEQ[:4]}\n{SEQ[4:8]}\n{SEQ[8:]}\n"
    path = _write(tmp_path, "seq.fasta", content)
    assert read_sequence_file(path) == SEQ_UPPER


def test_fasta_no_trailing_newline(tmp_path):
    content = f">header\n{SEQ}"
    path = _write(tmp_path, "seq.fasta", content)
    assert read_sequence_file(path) == SEQ_UPPER


# ---------------------------------------------------------------------------
# read_sequence_file — error cases
# ---------------------------------------------------------------------------

def test_unsupported_extension_raises(tmp_path):
    path = _write(tmp_path, "seq.fastq", SEQ)
    with pytest.raises(ValueError, match="Unsupported"):
        read_sequence_file(path)


def test_empty_txt_raises(tmp_path):
    path = _write(tmp_path, "seq.txt", "   \n  ")
    with pytest.raises(ValueError, match="No sequence"):
        read_sequence_file(path)


def test_empty_fasta_raises(tmp_path):
    path = _write(tmp_path, "seq.fasta", ">header only\n")
    with pytest.raises(ValueError, match="No sequence"):
        read_sequence_file(path)


# ---------------------------------------------------------------------------
# parse_sequence_content — in-memory variant
# ---------------------------------------------------------------------------

def test_parse_txt_content():
    assert parse_sequence_content(TXT_CONTENT, "seq.txt") == SEQ_UPPER


def test_parse_fasta_content():
    assert parse_sequence_content(FASTA_SINGLE, "seq.fasta") == SEQ_UPPER


def test_parse_fa_content():
    assert parse_sequence_content(FASTA_SINGLE, "seq.fa") == SEQ_UPPER


def test_parse_unsupported_extension_raises():
    with pytest.raises(ValueError, match="Unsupported"):
        parse_sequence_content(SEQ, "seq.gz")


def test_parse_empty_fasta_raises():
    with pytest.raises(ValueError, match="No sequence"):
        parse_sequence_content(">header\n\n", "seq.fasta")
