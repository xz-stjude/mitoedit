"""
io.py — Sequence file reading utilities for MitoEdit.

Supported formats:
  .txt          Plain text — the entire file content is the sequence (no header).
  .fasta / .fa  FASTA format — one or more header lines starting with '>' are
                ignored; all remaining lines are joined into the sequence.
"""


def read_sequence_file(path: str) -> str:
    """Read a DNA sequence from *path* and return it as a clean uppercase string.

    Format is inferred from the file extension:
    - ``.txt``           → plain text, no header expected.
    - ``.fasta`` / ``.fa`` → FASTA; header lines (starting with ``>``) are skipped.

    Parameters
    ----------
    path : str
        Absolute or relative path to the sequence file.

    Returns
    -------
    str
        The DNA sequence with all whitespace removed and letters uppercased.

    Raises
    ------
    ValueError
        If the file extension is not recognised or the parsed sequence is empty.
    """
    with open(path, "r") as fh:
        content = fh.read()

    lower = path.lower()
    if lower.endswith(".fasta") or lower.endswith(".fa"):
        seq = _parse_fasta(content)
    elif lower.endswith(".txt"):
        seq = _parse_plain(content)
    else:
        raise ValueError(
            f"Unsupported sequence file extension for '{path}'. "
            "Use .txt for plain text or .fasta / .fa for FASTA format."
        )

    if not seq:
        raise ValueError(f"No sequence found in '{path}'.")
    return seq


def parse_sequence_content(content: str, filename: str) -> str:
    """Parse a DNA sequence from *content* using *filename* to detect format.

    This is the in-memory equivalent of :func:`read_sequence_file`, intended
    for use when the file content has already been read (e.g. from an uploaded
    web form).

    Parameters
    ----------
    content : str
        Raw file content as a string.
    filename : str
        Original filename (used only for extension detection).

    Returns
    -------
    str
        The DNA sequence with all whitespace removed and letters uppercased.

    Raises
    ------
    ValueError
        If the extension is not recognised or the parsed sequence is empty.
    """
    lower = filename.lower()
    if lower.endswith(".fasta") or lower.endswith(".fa"):
        seq = _parse_fasta(content)
    elif lower.endswith(".txt"):
        seq = _parse_plain(content)
    else:
        raise ValueError(
            f"Unsupported sequence file extension '{filename}'. "
            "Use .txt for plain text or .fasta / .fa for FASTA format."
        )

    if not seq:
        raise ValueError(f"No sequence found in '{filename}'.")
    return seq


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_plain(content: str) -> str:
    """Strip whitespace from plain-text sequence content."""
    return "".join(content.split()).upper()


def _parse_fasta(content: str) -> str:
    """Join all non-header lines from FASTA content into a single sequence."""
    lines = []
    for line in content.splitlines():
        stripped = line.strip()
        if stripped.startswith(">") or not stripped:
            continue
        lines.append(stripped.upper())
    return "".join(lines)
