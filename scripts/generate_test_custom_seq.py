"""
generate_test_custom_seq.py — Generate custom-sequence test fixture files.

Produces a version of the human reference mtDNA with 200 bases deleted from
positions 1001-1200 (1-indexed) and writes it to tests/input/ in two formats:

  tests/input/custom_seq_200bp_deletion.txt   — plain text, no header
  tests/input/custom_seq_200bp_deletion.fa    — single-entry FASTA

These files are used by tests/test_alignment_liftover.py.

Usage (from the repo root):
    python scripts/generate_test_custom_seq.py
"""

import sys
from importlib.resources import files
from pathlib import Path

# Ensure the repo root is on sys.path when the script is run directly
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))


# ---------------------------------------------------------------------------
# Deletion parameters — must stay in sync with test_alignment_liftover.py
# ---------------------------------------------------------------------------
DELETION_START = 1000   # 0-indexed, inclusive
DELETION_END = 1200     # 0-indexed, exclusive  (removes exactly 200 bases)

FASTA_HEADER = (
    f">custom_mtDNA_200bp_deletion_pos{DELETION_START + 1}-{DELETION_END} "
    f"[human mtDNA NC_012920.1 with bases {DELETION_START + 1}-{DELETION_END} deleted]"
)


def generate(output_dir: Path | None = None) -> tuple[Path, Path]:
    """Build the deletion sequence and write both output files.

    Parameters
    ----------
    output_dir : Path, optional
        Directory to write the files into.  Defaults to ``tests/input/``
        relative to the repository root (i.e. two levels above this script).

    Returns
    -------
    txt_path, fa_path : tuple[Path, Path]
        Paths to the written files.
    """
    if output_dir is None:
        repo_root = Path(__file__).resolve().parent.parent
        output_dir = repo_root / "tests" / "input"

    output_dir.mkdir(parents=True, exist_ok=True)

    # Load reference sequence from the installed package resource
    ref_seq = (
        files("mitoedit.resources")
        .joinpath("mito.txt")
        .read_text()
        .replace("\n", "")
    )

    # Apply the deletion
    custom_seq = ref_seq[:DELETION_START] + ref_seq[DELETION_END:]

    # Plain-text .txt
    txt_path = output_dir / "custom_seq_200bp_deletion.txt"
    txt_path.write_text(custom_seq)

    # FASTA .fa  (60-character line width is conventional)
    fa_path = output_dir / "custom_seq_200bp_deletion.fa"
    lines = [FASTA_HEADER]
    for i in range(0, len(custom_seq), 60):
        lines.append(custom_seq[i : i + 60])
    fa_path.write_text("\n".join(lines) + "\n")

    print(f"Reference length : {len(ref_seq):,} bp")
    print(f"Custom length    : {len(custom_seq):,} bp  (deleted {DELETION_END - DELETION_START} bases)")
    print(f"Written .txt     : {txt_path}")
    print(f"Written .fa      : {fa_path}")

    return txt_path, fa_path


if __name__ == "__main__":
    generate()
