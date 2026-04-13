"""
alignment.py — Global pairwise alignment utilities for MitoEdit.

Aligns a custom mtDNA sequence to the reference and builds bidirectional
1-indexed position maps between the two sequences.
"""

import logging
import re

from Bio.Align import PairwiseAligner

logger = logging.getLogger(__name__)

_MATCH_SCORE = 2
_MISMATCH_SCORE = -1
_OPEN_GAP_SCORE = -5
_EXTEND_GAP_SCORE = -0.5


def build_position_maps(custom_seq: str, ref_seq: str) -> tuple:
    """Align custom_seq to ref_seq and return bidirectional 1-indexed position maps.

    Parameters
    ----------
    custom_seq : str
        User-supplied mtDNA sequence (whitespace already stripped).
    ref_seq : str
        Reference mtDNA sequence from resources/mito.txt.

    Returns
    -------
    custom_to_ref : dict
        1-indexed custom position → 1-indexed ref position, or None (insertion).
    ref_to_custom : dict
        1-indexed ref position → 1-indexed custom position, or None (deletion).
    """
    custom_seq = custom_seq.upper()
    ref_seq = ref_seq.upper()

    logger.info(
        "Building alignment position maps: custom_len=%d, ref_len=%d",
        len(custom_seq), len(ref_seq),
    )

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = _MATCH_SCORE
    aligner.mismatch_score = _MISMATCH_SCORE
    aligner.open_gap_score = _OPEN_GAP_SCORE
    aligner.extend_gap_score = _EXTEND_GAP_SCORE

    # Lazy iterator — only materialise the single best alignment
    alignment = next(iter(aligner.align(custom_seq, ref_seq)))
    logger.info("Alignment score: %.1f", alignment.score)

    # alignment.aligned[0] — ungapped blocks on the query  (custom_seq, 0-based half-open)
    # alignment.aligned[1] — ungapped blocks on the target (ref_seq,    0-based half-open)
    custom_to_ref: dict = {}
    ref_to_custom: dict = {}

    for (q_start, q_end), (t_start, t_end) in zip(
        alignment.aligned[0], alignment.aligned[1]
    ):
        for offset in range(q_end - q_start):
            c_pos = q_start + offset + 1   # 1-indexed
            r_pos = t_start + offset + 1   # 1-indexed
            custom_to_ref[c_pos] = r_pos
            ref_to_custom[r_pos] = c_pos

    # Positions absent from any block are insertions / deletions
    for c_pos in range(1, len(custom_seq) + 1):
        custom_to_ref.setdefault(c_pos, None)

    for r_pos in range(1, len(ref_seq) + 1):
        ref_to_custom.setdefault(r_pos, None)

    n_ins = sum(1 for v in custom_to_ref.values() if v is None)
    n_del = sum(1 for v in ref_to_custom.values() if v is None)
    logger.info("Position map built: %d insertions, %d deletions in custom seq.", n_ins, n_del)

    return custom_to_ref, ref_to_custom


def window_has_indels(window_row, custom_to_ref: dict) -> bool:
    """Return True if the window's aligned span contains any indel.

    Derives the window start from the marked Window Sequence string so the
    logic is pipeline-agnostic (works for Mok2020 and Cho formats alike).

    Parameters
    ----------
    window_row : pandas Series
        Row from windows_df with columns 'Position', 'Window Size',
        and 'Window Sequence'.
    custom_to_ref : dict
        custom → ref position map from build_position_maps().
    """
    window_seq = str(window_row["Window Sequence"])
    pos = int(window_row["Position"])
    window_size = int(str(window_row["Window Size"]).replace("bp", ""))

    # Count alpha characters before the '[' target marker to get the
    # 0-indexed offset of the target within the window.
    target_marker = window_seq.index("[")
    target_0idx = len(re.sub(r"[^a-zA-Z]", "", window_seq[:target_marker]))

    # 1-indexed first custom position covered by this window
    window_start = pos - target_0idx

    ref_positions = [custom_to_ref.get(p) for p in range(window_start, window_start + window_size)]

    # Insertion: any custom position has no reference equivalent
    if any(r is None for r in ref_positions):
        return True

    # Deletion: reference positions are not all consecutive
    non_none = sorted(r for r in ref_positions if r is not None)
    for i in range(1, len(non_none)):
        if non_none[i] != non_none[i - 1] + 1:
            return True

    return False
