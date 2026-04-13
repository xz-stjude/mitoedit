"""
Evaluate MitoEdit (current) bystander predictions against:
  - actual observed bystander counts (from curated dataset)
  - predictions from MitoEdit 1.0 (from curated dataset)

Matching is done on the exact window sequence from the dataset.

Output is saved to tests/output/bystander_evaluations.csv.

Usage:
    python scripts/evaluate_bystanders.py
"""

import logging
import re
import sys
from pathlib import Path

import pandas as pd

# Ensure repo root is on sys.path when run directly
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from mitoedit import process_mitoedit

logging.basicConfig(level=logging.WARNING, stream=sys.stderr)

DATASET_PATH = "/Users/ishaangupta/Lab/SJHack/data/mitoedit_curated_dataset/mitoedit_curated_dataset.csv"

# Map dataset pipeline names → new code pipeline names
PIPELINE_MAP = {
    "Mok2020_G1397": "Mok2020_Unified_G1397",
    "Mok2020_G1333": "Mok2020_Unified_G1333",
    "Mok2020_DddA11": "Mok2020_Unified_DddA11",
}

# In the dataset, alt_base is the REFERENCE base being edited.
REF_TO_MUTANT = {
    "C": "T",
    "G": "A",
    "T": "C",
    "A": "G",
}


def strip_markers(seq):
    """Remove [x], {x} bracket markers and lowercase for comparison."""
    return re.sub(r"[\[\]{}]", "", seq).lower()


TALEN_PARAMS = {
    "min_spacer": 13,
    "max_spacer": 18,
    "array_min": 14,
    "array_max": 18,
}


def run_mitoedit_cached(position, mutant_base, cache):
    """Run process_mitoedit and cache by (position, mutant_base)."""
    key = (position, mutant_base)
    if key not in cache:
        try:
            result = process_mitoedit(position=position, mutant_base=mutant_base, talen_params=TALEN_PARAMS)
            cache[key] = result["windows_df"]
        except Exception as e:
            print(f"  ERROR running mitoedit for pos={position} mut={mutant_base}: {e}", file=sys.stderr)
            cache[key] = None
    return cache[key]


def format_edit_list(total_count, positions):
    """Format as 'n - pos1, pos2, ...' matching the dataset convention."""
    if not positions:
        return str(int(total_count))
    pos_str = ", ".join(str(p) for p in sorted(positions))
    return f"{int(total_count)} - {pos_str}"


def lookup_by_window(windows_df, new_pipeline_name, win_bp, target_window_seq):
    """
    Find the row in windows_df matching pipeline, window size, and exact window sequence.
    Returns (bystander_count, raw_bystander_positions, marked_seq) or (None, None, None).
    """
    win_str = f"{win_bp}bp"
    target_clean = target_window_seq.strip().lower()

    candidates = windows_df[
        (windows_df["Pipeline"] == new_pipeline_name) &
        (windows_df["Window Size"] == win_str)
    ]

    for _, r in candidates.iterrows():
        if strip_markers(r["Window Sequence"]) == target_clean:
            count = int(r["Number of Bystanders"])
            positions = r["Position of Bystanders"]
            if not isinstance(positions, list):
                positions = []
            marked_seq = r["Window Sequence"]
            return count, positions, marked_seq

    return None, None, None


def extract_target_idx_from_marked(marked_window):
    """Return 0-based sequence index of the [X] target in a marked window."""
    seq_idx = 0
    i = 0
    while i < len(marked_window):
        c = marked_window[i]
        if c == '[':
            return seq_idx
        elif c in ('{', '}', ']'):
            i += 1
        else:
            seq_idx += 1
            i += 1
    return None


def parse_position_list(list_str):
    """Parse 'n - pos1 (x%), pos2, ...' into a set of int genomic positions, ignoring any parenthetical annotations."""
    s = str(list_str).strip()
    if s in ('', 'nan', 'None'):
        return set()
    if re.search(r' -+ ', s):
        _, pos_part = re.split(r' -+ ', s, maxsplit=1)
        positions = set()
        for token in pos_part.split(','):
            m = re.match(r'\s*(\d+)', token)
            if m:
                positions.add(int(m.group(1)))
        return positions
    return set()


def build_marked_window(raw_window, window_start, target_pos, all_edited_positions):
    """Build marked window: [X] for target, {X} for non-target edits."""
    bystander_set = all_edited_positions - {target_pos}
    output = []
    for i, base in enumerate(raw_window.upper()):
        gpos = window_start + i
        if gpos == target_pos:
            output.append(f"[{base}]")
        elif gpos in bystander_set:
            output.append(f"{{{base}}}")
        else:
            output.append(base)
    return "".join(output)


def main():
    df = pd.read_csv(DATASET_PATH)
    print(f"Loaded {len(df)} rows from dataset.\n")

    cache = {}
    results = []

    for _, row in df.iterrows():
        position = int(row["position"])
        alt_base = str(row["alt_base"]).strip()
        win_bp = int(row["win_bp"])
        pipeline = str(row["pipeline"]).strip()
        window_seq = str(row["window"]).strip()
        v1_pred = row["predicted_by_mitoEdit1.0"]
        v1_list = str(row.get("predicted_by_mitoEdit1.0_list", "")).strip()
        actual = row["actual"]
        actual_list = str(row.get("actual_list", "")).strip()

        new_pipeline = PIPELINE_MAP.get(pipeline)
        if new_pipeline is None:
            print(f"WARNING: Unknown pipeline '{pipeline}' for pos={position}, skipping.", file=sys.stderr)
            continue

        mutant_base = REF_TO_MUTANT.get(alt_base.upper())
        if mutant_base is None:
            print(f"WARNING: Cannot infer mutant for alt_base='{alt_base}', skipping.", file=sys.stderr)
            continue

        # actual includes the target position itself; subtract 1 for bystander-only count
        actual_bystanders = actual - 1
        v1_bystanders = v1_pred - 1 if v1_pred is not None else None

        print(f"Running pos={position} ref={alt_base} mut={mutant_base} ...", end=" ", flush=True)
        windows_df = run_mitoedit_cached(position, mutant_base, cache)

        if windows_df is None:
            print("FAILED")
            results.append({
                "position": position,
                "alt_base": alt_base,
                "win_bp": win_bp,
                "pipeline": pipeline,
                "window": window_seq,
                "actual_num_bystanders": actual_bystanders,
                "actual_list": actual_list,
                "v1.0_predicted_num_bystanders": v1_bystanders,
                "v1.0_predicted_list": v1_list,
                "v2_predicted_num_bystanders": None,
                "v2_predicted_list": None,
                "v2_marked_window": None,
                "v1.0_marked_window": None,
                "actual_marked_window": None,
                "window_found": False,
                "v2_matches_actual": None,
                "v1_matches_actual": v1_bystanders == actual_bystanders,
                "v1_matches_v2": None,
            })
            continue

        v2_count, v2_positions, v2_marked = lookup_by_window(windows_df, new_pipeline, win_bp, window_seq)

        # Build v2_predicted_list including the intended (target) edit
        if v2_count is not None:
            all_v2_edits = sorted([position] + list(v2_positions))
            v2_list = format_edit_list(v2_count + 1, all_v2_edits)
        else:
            v2_list = None

        # Build marked windows for v1.0 and actual using the target index from v2_marked_window
        v1_marked = None
        actual_marked = None
        if v2_marked is not None:
            tidx = extract_target_idx_from_marked(v2_marked)
            if tidx is not None:
                window_start = position - tidx
                raw = window_seq.upper()

                v1_positions = parse_position_list(v1_list)
                v1_marked = build_marked_window(raw, window_start, position, v1_positions)

                actual_positions = parse_position_list(actual_list)
                actual_marked = build_marked_window(raw, window_start, position, actual_positions)

        if v2_count is None:
            print(f"NO MATCH (window not found in v2 output)")
        else:
            print(f"done → v2={v2_count}, v1.0={v1_bystanders}, actual={actual_bystanders}")

        results.append({
            "position": position,
            "alt_base": alt_base,
            "win_bp": win_bp,
            "pipeline": pipeline,
            "window": window_seq,
            "actual_num_bystanders": actual_bystanders,
            "actual_list": actual_list,
            "v1.0_predicted_num_bystanders": v1_bystanders,
            "v1.0_predicted_list": v1_list,
            "v2_predicted_num_bystanders": v2_count,
            "v2_predicted_list": v2_list,
            "v2_marked_window": v2_marked,
            "v1.0_marked_window": v1_marked,
            "actual_marked_window": actual_marked,
            "window_found": v2_count is not None,
            "v2_matches_actual": (v2_count == actual_bystanders) if v2_count is not None else None,
            "v1_matches_actual": v1_bystanders == actual_bystanders,
            "v1_matches_v2": (v2_list == v1_list) if v2_list is not None else None,
        })

    results_df = pd.DataFrame(results)
    out_dir = _REPO_ROOT / "tests" / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "bystander_evaluations.csv"
    results_df.to_csv(out_path, index=False)

    print(f"\n{'='*70}")
    print("RESULTS SUMMARY")
    print(f"{'='*70}")
    print(results_df[[
        "position", "alt_base", "win_bp", "pipeline",
        "actual_num_bystanders", "actual_list",
        "v1.0_predicted_num_bystanders", "v1.0_predicted_list",
        "v2_predicted_num_bystanders", "v2_predicted_list",
        "window_found"
    ]].to_string(index=False))

    print(f"\n{'='*70}")
    print("ACCURACY SUMMARY")
    print(f"{'='*70}")

    found = results_df[results_df["window_found"]].copy()
    not_found = results_df[~results_df["window_found"]]
    print(f"Windows matched in v2 output : {len(found)}/{len(results_df)}")
    if len(not_found):
        print(f"Windows NOT found            : {len(not_found)} rows")
        for _, r in not_found.iterrows():
            print(f"  pos={r['position']} {r['pipeline']} {r['win_bp']}bp  seq={r['window']}")

    if len(found) > 0:
        v1_acc = found["v1_matches_actual"].mean() * 100
        v2_acc = found["v2_matches_actual"].mean() * 100
        v1_mae = (found["v1.0_predicted_num_bystanders"] - found["actual_num_bystanders"]).abs().mean()
        v2_mae = (found["v2_predicted_num_bystanders"] - found["actual_num_bystanders"]).abs().mean()

        print(f"\nOn matched rows ({len(found)}):")
        print(f"  Exact match rate  v1.0 : {v1_acc:.1f}%   v2 : {v2_acc:.1f}%")
        print(f"  MAE               v1.0 : {v1_mae:.3f}  v2 : {v2_mae:.3f}")

        v1_v2_match_rate = found["v1_matches_v2"].mean() * 100
        print(f"  v1 matches v2          : {v1_v2_match_rate:.1f}%")

        under = found[found["v2_predicted_num_bystanders"] < found["actual_num_bystanders"]]
        over  = found[found["v2_predicted_num_bystanders"] > found["actual_num_bystanders"]]

        if len(under):
            print(f"\nUNDER-PREDICTIONS ({len(under)} rows)  — v2 misses bystanders:")
            print(under[[
                "position", "pipeline", "win_bp",
                "actual_num_bystanders", "actual_list",
                "v1.0_predicted_num_bystanders", "v1.0_predicted_list",
                "v2_predicted_num_bystanders", "v2_predicted_list", "v2_marked_window"
            ]].to_string(index=False))

        if len(over):
            print(f"\nOVER-PREDICTIONS ({len(over)} rows)  — v2 finds extra bystanders:")
            print(over[[
                "position", "pipeline", "win_bp",
                "actual_num_bystanders", "actual_list",
                "v1.0_predicted_num_bystanders", "v1.0_predicted_list",
                "v2_predicted_num_bystanders", "v2_predicted_list", "v2_marked_window"
            ]].to_string(index=False))

    print(f"\nFull results saved to: {out_path.resolve()}")


if __name__ == "__main__":
    main()
