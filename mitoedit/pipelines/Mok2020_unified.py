import logging

logger = logging.getLogger(__name__)
from .base_pipeline import BasePipeline

WINDOW_SIZE_MIN = 14
WINDOW_SIZE_MAX = 19


DdCBE_G1397 = "DdCBE_G1397"
DdCBE_G1333 = "DdCBE_G1333"
DddA11_G1397 = "DddA11_G1397"
STRATEGIES = [DdCBE_G1397, DdCBE_G1333, DddA11_G1397]


def revcomp(seq):
    return "".join(
        [{"A": "T", "T": "A", "C": "G", "G": "C"}[c] for c in seq.upper()][::-1]
    )


class Mok2020UnifiedPipeline(BasePipeline):
    """
    Unified Mok2020 Base Editing Pipeline (All Variants Combined)

    Combines all Mok2020 base editing approaches (G1397, G1333, and DddA11) into a single
    comprehensive pipeline. This unified approach generates windows using all positioning
    strategies and supports all sequence contexts from the original three pipelines.

    Supported Editing Contexts:
    - C→T edits in: TC, AC, CC contexts (5'-XC-3' where X = T, A, or C)
    - G→A edits in: GA, GT, GG contexts (5'-GX-3' where X = A, T, or G)

    Positioning Strategies:
    - G1397 strategy: positions 4 to window_size-4 (conservative)
    - G1333 strategy: positions 3 to window_size-3 (aggressive)
    - DddA11 strategy: positions 4 to window_size-4 with extended contexts

    Key Features:
    - Combines all three Mok2020 variants in a single pipeline
    - Generates windows using all positioning strategies
    - Supports all sequence contexts from original pipelines
    - Returns combined results from all sub-approaches
    - Single adjacent_bases output (same for all variants)
    """

    def __init__(
        self,
        min_window_size: int = WINDOW_SIZE_MIN,
        max_window_size: int = WINDOW_SIZE_MAX,
    ):
        super().__init__(min_window_size, max_window_size)
        self.pipeline_name = "Mok2020"

    def _get_position_range(self, strategy_name, window_size):
        """Get the 1-based position range (inclusive) for target site within the spacer (from 5' end)."""
        position_ranges = {
            DdCBE_G1397: (window_size - 7 + 1, window_size - 4 + 1),
            DdCBE_G1333: (4, 10),
            DddA11_G1397: (window_size - 7 + 1, window_size - 4 + 1),
        }
        return position_ranges[strategy_name]

    def _get_context(self, strategy_name):
        """Get the dinucleotide context for target and bystanders."""
        CONTEXTS = {
            DdCBE_G1397: (["TC"], ["GA"]),
            DdCBE_G1333: (["TC"], ["GA"]),
            DddA11_G1397: (["TC", "AC", "CC"], ["GA", "GT", "GG"]),
        }
        return CONTEXTS[strategy_name]

    def _get_C_context(self, strategy_name):
        """Get the 5' dinucleotide context for target and bystanders."""
        return self._get_context(strategy_name)[0]

    def _get_G_context(self, strategy_name):
        """Get the 3' dinucleotide context for target and bystanders."""
        return self._get_context(strategy_name)[1]

    def _process_context_strategy(
        self,
        strategy_name: str,
        nospace_mtDNA,
        pos,
        ref_base,
        mut_base,
        edit_forward: bool,
    ):
        """Process a context using a particular positioning strategy."""
        all_windows = []
        circular_seq = nospace_mtDNA + nospace_mtDNA

        for window_size in range(self.min_window_size, self.max_window_size + 1):
            editable_start, editable_end = self._get_position_range(
                strategy_name, window_size
            )  # 1-based
            forward_pos_range = (editable_start - 1, editable_end)  # 0-based, exclusive
            reverse_pos_range = (
                window_size - editable_end,
                window_size - editable_start + 1,
            )  # 0-based, exclusive
            for target_pos in (
                range(*forward_pos_range) if edit_forward else range(*reverse_pos_range)
            ):
                start_pos = pos - target_pos - 1
                end_pos = start_pos + window_size

                if start_pos < 1 or end_pos > len(nospace_mtDNA):
                    continue

                window = circular_seq[start_pos:end_pos]
                # Find bystander positions in this window
                bystander_positions = []
                for ctx_pos in self._find_dinucs(
                    window[forward_pos_range[0] - 1 : forward_pos_range[1]],
                    self._get_C_context(strategy_name),
                    "C",
                    2,
                    True,
                    start_pos + forward_pos_range[0] - 1,
                ) + self._find_dinucs(
                    window[reverse_pos_range[0] : reverse_pos_range[1] + 1],
                    self._get_G_context(strategy_name),
                    "G",
                    1,
                    True,
                    start_pos + reverse_pos_range[0],
                ):
                    if ctx_pos != pos:
                        bystander_positions.append(ctx_pos)
                # Create window data tuple
                window_data = (
                    f"{self.pipeline_name}_{strategy_name}",  # Pipeline variant name
                    pos,  # Target position
                    ref_base,  # Reference base
                    mut_base,  # Mutant base
                    f"{window_size}bp",  # Window size
                    self._mark_bases(
                        window,
                        pos - start_pos,
                        [(pos - start_pos) for pos in bystander_positions],
                    ),  # Window sequence
                    f"Position {target_pos+1}",  # Target position in window
                    len(bystander_positions),  # Bystander count
                    bystander_positions,  # Bystander positions
                    f"{ref_base}→{mut_base} (5'-{nospace_mtDNA[pos-2:pos] if edit_forward else revcomp(nospace_mtDNA[pos-1:pos+1])} context)",  # Edit type description
                    strategy_name,  # Strategy identifier
                    None,
                    None,
                )
                all_windows.append(window_data)

        return all_windows

    def process_mtDNA(self, mtDNA_seq, pos):
        """Main function which processes the DNA using all Mok2020 variants."""
        logger.info(
            f"Processing mtDNA sequence for position {pos} using unified Mok2020 pipeline."
        )

        nospace_mtDNA = self._capitalize(self._remove_whitespace(mtDNA_seq))

        # Calculate adjacent bases (same for all variants)
        start_index = pos - 31
        end_index = pos + 30
        circular_seq = nospace_mtDNA + nospace_mtDNA
        adjacent_bases = circular_seq[start_index:end_index]

        all_windows = []
        for strategy_name in STRATEGIES:
            if nospace_mtDNA[pos - 2 : pos] in self._get_C_context(strategy_name):
                edit_forward = True
                logger.info(
                    f"Base at position {pos} is in a 5'-{nospace_mtDNA[pos-2:pos]} context."
                )
            elif nospace_mtDNA[pos - 1 : pos + 1] in self._get_G_context(strategy_name):
                edit_forward = False
                logger.info(
                    f"Base at position {pos} is in a 5'-{revcomp(nospace_mtDNA[pos-1:pos+1])} context."
                )
            else:
                continue
            all_windows.extend(
                self._process_context_strategy(
                    strategy_name,
                    nospace_mtDNA,
                    pos,
                    "C" if edit_forward else "G",
                    "T" if edit_forward else "A",
                    edit_forward,
                )
            )
        if len(all_windows) == 0:
            logger.warning(
                f"Base at position {pos} is not in any editable context for Mok2020 pipelines."
            )
        return all_windows, adjacent_bases
