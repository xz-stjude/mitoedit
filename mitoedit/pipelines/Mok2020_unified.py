import logging

logger = logging.getLogger(__name__)
from .base_pipeline import BasePipeline

WINDOW_SIZE_MIN = 14
WINDOW_SIZE_MAX = 19


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
        self.pipeline_name = "Mok2020_Unified"

    def _get_g1397_position_range(self, window_size):
        """Get the position range for G1397 window generation."""
        return range(4, window_size - 3)

    def _get_g1333_position_range(self, window_size):
        """Get the position range for G1333 window generation."""
        return range(3, window_size - 2)

    def _get_ddda11_position_range(self, window_size):
        """Get the position range for DddA11 window generation."""
        return range(4, window_size - 3)

    def _process_context_all_variants(
        self, nospace_mtDNA, pos, context_positions, ref_base, mut_base, edit_type
    ):
        """Process a context using all three positioning strategies."""
        all_windows = []

        # Generate circular sequence for window extraction
        circular_seq = nospace_mtDNA + nospace_mtDNA

        # Calculate adjacent bases (same for all variants)
        start_index = pos - 31
        end_index = pos + 30
        adjacent_bases = circular_seq[start_index:end_index]

        # Process with all three positioning strategies
        strategies = [
            ("G1397", self._get_g1397_position_range),
            ("G1333", self._get_g1333_position_range),
            ("DddA11", self._get_ddda11_position_range),
        ]

        for strategy_name, position_range_func in strategies:
            for window_size in range(self.min_window_size, self.max_window_size + 1):
                position_range = position_range_func(window_size)

                for target_pos in position_range:
                    start_pos = pos - target_pos
                    end_pos = start_pos + window_size

                    if start_pos < 1 or end_pos > len(nospace_mtDNA):
                        continue

                    window = circular_seq[start_pos - 1 : end_pos - 1]

                    # Find bystander positions in this window
                    bystander_positions = []
                    for ctx_pos in context_positions:
                        if start_pos <= ctx_pos <= end_pos and ctx_pos != pos:
                            bystander_positions.append(ctx_pos)

                    # Create window data tuple
                    window_data = (
                        f"{self.pipeline_name}_{strategy_name}",  # Pipeline variant name
                        pos,  # Target position
                        ref_base,  # Reference base
                        mut_base,  # Mutant base
                        f"{window_size}bp",  # Window size
                        self._mark_bases(
                            window[1:],
                            pos - start_pos,
                            [(pos - start_pos) for pos in bystander_positions],
                        ),  # Window sequence
                        f"Position {target_pos}",  # Target position in window
                        len(bystander_positions),  # Bystander count
                        bystander_positions,  # Bystander positions
                        edit_type,  # Edit type description
                        strategy_name,  # Strategy identifier
                    )

                    all_windows.append(window_data)

        return all_windows, adjacent_bases

    def process_mtDNA(self, mtDNA_seq, pos):
        """Main function which processes the DNA using all Mok2020 variants."""
        logger.info(
            f"Processing mtDNA sequence for position {pos} using unified Mok2020 pipeline."
        )

        nospace_mtDNA = self._capitalize(self._remove_whitespace(mtDNA_seq))

        # Find all context positions (from all variants)
        C_CONTEXT = ["TC", "AC", "CC"]
        G_CONTEXT = ["GA", "GT", "GG"]

        C_context_positions = self._find_dinucs(nospace_mtDNA, C_CONTEXT, "C", 2)
        G_context_positions = self._find_dinucs(nospace_mtDNA, G_CONTEXT, "G", 1)
        logger.info(f"C_context_positions:{C_context_positions}")
        logger.info(f"G_context_positions:{G_context_positions}")
        if pos in C_context_positions:
            logger.info(
                f"Base at position {pos} is in a 5'-{nospace_mtDNA[pos-1:pos+1]} context."
            )
            return self._process_context_all_variants(
                nospace_mtDNA,
                pos,
                C_context_positions + G_context_positions,
                "C",
                "T",
                f"C→T (5'-{nospace_mtDNA[pos-1:pos+1]} context)",
            )

        if pos in G_context_positions:
            logger.info(
                f"Base at position {pos} is in a 5'-{nospace_mtDNA[pos-1:pos+1]} context."
            )
            return self._process_context_all_variants(
                nospace_mtDNA,
                pos,
                C_context_positions + G_context_positions,
                "G",
                "A",
                f"G→A (5'-{nospace_mtDNA[pos-1:pos+1]} context)",
            )

        consecutive_TC_positions = self._find_consecutive_TC_sequences(nospace_mtDNA)
        consecutive_AC_positions = self._find_consecutive_AC_sequences(nospace_mtDNA)
        consecutive_CC_positions = self._find_consecutive_CC_sequences(nospace_mtDNA)
        consecutive_GA_positions = self._find_consecutive_GA_sequences(nospace_mtDNA)
        consecutive_GT_positions = self._find_consecutive_GT_sequences(nospace_mtDNA)
        consecutive_GG_positions = self._find_consecutive_GG_sequences(nospace_mtDNA)

        # Check which context the position belongs to and process with all variants
        if pos in consecutive_TC_positions:
            logger.info(f"Base at position {pos} is in a 5'-TC context.")
            return self._process_context_all_variants(
                nospace_mtDNA,
                pos,
                consecutive_TC_positions,
                "C",
                "T",
                "C→T (TC context)",
            )
        elif pos in consecutive_AC_positions:
            logger.info(f"Base at position {pos} is in a 5'-AC context.")
            return self._process_context_all_variants(
                nospace_mtDNA,
                pos,
                consecutive_AC_positions,
                "C",
                "T",
                "C→T (AC context)",
            )
        elif pos in consecutive_CC_positions:
            logger.info(f"Base at position {pos} is in a 5'-CC context.")
            return self._process_context_all_variants(
                nospace_mtDNA,
                pos,
                consecutive_CC_positions,
                "C",
                "T",
                "C→T (CC context)",
            )
        elif pos in consecutive_GA_positions:
            logger.info(f"Base at position {pos} is in a 5'-GA context.")
            return self._process_context_all_variants(
                nospace_mtDNA,
                pos,
                consecutive_GA_positions,
                "G",
                "A",
                "G→A (GA context)",
            )
        elif pos in consecutive_GT_positions:
            logger.info(f"Base at position {pos} is in a 5'-GT context.")
            return self._process_context_all_variants(
                nospace_mtDNA,
                pos,
                consecutive_GT_positions,
                "G",
                "A",
                "G→A (GT context)",
            )
        elif pos in consecutive_GG_positions:
            logger.info(f"Base at position {pos} is in a 5'-GG context.")
            return self._process_context_all_variants(
                nospace_mtDNA,
                pos,
                consecutive_GG_positions,
                "G",
                "A",
                "G→A (GG context)",
            )
        else:
            logger.warning(
                f"Base at position {pos} is not in any editable context for Mok2020 pipelines."
            )
            return [], []
