import random
import os
import pandas as pd
from abc import ABC, abstractmethod
import logging

logger = logging.getLogger(__name__)


class BasePipeline(ABC):
    """Abstract base class for all MitoEdit pipelines"""

    def __init__(
        self,
        min_window_size,
        max_window_size,
    ):
        self.pipeline_name = None  # To be set by subclasses
        self.min_window_size = min_window_size
        self.max_window_size = max_window_size

    @abstractmethod
    def process_mtDNA(self, mtDNA_seq, pos):
        """Main function which processes the DNA - must be implemented by subclasses

        Returns:
            tuple: (all_windows, adjacent_bases) where all_windows is a list of window data
                   and adjacent_bases is the sequence context around the target position
        """

    def process_bystander_data(self, all_windows_df, additional_df=None):
        """Process bystander information from additional DataFrame and return DataFrames."""
        logger.info("Processing bystander information.")

        # Process additional bystander data if provided
        if additional_df is not None and not additional_df.empty:
            # Extract bystander positions from the DataFrame column
            bystander_positions = set()
            for positions_list in all_windows_df["Position of Bystanders"]:
                if isinstance(positions_list, list):
                    bystander_positions.update(positions_list)

            filtered_df = additional_df[
                additional_df["mtDNA_pos"].isin(bystander_positions)
            ]
            new_data = filtered_df[
                [
                    "mtDNA_pos",
                    "Ref. Allele",
                    "Mutant Allele",
                    "Location",
                    "Predicted Impact",
                    "Syn vs NonSyn",
                    "AA Variant",
                    "Func. Impact",
                    "MutationAssessor Score",
                ]
            ]
            new_data.columns = [
                "Bystander Position",
                "Reference Base",
                "Mutant Base",
                "Location On Genome",
                "Predicted Mutation Impact",
                "SNV Type",
                "AA Variant",
                "Functional Impact",
                "MutationAssessor Score",
            ]
        else:
            logger.info("No additional bystander data provided.")
            new_data = pd.DataFrame()

        logger.info("Successfully processed bystander information.")

        return all_windows_df, new_data

    def _mark_bases(self, sequence, target_position, off_target_positions):
        """Mark the target and bystander bases in the window"""
        target_position -= 1
        off_target_positions = set(p - 1 for p in off_target_positions)
        marked_sequence = []
        for index, char in enumerate(sequence):
            if index == target_position:
                marked_sequence.append(
                    f"[{char}]"
                )  # Target base with square brackets []
            elif index in off_target_positions:
                marked_sequence.append(
                    f"{{{char}}}"
                )  # Off-target base with curly braces {}
            else:
                marked_sequence.append(char)  # No special marking
        return "".join(marked_sequence)

    def _mark_base_at_position(self, sequence, target_position):
        """Mark the base at the target position --> to mark the target base in the window"""
        logger.debug("Marking base at the specific position")
        marked_sequence = []
        for index, char in enumerate(sequence):
            if index == target_position:
                marked_sequence.append(
                    f"{{{char}}}"
                )  # Target base with curly brackets {}
            else:
                marked_sequence.append(char)  # No special marking
        return "".join(marked_sequence)

    def _reverse_complement(self, sequence):
        """Get the reverse complement of a DNA sequence"""
        logger.debug("Generating reverse complement.")
        complement_dict = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "}": "{",
            "{": "}",
            "[": "]",
            "]": "[",
        }
        reverse_complement_sequence = "".join(
            [
                complement_dict.get(
                    base, random.choice(["A", "T", "C", "G"]) if base == "N" else base
                )
                for base in sequence
            ]
        )
        return reverse_complement_sequence[::-1]

    def _complementing(self, sequence):
        """Get the complement of a DNA sequence."""
        logger.debug("Generating complement.")
        complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
        complement_sequence = "".join(
            [
                complement_dict.get(
                    base, random.choice(["A", "T", "C", "G"]) if base == "N" else base
                )
                for base in sequence
            ]
        )
        return complement_sequence

    def _create_window(self, mtDNA_seq, pos, start_index, end_index):
        """To create the window"""
        logger.debug(f"Creating window from position {pos}.")
        circular_seq = mtDNA_seq + mtDNA_seq
        window = circular_seq[start_index:end_index]
        return window

    def _remove_whitespace(self, sequence):
        """Remove whitespace from the sequence"""
        logger.debug("Removing whitespace from the sequence.")
        return sequence.replace(" ", "").replace("\t", "").replace("\n", "")

    def _capitalize(self, sequence):
        """Capitalize the sequence"""
        logger.debug("Capitalizing the sequence.")
        return sequence.upper()

    def _find_dinucs(
        self, sequence, dinucs, target_base, double_match_index=1, unique=True
    ):
        """Find all positions where the specified dinucleotides occur"""
        logger.debug(f"Finding dinucleotides {dinucs} in sequence.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i : i + 2] in dinucs:
                # add 1-based indexing for the target base
                if sequence[i] == target_base:
                    if sequence[i + 1] == target_base:  # double match
                        assert double_match_index in [
                            1,
                            2,
                        ], "double_match_index must be 1 or 2"
                        positions.append(
                            i + double_match_index
                        )  # if double match, return double_match_index (1-based)
                    else:
                        positions.append(i + 1)
                elif sequence[i + 1] == target_base:
                    positions.append(i + 2)
                else:
                    raise ValueError(
                        f"given dinucleotide {sequence[i:i+2]} does not contain the target base is {target_base}"
                    )
        return list(set(positions)) if unique else positions

    def _find_N_positions(self, window, start_position, pattern):
        """Find positions of a specific pattern in the window"""
        logger.debug(f"Finding positions of pattern {pattern} in sequence.")
        positions = []
        for i in range(len(window) - len(pattern) + 1):
            if window[i : i + len(pattern)] == pattern:
                positions.append(start_position + i)
        return positions
