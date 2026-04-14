import io
import os
import re
import tempfile
from importlib.resources import files
from pathlib import Path
from types import NoneType

import pandas as pd

from .pipelines.Cho_sTALEDs import ChosTALEDsPipeline
from .pipelines.Mok2020_unified import Mok2020UnifiedPipeline
from .alignment import build_position_maps, window_has_indels
from talenWF import FindTALTask

import logging

logger = logging.getLogger(__name__)

MIN_SPACER = 14
MAX_SPACER = 18
ARR_MIN = 14
ARR_MAX = 18


class ReferenceBaseError(Exception):
    """Raised when the provided reference base does not match the base at the given position in the sequence."""

    pass


class PipelineError(Exception):
    """Raised when"""

    pass


class CommandError(Exception):
    """Raised when"""

    pass


class MitoEditError(Exception):
    """Raised when"""

    pass


def process_mitoedit(
    position,
    mutant_base,
    mtdna_seq=None,
    reference_base=None,
    ref_annot_df=None,
    talen_params={
        "min_spacer": MIN_SPACER,
        "max_spacer": MAX_SPACER,
        "array_min": ARR_MIN,
        "array_max": ARR_MAX,
    },
):
    """
    Core MitoEdit processing function for programmatic use.

    Args:
        mtdna_seq (str): mtDNA sequence string
        position (int): Position of the base to be changed (1-based)
        mutant_base (str): Mutant base to be changed into
        ref_annot_df (pd.DataFrame, optional): Reference annotation DataFrame for bystander effect lookup (human mtDNA annotation). Auto-loaded from bundled resources when a custom sequence is provided.
        talen_params (dict, optional): TALE-NT parameters for findTAL analysis

    Returns:
        dict: Results containing windows_df, bystanders_df, adjacent_bases, fasta_content, and talen_output_df
    """
    is_custom_sequence = mtdna_seq is not None

    if not is_custom_sequence:
        logger.info("Using default mtDNA sequence from resources/mito.txt")
        try:
            mtdna_seq = (
                files("mitoedit.resources")
                .joinpath("mito.txt")
                .read_text()
                .replace("\n", "")
            )
        except FileNotFoundError:
            logger.error("Default mtDNA sequence file not found in resources/mito.txt")
            raise

    if is_custom_sequence and ref_annot_df is None:
        try:
            with files("mitoedit.resources").joinpath(
                "annotated_human_mtDNA_10022024_for_bystanders_EDITED.xlsx"
            ).open("rb") as fh:
                ref_annot_df = pd.read_excel(fh)
            logger.info("Auto-loaded bundled bystander annotation: %d rows.", len(ref_annot_df))
        except Exception as exc:
            logger.warning("Could not auto-load bundled bystander annotation: %s", exc)
    mutant_base = mutant_base.upper()

    reference_base_in_seq = mtdna_seq[position - 1].upper()
    logger.info(f"Reference base at position {position} is {reference_base_in_seq}")

    if reference_base is not None:
        if reference_base.upper() != reference_base_in_seq:
            raise ReferenceBaseError(
                f"Reference base mismatch at position {position}: "
                f"expected '{reference_base.upper()}' but found '{reference_base_in_seq}' in the sequence."
            )

    reference_base = reference_base_in_seq

    custom_to_ref = None
    ref_to_custom = None
    if is_custom_sequence:
        _ref_seq = (
            files("mitoedit.resources")
            .joinpath("mito.txt")
            .read_text()
            .replace("\n", "")
        )
        logger.info("Running global alignment: custom sequence → reference mtDNA.")
        custom_to_ref, ref_to_custom = build_position_maps(mtdna_seq, _ref_seq)

    if (reference_base, mutant_base) in [("C", "T"), ("G", "A")]:
        pipeline_class = Mok2020UnifiedPipeline
        pipeline_name = "Mok2020_Unified"
    elif (reference_base, mutant_base) in [("A", "G"), ("T", "C")]:
        pipeline_class = ChosTALEDsPipeline
        pipeline_name = "Cho_sTALEDs"
    else:
        raise ValueError(
            f"No pipeline found for reference base {reference_base} and the mutant base {mutant_base}"
        )

    logger.info(f"Selected pipeline: {pipeline_name}")

    pipeline_instance = pipeline_class(
        talen_params["min_spacer"], talen_params["max_spacer"]
    )

    logger.info(f"Processing mtDNA sequence for position {position}.")
    all_windows, adjacent_bases = pipeline_instance.process_mtDNA(mtdna_seq, position)

    if not adjacent_bases:
        raise ValueError(f"The base found at position {position} cannot be edited.")

    windows_df = pd.DataFrame(
        all_windows,
        columns=[
            "Pipeline",
            "Position",
            "Reference Base",
            "Mutant Base",
            "Window Size",
            "Window Sequence",
            "Target Location",
            "Number of Bystanders",
            "Position of Bystanders",
            "Edit type",
            "Strategy Name",
            "Optimal Flanking TALEs",
            "Flag (CheckBystanderEffect)",
        ],
    )

    fasta_content = f">Adjacent_bases_position_{position}\n{adjacent_bases}\n"

    logger.info("Processing pipeline data.")
    if is_custom_sequence and custom_to_ref is not None:
        # Lift bystander positions to reference coordinates
        windows_df["Reference Position of Bystanders"] = windows_df[
            "Position of Bystanders"
        ].apply(
            lambda lst: [custom_to_ref.get(p) for p in lst] if isinstance(lst, list) else []
        )
        # Flag windows that span an indel
        windows_df["Window Has Indels"] = windows_df.apply(
            lambda row: window_has_indels(row, custom_to_ref), axis=1
        )
        # Look up bystander effects using reference positions
        if ref_annot_df is not None and not ref_annot_df.empty:
            ref_positions = {
                rp
                for lst in windows_df["Reference Position of Bystanders"]
                for rp in lst
                if rp is not None
            }
            filtered = ref_annot_df[ref_annot_df["mtDNA_pos"].isin(ref_positions)].copy()
            bystanders_df = filtered[
                [
                    "mtDNA_pos", "Ref. Allele", "Mutant Allele", "Location",
                    "Predicted Impact", "Syn vs NonSyn", "AA Variant",
                    "Func. Impact", "MutationAssessor Score",
                ]
            ].copy()
            bystanders_df.columns = [
                "Reference Bystander Position", "Reference Base", "Mutant Base",
                "Location On Genome", "Predicted Mutation Impact", "SNV Type",
                "AA Variant", "Functional Impact", "MutationAssessor Score",
            ]
            bystanders_df.insert(
                0,
                "Custom Bystander Position",
                bystanders_df["Reference Bystander Position"].map(
                    lambda rp: ref_to_custom.get(rp)
                ),
            )
            bystanders_df = bystanders_df.reset_index(drop=True)
        else:
            bystanders_df = pd.DataFrame()
    else:
        # Original reference-sequence path — untouched
        windows_df, bystanders_df = pipeline_instance.process_bystander_data(
            windows_df, ref_annot_df
        )

    logger.info("Pipeline processing completed successfully.")
    logger.info("Pipeline processing completed successfully.")

    if talen_params is None:
        talen_params = {
            "min_spacer": MIN_SPACER,
            "max_spacer": MAX_SPACER,
            "array_min": ARR_MIN,
            "array_max": ARR_MAX,
        }

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as temp_fasta:
        temp_fasta.write(fasta_content)
        temp_fasta_path = temp_fasta.name

    try:
        logger.info("Running talenWF findTAL analysis")
        talen_output_df = FindTALTask(
            fasta=temp_fasta_path,
            min_spacer=talen_params["min_spacer"],
            max_spacer=talen_params["max_spacer"],
            array_min=talen_params["array_min"],
            array_max=talen_params["array_max"],
            outpath=None,
        ).run()
        logger.info("talenWF analysis completed successfully")

        if "Plus strand sequence" not in talen_output_df.columns:
            logger.warning("Column 'Plus strand sequence' not found in the TALEN file")
            talen_output_df = pd.DataFrame()
        else:
            plus_strand_sequences = talen_output_df["Plus strand sequence"].tolist()
            txt_bases = set(re.findall(r"[a-z]+", " ".join(plus_strand_sequences)))

            for index, row in windows_df.iterrows():
                cleaned_sequence = re.sub(
                    r"[{}\[\]]", "", row["Window Sequence"]
                ).lower()
                windows_df.at[index, "Optimal Flanking TALEs"] = (
                    cleaned_sequence in txt_bases
                )
                for txt_base in txt_bases:
                    if txt_base == cleaned_sequence:
                        left_index = 1
                        right_index = 1
                        matching_sequences = [
                            sequence
                            for sequence in plus_strand_sequences
                            if re.sub(r"[^a-z]", "", sequence) == txt_base
                        ]
                        for sequence in matching_sequences:
                            lower_indices = [
                                i for i, char in enumerate(sequence) if char.islower()
                            ]
                            if len(lower_indices) >= 2:
                                spacer_start = lower_indices[0]
                                spacer_end = lower_indices[-1]
                                left_tale = sequence[:spacer_start].upper()
                                right_tale = sequence[spacer_end + 1 :].upper()
                                windows_df.at[index, f"Left TALE {left_index}"] = (
                                    left_tale
                                )
                                windows_df.at[index, f"Right TALE {right_index}"] = (
                                    right_tale
                                )
                                left_index += 1
                                right_index += 1

    finally:
        os.unlink(temp_fasta_path)

    logger.info("All processing completed successfully.")

    return {
        "windows_df": windows_df,
        "bystanders_df": bystanders_df,
        "adjacent_bases": adjacent_bases,
        "fasta_content": fasta_content,
        "talen_output_df": talen_output_df,
    }
