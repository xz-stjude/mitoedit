import os
import re
import tempfile
from importlib.resources import files

import pandas as pd

from .pipelines.Cho_sTALEDs import ChosTALEDsPipeline
from .pipelines.Mok2020_unified import Mok2020UnifiedPipeline
from talenWF import FindTALTask

import logging

logger = logging.getLogger(__name__)

MIN_SPACER = 14
MAX_SPACER = 18
ARR_MIN = 14
ARR_MAX = 18


class ReferenceBaseError(Exception):
    """Raised when"""

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
    bystander_df=None,
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
        bystander_df (pd.DataFrame, optional): DataFrame containing bystander effect annotations
        talen_params (dict, optional): TALE-NT parameters for findTAL analysis

    Returns:
        dict: Results containing windows_df, bystanders_df, adjacent_bases, fasta_content, and talen_output_df
    """
    if mtdna_seq is None:
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
    mutant_base = mutant_base.upper()

    reference_base = mtdna_seq[position - 1].upper()
    logger.info(f"Reference base at position {position} is {reference_base}")

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
            "Optimal Flanking TALEs",
            "Flag (CheckBystanderEffect)",
        ],
    )

    fasta_content = f">Adjacent_bases_position_{position}\n{adjacent_bases}\n"

    logger.info("Processing pipeline data.")
    windows_df, bystanders_df = pipeline_instance.process_bystander_data(
        windows_df, bystander_df
    )

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
                windows_df.at[index, "Matching TALEs"] = cleaned_sequence in txt_bases
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
