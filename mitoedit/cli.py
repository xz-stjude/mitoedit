import argparse
import os
from importlib.resources import files

import pandas as pd

from . import process_mitoedit
from . import MIN_SPACER, MAX_SPACER, ARR_MIN, ARR_MAX
from .io import read_sequence_file
import logging
import sys


def main():
    """CLI entry point for MitoEdit."""
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s",
        stream=sys.stdout,
    )
    logger = logging.getLogger("mitoedit")

    parser = argparse.ArgumentParser(
        description="Process DNA sequence for base editing."
    )
    # yapf: disable
    parser.add_argument('--mtdna_seq_path', '-i', type=str, default=None,       help='File containing the mtDNA sequence as plain text.')
    parser.add_argument('--bystander_file'      , type=str,                     help='Excel file containing bystander effect annotations (optional, for human mtDNA analysis)',
                            default=os.path.join(
                                os.path.dirname(os.path.abspath(__file__)),
                                'resources',
                                "annotated_human_mtDNA_10022024_for_bystanders_EDITED.xlsx"
                            ))
    parser.add_argument('--output_prefix', '-o' , type=str, default='output',   help='Prefix for output CSV files (default: output)')
    parser.add_argument('--min_spacer'          , type=int, default=MIN_SPACER, help=f'Minimum spacer length for TALE-NT (default: {MIN_SPACER})')
    parser.add_argument('--max_spacer'          , type=int, default=MAX_SPACER, help=f'Maximum spacer length for TALE-NT (default: {MAX_SPACER})')
    parser.add_argument('--array_min'           , type=int, default=ARR_MIN,    help=f'Minimum array length for TALE-NT (default: {ARR_MIN})')
    parser.add_argument('--array_max'           , type=int, default=ARR_MAX,    help=f'Maximum array length for TALE-NT (default: {ARR_MAX})')
    parser.add_argument('position'              , type=int,                     help='Position of the base to be changed')
    parser.add_argument('mutant_base'           , type=str,                     help='Mutant base to be changed into')
    parser.add_argument('--write_excel'           , action="store_true",          help='Save windows and bystander effects in spreadsheeets of an excel file (final_{position}_{mutant base})')
    # yapf: enable
    args = parser.parse_args()

    if args.mtdna_seq_path is None:
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
    else:
        logger.info(f"Reading mtDNA sequence from file: {args.mtdna_seq_path}")
        mtdna_seq = read_sequence_file(args.mtdna_seq_path)

    bystander_df = None
    logger.info(f"Reading Bystander information from file: {args.bystander_file}")
    if args.bystander_file:
        bystander_file = os.path.abspath(args.bystander_file)
        if os.path.isfile(bystander_file):
            logger.info(f"Loading bystander data from {bystander_file}")
            bystander_df = pd.read_excel(bystander_file)
        else:
            logger.warning(
                f"Bystander file {bystander_file} does not exist. Skipping bystander information."
            )
        if bystander_df is not None and len(bystander_df):
            logger.info(
                f"Successfully got Bystander information for {len(bystander_df)} mutations."
            )

    talen_params = {
        "min_spacer": args.min_spacer,
        "max_spacer": args.max_spacer,
        "array_min": args.array_min,
        "array_max": args.array_max,
    }

    results = process_mitoedit(
        mtdna_seq=mtdna_seq,
        position=args.position,
        mutant_base=args.mutant_base,
        reference_base=args.reference_base,
        bystander_df=bystander_df,
        talen_params=talen_params,
    )

    if results["windows_df"].empty:
        logger.warning("No results generated. Exiting.")
        return

    os.makedirs(args.output_prefix, exist_ok=True)
    logger.info(f"Output directory created/verified: {args.output_prefix}")

    logger.info("Writing adjacent bases to FASTA file.")
    fasta_file = f"{args.output_prefix}/adjacent_bases.fasta"
    with open(fasta_file, "w") as file:
        file.write(results["fasta_content"])
        logger.info(f"Finished writing FASTA file to {fasta_file}.")

    pipeline_windows_csv = f"{args.output_prefix}/pipeline_windows.csv"
    pipeline_bystanders_csv = f"{args.output_prefix}/pipeline_bystanders.csv"

    logger.info(f"Writing pipeline windows data to {pipeline_windows_csv}.")
    results["windows_df"].to_csv(pipeline_windows_csv, index=False)

    if not results["bystanders_df"].empty:
        logger.info(f"Writing pipeline bystanders data to {pipeline_bystanders_csv}.")
        results["bystanders_df"].to_csv(pipeline_bystanders_csv, index=False)
    else:
        logger.info("No bystanders information available to write.")

    combined_windows_csv = f"{args.output_prefix}/all_windows.csv"
    combined_bystanders_csv = f"{args.output_prefix}/all_bystanders.csv"
    if args.write_excel:

        combined_excel = (
            f"{args.output_prefix}/final_{args.position}_{args.mutant_base}.xlsx"
        )
        excel_writer = pd.ExcelWriter(combined_excel, engine="xlsxwriter")
    logger.info(f"Saving combined windows data to {combined_windows_csv}.")
    results["windows_df"].to_csv(combined_windows_csv, index=False)
    if args.write_excel:
        results["windows_df"].to_excel(
            excel_writer, sheet_name="all_windows", index=False
        )

    if not results["bystanders_df"].empty:
        logger.info(f"Saving combined bystanders data to {combined_bystanders_csv}.")
        results["bystanders_df"].to_csv(combined_bystanders_csv, index=False)
        if args.write_excel:
            results["bystanders_df"].to_excel(
                excel_writer, sheet_name="bystanders_effects", index=False
            )
    else:
        logger.info("No combined bystanders data to save.")

    if not results["talen_output_df"].empty:
        talen_output_path = f"{args.output_prefix}/talen_output.csv"
        logger.info(f"Saving talen output to {talen_output_path}.")
        results["talen_output_df"].to_csv(talen_output_path, index=False)
        if args.write_excel:
            results["talen_output_df"].to_excel(
                excel_writer, sheet_name="talen_output", index=False
            )
    if args.write_excel:
        excel_writer.close()
    logger.info("MitoEdit processing completed successfully.")
