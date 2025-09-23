import argparse
import os
from importlib.resources import files

import pandas as pd

from . import process_mitoedit

import logging
import sys

MIN_SPACER = 14
MAX_SPACER = 18
ARR_MIN = 14
ARR_MAX = 18
FILTER = 1
CUT_POS = 31


def main():
    """CLI entry point for MitoEdit."""
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s', stream=sys.stdout)
    logger = logging.getLogger('mitoedit')
    
    parser = argparse.ArgumentParser(description='Process DNA sequence for base editing.')
    # yapf: disable
    parser.add_argument('--mtdna_seq_path', '-i', type=str, default=None,       help='File containing the mtDNA sequence as plain text.')
    parser.add_argument('--bystander_file'      , type=str,                     help='Excel file containing bystander effect annotations (optional, for human mtDNA analysis)')
    parser.add_argument('--output_prefix', '-o' , type=str, default='output',   help='Prefix for output CSV files (default: output)')
    parser.add_argument('--min_spacer'          , type=int, default=MIN_SPACER, help=f'Minimum spacer length for TALE-NT (default: {MIN_SPACER})')
    parser.add_argument('--max_spacer'          , type=int, default=MAX_SPACER, help=f'Maximum spacer length for TALE-NT (default: {MAX_SPACER})')
    parser.add_argument('--array_min'           , type=int, default=ARR_MIN,    help=f'Minimum array length for TALE-NT (default: {ARR_MIN})')
    parser.add_argument('--array_max'           , type=int, default=ARR_MAX,    help=f'Maximum array length for TALE-NT (default: {ARR_MAX})')
    parser.add_argument('--filter'              , type=int, default=FILTER,     help=f'TALE-NT filter setting (default: {FILTER})')
    parser.add_argument('--cut_pos'             , type=int, default=CUT_POS,    help=f'TALE-NT cut position (default: {CUT_POS})')
    parser.add_argument('position'              , type=int,                     help='Position of the base to be changed')
    parser.add_argument('mutant_base'           , type=str,                     help='Mutant base to be changed into')
    # yapf: enable
    args = parser.parse_args()

    if args.mtdna_seq_path is None:
        logger.info("Using default mtDNA sequence from resources/mito.txt")
        try:
            mtdna_seq = files('mitoedit.resources').joinpath('mito.txt').read_text().replace("\n", "")
        except FileNotFoundError:
            logger.error("Default mtDNA sequence file not found in resources/mito.txt")
            raise
    else:
        logger.info(f"Reading mtDNA sequence from file: {args.mtdna_seq_path}")
        with open(args.mtdna_seq_path, "r") as fh:
            mtdna_seq = fh.read().replace("\n", "")

    bystander_df = None
    if args.bystander_file:
        bystander_file = os.path.abspath(args.bystander_file)
        if os.path.isfile(bystander_file):
            logger.info(f"Loading bystander data from {bystander_file}")
            bystander_df = pd.read_excel(bystander_file)
        else:
            logger.warning(f"Bystander file {bystander_file} does not exist. Skipping bystander information.")

    tale_nt_params = {
        'min_spacer': args.min_spacer,
        'max_spacer': args.max_spacer,
        'array_min': args.array_min,
        'array_max': args.array_max,
        'filter': args.filter,
        'cut_pos': args.cut_pos
    }

    results = process_mitoedit(mtdna_seq=mtdna_seq,
                               position=args.position,
                               mutant_base=args.mutant_base,
                               bystander_df=bystander_df,
                               tale_nt_params=tale_nt_params)

    if results['windows_df'].empty:
        logger.warning("No results generated. Exiting.")
        return

    os.makedirs(args.output_prefix, exist_ok=True)
    logger.info(f"Output directory created/verified: {args.output_prefix}")

    logger.info("Writing adjacent bases to FASTA file.")
    fasta_file = f'{args.output_prefix}/adjacent_bases.fasta'
    with open(fasta_file, 'w') as file:
        file.write(results['fasta_content'])
        logger.info(f"Finished writing FASTA file to {fasta_file}.")

    pipeline_windows_csv = f'{args.output_prefix}/pipeline_windows.csv'
    pipeline_bystanders_csv = f'{args.output_prefix}/pipeline_bystanders.csv'

    logger.info(f"Writing pipeline windows data to {pipeline_windows_csv}.")
    results['windows_df'].to_csv(pipeline_windows_csv, index=False)

    if not results['bystanders_df'].empty:
        logger.info(f"Writing pipeline bystanders data to {pipeline_bystanders_csv}.")
        results['bystanders_df'].to_csv(pipeline_bystanders_csv, index=False)
    else:
        logger.info("No bystanders information available to write.")

    combined_windows_csv = f'{args.output_prefix}/all_windows.csv'
    combined_bystanders_csv = f'{args.output_prefix}/all_bystanders.csv'

    logger.info(f"Saving combined windows data to {combined_windows_csv}.")
    results['windows_df'].to_csv(combined_windows_csv, index=False)

    if not results['bystanders_df'].empty:
        logger.info(f"Saving combined bystanders data to {combined_bystanders_csv}.")
        results['bystanders_df'].to_csv(combined_bystanders_csv, index=False)
    else:
        logger.info("No combined bystanders data to save.")

    if not results['talen_output_df'].empty:
        talen_output_path = f'{args.output_prefix}/talen_output.txt'
        logger.info(f"Saving TALE-NT output to {talen_output_path}.")
        results['talen_output_df'].to_csv(talen_output_path, sep='\t', index=False)

    logger.info("MitoEdit processing completed successfully.")