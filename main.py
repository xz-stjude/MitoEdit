#!/usr/bin/env python

'''
Author - Devansh Shah
'''
import argparse
import os
import pandas as pd
import logging
import sys
import subprocess
import re
import random

# Setting up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', filename='logging_main.log')
logger = logging.getLogger(__name__)

# Constants for the TALE-NT Tool
MIN_SPACER = 14
MAX_SPACER = 18
ARR_MIN = 14
ARR_MAX = 18
CUT_POS = 31

def import_pipeline(pipeline_name):
    """Import the selected pipeline module from the pipelines folder."""
    try:
        module = __import__(f"pipelines.{pipeline_name}", fromlist=["process_mtDNA", "append_to_excel", "list_to_fasta"])
        process_mtDNA = module.process_mtDNA
        append_to_excel = module.append_to_excel
        list_to_fasta = module.list_to_fasta
    except ImportError as e:
        logger.error("Failed to import pipeline: %s", e)
        raise ValueError(f"Unknown pipeline: {pipeline_name}")

    return process_mtDNA, append_to_excel, list_to_fasta

def validate_input_files(input_file, additional_file):
    """Validate that the input files exist."""
    for file in [input_file, additional_file]:
        if not os.path.isfile(file):
            logger.error("The file - %s does not exist.", file)
            raise FileNotFoundError(f"The file - {file} does not exist.")

def run_findTAL(TALEN_dir, FASTA_fn, OUT_fn):
    """Run the findTAL script using subprocess."""
    findTAL_path = os.path.join(TALEN_dir, "findTAL.py")
    
    cmd = [
        "conda", "run", "-n", "run_talen_env",
        "python", findTAL_path,
        "--fasta", FASTA_fn,
        "--min", str(MIN_SPACER),
        "--max", str(MAX_SPACER),
        "--arraymin", str(ARR_MIN),
        "--arraymax", str(ARR_MAX),
        "--outpath", OUT_fn,
        "--filter", "1",
        "--filterbase", str(CUT_POS)
    ]
    
    logger.info("[cmd] %s", ' '.join(cmd))  # Log the command for debugging
    try:
        result = subprocess.call(cmd)
        if result != 0:
            raise RuntimeError(f"Command failed with return code {result}.")
    except Exception as e:
        logger.error("Error executing command: %s", e)
        sys.exit(1)

def process_talen_output(excel_file_path, txt_file_path, output_excel_file_path):
    """Process the output from the TALEN tool."""
    logger.info("Loading Excel file from %s", excel_file_path)
    
    try:
        excel_df = pd.read_excel(excel_file_path, sheet_name='All_Windows')
        second_sheet_df = pd.read_excel(excel_file_path, sheet_name='Bystanders_Info')  
        logger.info("Successfully loaded Excel file.")
    except Exception as e:
        logger.error("Failed to load Excel file %s: %s", excel_file_path, e)
        raise

    logger.info("Loading TXT file from %s", txt_file_path)
    
    try:
        txt_df = pd.read_csv(txt_file_path, delimiter='\t', skiprows=2)
        logger.info("Successfully loaded TXT file.")
    except pd.errors.ParserError:
        logger.error("Error reading the TXT file. Ensure it is tab-delimited and correctly formatted.")
        raise
    except Exception as e:
        logger.error("Failed to load TXT file %s: %s", txt_file_path, e)
        raise

    if 'Plus strand sequence' not in txt_df.columns:
        logger.error("Column 'Plus strand sequence' not found in the TALEN file.")
        raise ValueError("Column 'Plus strand sequence' not found in the TALEN file")

    plus_strand_sequences = txt_df['Plus strand sequence'].tolist()
    txt_bases = set(re.findall(r'[a-z]+', ' '.join(plus_strand_sequences)))

    logger.info("Comparing CSV file sequences with TALEN_Tool TXT file bases.")
    for index, row in excel_df.iterrows():
        cleaned_sequence = row['Window Sequence'].translate(str.maketrans('', '', '{}[]'))
        excel_df.at[index, 'Matching TALEs'] = cleaned_sequence.lower() in txt_bases
    
    logger.info("Saving updated DataFrame to output Excel file %s.", output_excel_file_path)
    with pd.ExcelWriter(output_excel_file_path, engine='openpyxl') as writer:
        excel_df.to_excel(writer, sheet_name='All_Windows', index=False)
        second_sheet_df.to_excel(writer, sheet_name='Bystanders_Info', index=False)
        combined_df = pd.concat([excel_df, second_sheet_df], ignore_index=True)
        combined_df.to_excel(writer, sheet_name='Combined_Information', index=False)
    
    logger.info("Updated final Excel file saved successfully.")

def setup_directories(parent_directory, pipeline_name):
    """Create necessary directories for output files."""
    directories = {
        'fasta': os.path.join(parent_directory, f'{pipeline_name}_fasta'),
        'all_windows': os.path.join(parent_directory, f'{pipeline_name}_all_windows'),
        'talen': os.path.join(parent_directory, f'{pipeline_name}_talent'),
        'final_output': os.path.join(parent_directory, f'{pipeline_name}_final_output')
    }

    for dir_path in directories.values():
        os.makedirs(dir_path, exist_ok=True)
    return directories

def main():
    parser = argparse.ArgumentParser(description='Process mtDNA sequence for base editing.')
    parser.add_argument('pipeline', type=str, choices=["Mok2020_G1397", "Mok2020_G1333", "Mok2022_DddA6", "Mok2022_DddA11", "Cho_sTALEDs"], help='Please specify the pipeline to use')
    #parser.add_argument('input_file', type=str, help='File containing the mtDNA sequence') #hard code?
    parser.add_argument('position', type=int, help='Position of the base to be changed')
    #parser.add_argument('additional_file', type=str, help='Excel file containing additional bystander information') #hard code?
    args = parser.parse_args()

    # Define file paths
    input_file = 'inputs/mito.txt'
    additional_file = 'inputs/annotated_mtDNA_10022024_for_bystanders.xlsx'

    # Import the selected pipeline functions
    process_mtDNA, append_to_excel, list_to_fasta = import_pipeline(args.pipeline)

    # Validate input files
    validate_input_files(input_file, additional_file)

    # Read the input DNA sequence
    logger.info("Reading the input DNA sequence %s.", input_file)
    with open(input_file, "r") as fh:
        mtDNA_seq = fh.read().replace("\n", "")

    while True:
        logging.info("Processing mtDNA sequence for position %d.", args.position)
        all_windows, adjacent_bases = process_mtDNA(mtDNA_seq, args.position)

        if not adjacent_bases:
            retry = input("Would you like to try a different position? (y/n): ").strip().lower()
            if retry == 'y':
                new_position = input("Enter a new position: ")
                try:
                    new_position = int(new_position)
                    args.position = new_position
                    continue
                except ValueError:
                    print("Invalid input. Please enter a valid integer.")
                    continue
            else:
                logging.info("Exiting the program.")
                return

        # Define paths for output files
        parent_directory = os.path.dirname(os.path.abspath(__file__))
        directories = setup_directories(parent_directory, args.pipeline)

        # Define output file paths
        fasta_file = os.path.join(directories['fasta'], f'{args.pipeline}_adjacent_bases_{args.position}.fasta')
        allw_file = os.path.join(directories['all_windows'], f'{args.pipeline}_all_windows_{args.position}.xlsx')

        if adjacent_bases:
            logging.info("Writing adjacent bases to FASTA file.")
            fasta_content = list_to_fasta(adjacent_bases, args.position)
            with open(fasta_file, 'w') as file:
                file.write(fasta_content)
                logging.info("Finished writing FASTA file to %s.", fasta_file)
        
        if all_windows:
            logging.info("Writing all windows and the bystander information to Excel file.")
            append_to_excel(all_windows, additional_file, allw_file)

        # Now call the TALEN tool
        #TALEN_dir = "C:\\Users\\dshah35\\OneDrive - St. Jude Children's Research Hospital\\Desktop\\files\\Github\\software\\talent_tools_master"
        TALEN_dir = "software/talent_tools_master"
        out_fn = os.path.join(directories['talen'], f'TALEN_output_{args.position}.txt')
        run_findTAL(TALEN_dir, fasta_file, out_fn)

        # Process the TALEN output
        output_excel_path = os.path.join(directories['final_output'], f'final_output_{args.position}.xlsx')
        process_talen_output(allw_file, out_fn, output_excel_path)
        break  # Exit the retry loop if processing was successful

if __name__ == "__main__":
    main()
