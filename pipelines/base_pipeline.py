import argparse
import random 
import os 
import pandas as pd
from abc import ABC, abstractmethod
import logging
logger = logging.getLogger(__name__)


class BasePipeline(ABC):
    """Abstract base class for all MitoEdit pipelines"""
    
    def __init__(self):
        self.pipeline_name = None  # To be set by subclasses
    
    @abstractmethod
    def process_mtDNA(self, mtDNA_seq, pos):
        """Main function which processes the DNA - must be implemented by subclasses
        
        Returns:
            tuple: (all_windows, adjacent_bases) where all_windows is a list of window data
                   and adjacent_bases is the sequence context around the target position
        """
        pass
    
    def process_bystander_data(self, all_windows, additional_file):
        """Process bystander information from additional file and return DataFrames."""
        logger.info("Appending additional bystanders information to the Excel file.")

        logger.debug("Number of columns in all_windows: %d", len(all_windows[0]) if all_windows else 0)
        
        # Create a DataFrame from all_windows
        all_windows_df = pd.DataFrame(all_windows, columns=[
            'Pipeline', 'Position', 'Reference Base', 'Mutant Base', 'Window Size', 
            'Window Sequence', 'Target Location', 'Number of Bystanders', 
            'Position of Bystanders', 'Optimal Flanking TALEs', 'Flag (CheckBystanderEffect)'
        ])

        # Only read and process the additional file if it is provided
        if additional_file:
            if not os.path.isfile(additional_file):
                logger.warning(f"The additional bystander file - {additional_file} does not exist. Skipping appending bystander information.")
                new_data = pd.DataFrame()  # Create an empty DataFrame
            else:
                bystander_positions = set(pos for _, _, _, _, _, _, _, _, _, positions, _, _ in all_windows for pos in positions)
                additional_df = pd.read_excel(additional_file)
                filtered_df = additional_df[additional_df['mtDNA_pos'].isin(bystander_positions)]
                new_data = filtered_df[['mtDNA_pos', 'Ref. Allele', 'Mutant Allele', 
                                         'Location', 'Predicted Impact', 'Syn vs NonSyn', 
                                         'AA Variant', 'Func. Impact', 'MutationAssessor Score']]
                new_data.columns = ['Bystander Position', 'Reference Base', 'Mutant Base', 
                                    'Location On Genome', 'Predicted Mutation Impact', 
                                    'SNV Type', 'AA Variant', 'Functional Impact', 
                                    'MutationAssessor Score']
        else:
            logger.warning("No additional file provided. Skipping bystander information.")
            new_data = pd.DataFrame()  # Create an empty DataFrame if no additional file is provided

        logger.info("Successfully processed bystander information, if available.")
        
        # Return the all_windows DataFrame for concatenation
        return all_windows_df, new_data

    def _mark_bases(self, sequence, target_position, off_target_positions):
        """Mark the target and bystander bases in the window"""
        logger.debug("Marking bases in the sequence.")
        target_position -= 1
        off_target_positions = set(p - 1 for p in off_target_positions)
        marked_sequence = []
        for index, char in enumerate(sequence):
            if index == target_position:
                marked_sequence.append(f"[{char}]")  # Target base with square brackets []
            elif index in off_target_positions:
                marked_sequence.append(f"{{{char}}}")  # Off-target base with curly braces {}
            else:
                marked_sequence.append(char)  # No special marking
        return ''.join(marked_sequence)

    def _mark_base_at_position(self, sequence, target_position):
        """Mark the base at the target position --> to mark the target base in the window"""
        logger.debug("Marking base at the specific position")
        marked_sequence = []
        for index, char in enumerate(sequence):
            if index == target_position:
                marked_sequence.append(f"{{{char}}}")  # Target base with curly brackets {}
            else:
                marked_sequence.append(char)  # No special marking
        return ''.join(marked_sequence)

    def _reverse_complement(self, sequence):
        """Get the reverse complement of a DNA sequence"""
        logger.debug("Generating reverse complement.")
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '}' : '{', '{' : '}', '[':']', ']':'['}
        reverse_complement_sequence = ''.join([complement_dict.get(base, random.choice(['A', 'T', 'C', 'G']) if base == 'N' else base) for base in sequence])
        return reverse_complement_sequence[::-1]

    def _complementing(self, sequence):
        """Get the complement of a DNA sequence."""
        logger.debug("Generating complement.")
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        complement_sequence = ''.join([complement_dict.get(base, random.choice(['A', 'T', 'C', 'G']) if base == 'N' else base) for base in sequence])
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

    def _find_consecutive_GA_sequences(self, sequence):
        """Find all positions where 'GA' sequences occur"""
        logger.debug("Finding consecutive GA sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'GA':
                positions.append(i + 1)  # 1-based indexing for G position
        return positions

    def _find_consecutive_GT_sequences(self, sequence):
        """Find all positions where 'GT' sequences occur"""
        logger.debug("Finding consecutive GT sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'GT':
                positions.append(i + 1)  # 1-based indexing for G position
        return positions

    def _find_consecutive_GG_sequences(self, sequence):
        """Find all positions where 'GG' sequences occur"""
        logger.debug("Finding consecutive GG sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'GG':
                positions.append(i + 1)  # 1-based indexing for first G position
        return positions

    def _find_consecutive_TC_sequences(self, sequence):
        """Find all positions where 'TC' sequences occur"""
        logger.debug("Finding consecutive TC sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'TC':
                positions.append(i + 2)  # 1-based indexing for C position
        return positions

    def _find_consecutive_AC_sequences(self, sequence):
        """Find all positions where 'AC' sequences occur"""
        logger.debug("Finding consecutive AC sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'AC':
                positions.append(i + 2)  # 1-based indexing for C position
        return positions

    def _find_consecutive_AG_sequences(self, sequence):
        """Find all positions where 'AG' sequences occur"""
        logger.debug("Finding consecutive AG sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'AG':
                positions.append(i + 2)  # 1-based indexing for G position
        return positions

    def _find_consecutive_CC_sequences(self, sequence):
        """Find all positions where 'CC' sequences occur"""
        logger.debug("Finding consecutive CC sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'CC':
                positions.append(i + 2)  # 1-based indexing for second C position
        return positions

    def _find_consecutive_TG_sequences(self, sequence):
        """Find all positions where 'TG' sequences occur"""
        logger.debug("Finding consecutive TG sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'TG':
                positions.append(i + 2)  # 1-based indexing for G position
        return positions

    def _find_consecutive_CT_sequences(self, sequence):
        """Find all positions where 'CT' sequences occur"""
        logger.debug("Finding consecutive CT sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'CT':
                positions.append(i + 1)  # 1-based indexing for C position
        return positions

    def _find_consecutive_CA_sequences(self, sequence):
        """Find all positions where 'CA' sequences occur"""
        logger.debug("Finding consecutive CA sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'CA':
                positions.append(i + 1)  # 1-based indexing for C position
        return positions

    def _find_N_positions(self, window, start_position, pattern):
        """Find positions of a specific pattern in the window"""
        logger.debug(f"Finding positions of pattern {pattern} in sequence.")
        positions = []
        for i in range(len(window) - len(pattern) + 1):
            if window[i:i+len(pattern)] == pattern:
                positions.append(start_position + i)
        return positions