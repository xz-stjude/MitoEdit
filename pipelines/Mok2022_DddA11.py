import argparse
import random 
import os 
import pandas as pd
from ..logging import logger


class Mok2022DddA11Pipeline:
    """Pipeline for Mok2022 DddA11 base editing approach (multiple contexts: TC, AC, CC, GA, GT, GG)"""
    
    def __init__(self):
        self.pipeline_name = "Mok2022_G1397_DddA11"

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
        logger.debug("Creating window from position %d.", pos)
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
        """Find GA contexts in the entire sequence"""
        logger.debug("Finding consecutive GA sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'GA':
                positions.append(i + 2)  # 1-based indexing for A position
        return positions

    def _find_consecutive_GT_sequences(self, sequence):
        """Find GT contexts in the entire sequence."""
        logger.debug("Finding consecutive GT sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'GT':
                positions.append(i + 2)  # 1-based indexing for T position
        return positions

    def _find_consecutive_GG_sequences(self, sequence):
        """Find GG contexts in the entire sequence."""
        logger.debug("Finding consecutive GG sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'GG':
                positions.append(i + 2)  # 1-based indexing for second G position
        return positions

    def _find_consecutive_TC_sequences(self, sequence):
        """Find TC contexts in the entire sequence"""
        logger.debug("Finding consecutive TC sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'TC':
                positions.append(i + 2)  # 1-based indexing for C position
        return positions

    def _find_consecutive_AC_sequences(self, sequence):
        """Find AC contexts in the entire sequence."""
        logger.debug("Finding consecutive AC sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'AC':
                positions.append(i + 2)  # 1-based indexing for C position
        return positions

    def _find_consecutive_CC_sequences(self, sequence):
        """Find CC contexts in the entire sequence."""
        logger.debug("Finding consecutive CC sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'CC':
                positions.append(i + 2)  # 1-based indexing for second C position
        return positions

    def _find_N_positions(self, window, start_position, context):
        """Find positions of a specific context in the window relative to the mtDNA sequence."""
        logger.debug("Finding %s contexts in the window.", context)
        positions = []
        start_index = 0
        while True:
            index_context = window.find(context, start_index)
            if index_context == -1:
                break
            if context in ['GA', 'GT', 'GG']:
                positions.append(start_position + index_context + 1)  # G position
            else:  # TC, AC, CC
                positions.append(start_position + index_context + 2)  # C/T position
            start_index = index_context + 1
        return positions

    def _process_context(self, nospace_mtDNA, pos, context_positions, ref, mut, context_name):
        """Process a specific context (TC, AC, CC, GA, GT, GG)"""
        all_windows = []
        
        # Create 60bp window around the target position
        start_index = pos - 31
        end_index = pos + 29
        adjacent_bases = self._create_window(nospace_mtDNA, pos, start_index, end_index)
        
        # Generate editing windows of different sizes
        for window_size in range(14, 21):  # 14-20bp windows
            for position_in_window in range(4, window_size - 3):  # Target can be at positions 4 to window_size-4
                window_start = pos - position_in_window
                window_end = pos + (window_size - position_in_window)
                window = self._create_window(nospace_mtDNA, pos, window_start, window_end)
                
                # Find bystander positions for all contexts
                bystander_positions = []
                for context in ['TC', 'AC', 'CC', 'GA', 'GT', 'GG']:
                    context_pos = self._find_N_positions(window, window_start, context)
                    bystander_positions.extend([p for p in context_pos if p != pos])
                
                bystander_positions = sorted(list(set(bystander_positions)))
                
                # Mark the window
                marked_window = self._mark_bases(window, position_in_window + 1, 
                                               [p - window_start + 1 for p in bystander_positions])
                
                all_windows.append((
                    self.pipeline_name, f"{context_name} editing", pos, ref, mut, f"{window_size}bp",
                    marked_window, f"Position {position_in_window + 1}", 
                    len(bystander_positions), bystander_positions, True, None
                ))
        
        return all_windows, adjacent_bases

    def process_mtDNA(self, mtDNA_seq, pos):
        """Main function which processes the DNA"""
        logger.info("Processing mtDNA sequence for position %d.", pos)
        
        nospace_mtDNA = self._capitalize(self._remove_whitespace(mtDNA_seq))
        
        # Find all context positions
        consecutive_TC_positions = self._find_consecutive_TC_sequences(nospace_mtDNA)
        consecutive_AC_positions = self._find_consecutive_AC_sequences(nospace_mtDNA)
        consecutive_CC_positions = self._find_consecutive_CC_sequences(nospace_mtDNA)
        consecutive_GA_positions = self._find_consecutive_GA_sequences(nospace_mtDNA)
        consecutive_GT_positions = self._find_consecutive_GT_sequences(nospace_mtDNA)
        consecutive_GG_positions = self._find_consecutive_GG_sequences(nospace_mtDNA)
        
        # Check which context the position belongs to and process accordingly
        if pos in consecutive_TC_positions:
            logger.info("Base at position %d is in a 5'-TC context.", pos)
            return self._process_context(nospace_mtDNA, pos, consecutive_TC_positions, 'C', 'T', "C→T (TC context)")
            
        elif pos in consecutive_AC_positions:
            logger.info("Base at position %d is in a 5'-AC context.", pos)
            return self._process_context(nospace_mtDNA, pos, consecutive_AC_positions, 'C', 'T', "C→T (AC context)")
            
        elif pos in consecutive_CC_positions:
            logger.info("Base at position %d is in a 5'-CC context.", pos)
            return self._process_context(nospace_mtDNA, pos, consecutive_CC_positions, 'C', 'T', "C→T (CC context)")
            
        elif pos in consecutive_GA_positions:
            logger.info("Base at position %d is in a 5'-GA context.", pos)
            return self._process_context(nospace_mtDNA, pos, consecutive_GA_positions, 'G', 'A', "G→A (GA context)")
            
        elif pos in consecutive_GT_positions:
            logger.info("Base at position %d is in a 5'-GT context.", pos)
            return self._process_context(nospace_mtDNA, pos, consecutive_GT_positions, 'G', 'A', "G→A (GT context)")
            
        elif pos in consecutive_GG_positions:
            logger.info("Base at position %d is in a 5'-GG context.", pos)
            return self._process_context(nospace_mtDNA, pos, consecutive_GG_positions, 'G', 'A', "G→A (GG context)")
            
        else:
            logger.warning("Base at position %d is not in an editable context and cannot be edited by the %s pipeline.", pos, self.pipeline_name)
            print(f"Position {pos} is not in an editable context and cannot be edited by the {self.pipeline_name}.")
            return [], []

    def process_bystander_data(self, all_windows, additional_file):
        """Process bystander information from additional file and return DataFrames."""
        logger.info("Appending additional bystanders information to the Excel file.")
        
        # Create a DataFrame from all_windows
        all_windows_df = pd.DataFrame(all_windows, columns=[
            'Pipeline', 'Editing Type', 'Position', 'Reference Base', 'Mutant Base', 'Window Size', 
            'Window Sequence', 'Target Location', 'Number of Bystanders', 
            'Position of Bystanders', 'Optimal Flanking TALEs', 'Flag (CheckBystanderEffect)'
        ])

        # Only read and process the additional file if it is provided
        if additional_file:
            if not os.path.isfile(additional_file):
                logger.warning("The additional bystander file - %s does not exist. Skipping appending bystander information.", additional_file)
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

    def list_to_fasta(self, dna_list, pos):
        """Converting the adjacent_bases into FASTA format"""
        logger.debug("converting to FASTA format")
        fasta_str = ""
        sequence = dna_list
        header = f">chrM_{pos}"
        fasta_str = f"{header}\n{sequence}\n"
        return fasta_str


def main():
    parser = argparse.ArgumentParser(description='Process mtDNA sequence for base editing.')
    parser.add_argument('input_file', type=str, help='File containing the mtDNA sequence')
    parser.add_argument('position', type=int, help='Position of the base to be changed (between 1 and 16569)')
    parser.add_argument('additional_file', type=str, help='Excel file containing additional bystander information')
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        logger.error("The file %s does not exist.", args.input_file)
        return
    
    if not os.path.isfile(args.additional_file):
        logger.error("The additional bystander file - %s does not exist.", args.additional_file)
        return
    
    logger.info("Reading the input sequence %s.", args.input_file)
    with open(args.input_file, "r") as fh:
        mtDNA_seq = fh.read().replace("\n", "")

    pipeline = Mok2022DddA11Pipeline()

    while True:  # Added retry loop
        logger.info("Processing mtDNA sequence for position %d.", args.position)
        all_windows, adjacent_bases = pipeline.process_mtDNA(mtDNA_seq, args.position)

        # Check if editing is possible
        if not adjacent_bases:
            retry = input("Would you like to try a different position? (y/n): ").strip().lower()
            if retry == 'y':
                new_position = input("Enter a new position (between 1 and 16569): ")
                try:
                    new_position = int(new_position)
                    if new_position < 1 or new_position > 16569:
                        print("Position must be between 1 and 16569.")
                        continue  # Prompt for a new position
                    args.position = new_position  # Update the position
                    continue  # Restart the loop with the new position
                except ValueError:
                    print("Invalid input. Please enter a valid integer.")
                    continue  # Prompt for a new position
            else:
                logger.info("Exiting the program.")
                return  # Exit the program

        # Define paths for output files
        parent_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        fasta_directory = os.path.join(parent_directory, f'{pipeline.pipeline_name}_fasta')
        csv_directory = os.path.join(parent_directory, f'{pipeline.pipeline_name}_all_windows')

        # Defining the files
        file_name = f'{pipeline.pipeline_name}_adjacent_bases_{args.position}.fasta'
        allw_name = f'{pipeline.pipeline_name}_all_windows_{args.position}.xlsx'

        file_path = os.path.join(fasta_directory, file_name)
        allw_path = os.path.join(csv_directory, allw_name)

        # Making directory if it doesn't exist
        os.makedirs(fasta_directory, exist_ok=True)
        os.makedirs(csv_directory, exist_ok=True)

        if adjacent_bases:
            logger.info("Writing adjacent bases to FASTA file.")
            fasta_content = pipeline.list_to_fasta(adjacent_bases, args.position)
            with open(file_path, 'w') as file:
                file.write(fasta_content)
            logger.info("FASTA file written to %s", file_path)

        if all_windows:
            logger.info("Processing bystander data and writing to Excel file.")
            all_windows_df, new_data = pipeline.process_bystander_data(all_windows, args.additional_file)
            
            with pd.ExcelWriter(allw_path, engine='openpyxl') as writer:
                all_windows_df.to_excel(writer, sheet_name='All Windows', index=False)
                if not new_data.empty:
                    new_data.to_excel(writer, sheet_name='Bystander Information', index=False)
            
            logger.info("Excel file written to %s", allw_path)

        break  # Exit the loop after successful processing
