#!/usr/bin/env python

import random 
import pandas as pd
from ..logging import logger

class UnifiedPipeline:
    def __init__(self, pipeline_type):
        self.pipeline_configs = {
            "Cho_sTALEDs": {
                "name": "Cho_G1397_sTALEDs",
                "contexts": {"T": ["CT", "GT", "TG", "TC"], "A": ["AC", "AG", "CA", "GA"]},
                "mutations": {"T": "C", "A": "G"},
                "window_generator": self._generate_staled_windows,
                "complex_marking": True
            },
            "Mok2020_G1397": {
                "name": "Mok2020_G1397",
                "contexts": {"C": ["TC"], "G": ["GA"]},
                "mutations": {"C": "T", "G": "A"},
                "window_generator": self._generate_simple_windows,
                "complex_marking": False
            },
            "Mok2020_G1333": {
                "name": "Mok2020_G1333",
                "contexts": {"C": ["TC"], "G": ["GA"]},
                "mutations": {"C": "T", "G": "A"},
                "window_generator": self._generate_g1333_windows,
                "complex_marking": False
            },
            "Mok2022_DddA11": {
                "name": "Mok2022_G1397_DddA11",
                "contexts": {"C": ["TC", "AC", "CC"], "G": ["GA", "GT", "GG"]},
                "mutations": {"C": "T", "G": "A"},
                "window_generator": self._generate_ddda11_windows,
                "complex_marking": False
            }
        }
        
        if pipeline_type not in self.pipeline_configs:
            raise ValueError(f"Unknown pipeline type: {pipeline_type}")
        
        self.config = self.pipeline_configs[pipeline_type]
        self.pipeline_name = self.config["name"]

    def mark_bases(self, sequence, target_position, off_target_positions):
        target_position -= 1
        off_target_positions = set(p - 1 for p in off_target_positions)
        marked_sequence = []
        for index, char in enumerate(sequence):
            if index == target_position:
                marked_sequence.append(f"[{char}]")
            elif index in off_target_positions:
                marked_sequence.append(f"{{{char}}}")
            else:
                marked_sequence.append(char)
        return ''.join(marked_sequence)

    def mark_base_at_position(self, sequence, target_position):
        marked_sequence = []
        for index, char in enumerate(sequence):
            if index == target_position:
                marked_sequence.append(f"{{{char}}}")
            else:
                marked_sequence.append(char)
        return ''.join(marked_sequence)

    def reverse_complement(self, sequence):
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '}' : '{', '{' : '}', '[':']', ']':'['}
        reverse_complement_sequence = ''.join([complement_dict.get(base, random.choice(['A', 'T', 'C', 'G']) if base == 'N' else base) for base in sequence])
        return reverse_complement_sequence[::-1]

    def complementing(self, sequence):
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        complement_sequence = ''.join([complement_dict.get(base, random.choice(['A', 'T', 'C', 'G']) if base == 'N' else base) for base in sequence])
        return complement_sequence

    def create_window(self, mtDNA_seq, pos, start_index, end_index):
        return mtDNA_seq[start_index:end_index]

    def remove_whitespace(self, seq):
        return ''.join(seq.split())

    def capitalize(self, seq):
        return seq.upper()

    def find_consecutive_sequences(self, sequence, context):
        positions = []
        start_index = 0
        while True:
            index = sequence.find(context, start_index)
            if index == -1:
                break
            if context in ['TC', 'GA', 'CA', 'CT', 'AC', 'CC']:
                positions.append(index + 2)
            else:
                positions.append(index + 1)
            start_index = index + 1
        return positions

    def count_sequences(self, sequence, context):
        count = 0
        start_index = 0
        while True:
            index = sequence.find(context, start_index)
            if index == -1:
                break
            count += 1
            start_index = index + 1
        return count

    def find_positions_in_window(self, window, start_position, context):
        positions = []
        start_index = 0
        while True:
            index = window.find(context, start_index)
            if index == -1:
                break
            positions.append(start_position + index + 1)
            start_index = index + 1
        return positions

    def list_to_fasta(self, dna_list, pos):
        sequence = dna_list
        header = f">chrM_{pos}"
        return f"{header}\n{sequence}\n"

    def _generate_staled_windows(self, mtDNA_seq, pos, window_size, base_type):
        windows = []
        if base_type == "T":
            for i in range(5, 14):
                if i != 13:
                    start_index = max(0, pos - window_size + i - 1)
                    end_index = min(len(mtDNA_seq), start_index + window_size)
                    if len(mtDNA_seq[start_index:end_index]) == window_size:
                        window = self.create_window(mtDNA_seq, pos, start_index, end_index)
                        windows.append(("left", window, i))
            
            for i in range(5, 14):
                if i != 13:
                    start_index = max(0, pos - i)
                    end_index = min(len(mtDNA_seq), start_index + window_size)
                    if len(mtDNA_seq[start_index:end_index]) == window_size:
                        window = self.create_window(mtDNA_seq, pos, start_index, end_index)
                        windows.append(("right", window, i))
        else:
            for i in range(5, 14):
                if i != 13:
                    start_index = max(0, pos - window_size + i - 1)
                    end_index = min(len(mtDNA_seq), start_index + window_size)
                    if len(mtDNA_seq[start_index:end_index]) == window_size:
                        window = self.create_window(mtDNA_seq, pos, start_index, end_index)
                        windows.append(("left", window, i))
            
            for i in range(5, 14):
                if i != 13:
                    start_index = max(0, pos - i)
                    end_index = min(len(mtDNA_seq), start_index + window_size)
                    if len(mtDNA_seq[start_index:end_index]) == window_size:
                        window = self.create_window(mtDNA_seq, pos, start_index, end_index)
                        windows.append(("right", window, i))
        return windows

    def _generate_simple_windows(self, mtDNA_seq, pos, window_size, base_type):
        windows = []
        for i in range(4, 9):
            if i != 8:
                if base_type == "C":
                    start_index = max(0, pos - window_size + i - 1)
                    end_index = min(len(mtDNA_seq), start_index + window_size)
                    if len(mtDNA_seq[start_index:end_index]) == window_size:
                        window = self.create_window(mtDNA_seq, pos, start_index, end_index)
                        windows.append(("TC", window, i))
                else:
                    start_index = max(0, pos - i)
                    end_index = min(len(mtDNA_seq), start_index + window_size)
                    if len(mtDNA_seq[start_index:end_index]) == window_size:
                        window = self.create_window(mtDNA_seq, pos, start_index, end_index)
                        windows.append(("GA", window[::-1], i))
        return windows

    def _generate_g1333_windows(self, mtDNA_seq, pos, window_size, base_type):
        windows = []
        for i in range(4, 12):
            if i != 11:
                if base_type == "C":
                    start_index = max(0, pos - i)
                    end_index = min(len(mtDNA_seq), start_index + window_size)
                    if len(mtDNA_seq[start_index:end_index]) == window_size:
                        window = self.create_window(mtDNA_seq, pos, start_index, end_index)
                        windows.append(("TC", window, i))
                else:
                    start_index = max(0, pos - window_size + i - 1)
                    end_index = min(len(mtDNA_seq), start_index + window_size)
                    if len(mtDNA_seq[start_index:end_index]) == window_size:
                        window = self.create_window(mtDNA_seq, pos, start_index, end_index)
                        windows.append(("GA", window[::-1], i))
        return windows

    def _generate_ddda11_windows(self, mtDNA_seq, pos, window_size, base_type):
        windows = []
        for i in range(4, 9):
            if i != 8:
                if base_type == "C":
                    start_index = max(0, pos - window_size + i - 1)
                    end_index = min(len(mtDNA_seq), start_index + window_size)
                    if len(mtDNA_seq[start_index:end_index]) == window_size:
                        window = self.create_window(mtDNA_seq, pos, start_index, end_index)
                        windows.append(("NC", window, i))
                else:
                    start_index = max(0, pos - i)
                    end_index = min(len(mtDNA_seq), start_index + window_size)
                    if len(mtDNA_seq[start_index:end_index]) == window_size:
                        window = self.create_window(mtDNA_seq, pos, start_index, end_index)
                        windows.append(("GN", window[::-1], i))
        return windows

    def process_mtDNA(self, mtDNA_seq, pos):
        logger.info("Processing mtDNA sequence for position %d using %s pipeline.", pos, self.pipeline_name)
        
        nospace_mtDNA = self.capitalize(self.remove_whitespace(mtDNA_seq))
        
        all_positions = {}
        for base, contexts in self.config["contexts"].items():
            all_positions[base] = []
            for context in contexts:
                positions = self.find_consecutive_sequences(nospace_mtDNA, context)
                all_positions[base].extend(positions)
            all_positions[base] = list(set(all_positions[base]))

        for base, positions in all_positions.items():
            if pos in positions:
                return self._process_base_editing(nospace_mtDNA, pos, base)
        
        logger.warning("Base at position %d is not in an editable context for %s pipeline.", pos, self.pipeline_name)
        return [], []

    def _process_base_editing(self, nospace_mtDNA, pos, base_type):
        ref = base_type
        mut = self.config["mutations"][base_type]
        all_windows = []
        
        circular_seq = nospace_mtDNA + nospace_mtDNA
        start_index = pos - 31
        end_index = pos + 30
        adjacent_bases = circular_seq[start_index:end_index]
        
        left_adjacent_bases = adjacent_bases[:30]
        right_adjacent_bases = adjacent_bases[31:]
        logger.info("The left and right adjacent bases are: %s and %s", left_adjacent_bases, right_adjacent_bases)
        
        dummy = 0
        dum = []
        FLAG = None
        
        if right_adjacent_bases[0] == base_type:
            dummy += 1
            dum.append(pos + 1)
            FLAG = True
        if left_adjacent_bases[-1] == base_type:
            dummy += 1
            dum.append(pos - 1)
            FLAG = True

        for window_size in range(14, 19):
            windows = self.config["window_generator"](circular_seq, pos, window_size, base_type)
            
            for window_type, window, num in windows:
                window_desc = self._get_window_description(window_type, num)
                ws = f"{window_size}bp"
                
                marked_window, off_target_sites, sorted_positions = self._process_window(
                    window, pos, num, window_size, window_type, dummy, dum
                )
                
                all_windows.append((
                    self.pipeline_name, window_type, pos, ref, mut, ws, 
                    marked_window, window_desc, off_target_sites, 
                    sorted_positions, False, FLAG
                ))
        
        return all_windows, adjacent_bases

    def _get_window_description(self, window_type, num):
        if "left" in window_type or "sTALED with AD on the left" in window_type:
            return f"Position {num} from the 5' end"
        elif "right" in window_type or "sTALED with AD on the right" in window_type:
            return f"Position {num} from the 3' end"
        else:
            return f"Position {num} from the 3' end"

    def _process_window(self, window, pos, num, window_size, window_type, dummy, dum):
        contexts = []
        for base_contexts in self.config["contexts"].values():
            contexts.extend(base_contexts)
        
        off_target_positions = []
        for context in contexts:
            positions = self.find_consecutive_sequences(window, context)
            off_target_positions.extend(positions)
        
        if window_type == "sTALED":
            target_pos = num if "left" in window_type else window_size - num + 1
        else:
            target_pos = window_size - num + 1 if window_type in ["TC", "NC"] else num
        
        marked_window = self.mark_bases(window, target_pos, off_target_positions)
        off_target_sites = marked_window.count('{') + dummy
        
        all_positions = []
        for context in contexts:
            positions = self.find_positions_in_window(window, pos - num + 1, context)
            all_positions.extend(positions)
        
        all_positions = [p for p in all_positions if p != pos]
        all_positions.extend(dum)
        sorted_positions = sorted(set(all_positions))
        
        return marked_window, off_target_sites, sorted_positions

    def process_bystander_data(self, all_windows, additional_data=None):
        logger.info("Processing bystander information.")
        
        columns = ['Pipeline', 'Type', 'Position', 'Reference Base', 'Mutant Base', 'Window Size', 
                  'Window Sequence', 'Target Location', 'Number of Bystanders', 
                  'Position of Bystanders', 'Optimal Flanking TALEs', 'Flag (CheckBystanderEffect)']
        
        all_windows_df = pd.DataFrame(all_windows, columns=columns)

        if additional_data is not None:
            ftc_fga_positions = set(pos for row in all_windows for pos in row[9])
            filtered_df = additional_data[additional_data['mtDNA_pos'].isin(ftc_fga_positions)]
            new_data = filtered_df[['mtDNA_pos', 'Ref. Allele', 'Mutant Allele', 
                                     'Location', 'Predicted Impact', 'Syn vs NonSyn', 
                                     'AA Variant', 'Func. Impact', 'MutationAssessor Score']]
            new_data.columns = ['Bystander Position', 'Reference Base', 'Mutant Base', 
                                'Location On Genome', 'Predicted Mutation Impact', 
                                'SNV Type', 'AA Variant', 'Functional Impact', 
                                'MutationAssessor Score']
        else:
            new_data = pd.DataFrame()

        logger.info("Successfully processed bystander information.")
        return all_windows_df, new_data

def process_pipeline(pipeline_type, mtDNA_seq, pos, additional_data=None):
    """Main entry point for processing any pipeline type."""
    pipeline = UnifiedPipeline(pipeline_type)
    all_windows, adjacent_bases = pipeline.process_mtDNA(mtDNA_seq, pos)
    
    if not adjacent_bases:
        return None, None, None
    
    all_windows_df, bystander_df = pipeline.process_bystander_data(all_windows, additional_data)
    fasta_content = pipeline.list_to_fasta(adjacent_bases, pos)
    
    return all_windows_df, bystander_df, fasta_content