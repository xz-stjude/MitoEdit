import argparse
import random
import os
import pandas as pd
from ..logging import logger


class ChosTALEDsPipeline:
    """Pipeline for Cho sTALEDs base editing approach (T→C and A→G edits)"""
    
    def __init__(self):
        self.pipeline_name = "Cho_G1397_sTALEDs"

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
        """Get the complement of a DNA sequence"""
        logger.debug("Generating complement.")
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '}' : '{', '{' : '}', '[':']', ']':'['}
        complement_sequence = ''.join([complement_dict.get(base, random.choice(['A', 'T', 'C', 'G']) if base == 'N' else base) for base in sequence])
        return complement_sequence

    def _remove_whitespace(self, sequence):
        """Remove whitespace from the sequence"""
        logger.debug("Removing whitespace from the sequence.")
        return sequence.replace(" ", "").replace("\t", "").replace("\n", "")

    def _capitalize(self, sequence):
        """Capitalize the sequence"""
        logger.debug("Capitalizing the sequence.")
        return sequence.upper()

    def _find_consecutive_CT_sequences(self, sequence):
        """Find all positions where 'CT' sequences occur"""
        logger.debug("Finding consecutive CT sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'CT':
                positions.append(i + 1)  # 1-based indexing for C position
        return positions

    def _find_consecutive_GT_sequences(self, sequence):
        """Find all positions where 'GT' sequences occur"""
        logger.debug("Finding consecutive GT sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'GT':
                positions.append(i + 1)  # 1-based indexing for G position
        return positions

    def _find_consecutive_TG_sequences(self, sequence):
        """Find all positions where 'TG' sequences occur"""
        logger.debug("Finding consecutive TG sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'TG':
                positions.append(i + 2)  # 1-based indexing for G position
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

    def _find_consecutive_CA_sequences(self, sequence):
        """Find all positions where 'CA' sequences occur"""
        logger.debug("Finding consecutive CA sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'CA':
                positions.append(i + 1)  # 1-based indexing for C position
        return positions

    def _find_consecutive_GA_sequences(self, sequence):
        """Find all positions where 'GA' sequences occur"""
        logger.debug("Finding consecutive GA sequences.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'GA':
                positions.append(i + 1)  # 1-based indexing for G position
        return positions

    def _find_N_positions(self, sequence, start_position, pattern):
        """Find positions of a specific pattern in the sequence"""
        logger.debug("Finding positions of pattern %s in sequence.", pattern)
        positions = []
        for i in range(len(sequence) - len(pattern) + 1):
            if sequence[i:i+len(pattern)] == pattern:
                positions.append(start_position + i)
        return positions

    def _generate_windows_from_right_sTALED(self, circular_seq, pos, window_size):
        """Generate windows from right sTALED"""
        logger.debug("Generating windows from right sTALED.")
        windows = []
        for i in range(5, window_size - 4):
            start_index = pos - i
            end_index = pos + (window_size - i)
            window = circular_seq[start_index:end_index]
            windows.append(window)
        return windows

    def _generate_windows_from_left_sTALED(self, circular_seq, pos, window_size):
        """Generate windows from left sTALED"""
        logger.debug("Generating windows from left sTALED.")
        windows = []
        for i in range(5, window_size - 4):
            start_index = pos - (window_size - i)
            end_index = pos + i
            window = circular_seq[start_index:end_index]
            windows.append(window)
        return windows

    def process_mtDNA(self, mtDNA_seq, pos):
        """Main function which processes the DNA"""
        logger.info("Processing mtDNA sequence for position %d.", pos)

        nospace_mtDNA = self._capitalize(self._remove_whitespace(mtDNA_seq))
        T_positions = []
        A_positions = []
        # list to store positions where 'T' is present in the mtDNA sequence
        all_T_positions = self._find_consecutive_CT_sequences(nospace_mtDNA) + self._find_consecutive_GT_sequences(nospace_mtDNA) + self._find_consecutive_TG_sequences(nospace_mtDNA) + self._find_consecutive_TC_sequences(nospace_mtDNA)
        [T_positions.append(x) for x in all_T_positions if x not in T_positions]
        # list to store positions where 'A' is present in the mtDNA sequence
        all_A_positions = self._find_consecutive_AC_sequences(nospace_mtDNA) + self._find_consecutive_AG_sequences(nospace_mtDNA) + self._find_consecutive_CA_sequences(nospace_mtDNA) + self._find_consecutive_GA_sequences(nospace_mtDNA)
        [A_positions.append(x) for x in all_A_positions if x not in A_positions]
        dummy = 0

        if pos in T_positions:
            ref, mut, all_windows,dum = 'T', 'C',[],[]
            circular_seq = nospace_mtDNA + nospace_mtDNA
            start_index = pos - (16 + 15)
            end_index = pos + (15 + 15)
            adjacent_bases = circular_seq[start_index:end_index]
            marked_adjacent = self._mark_bases(adjacent_bases, 31, self._find_consecutive_GA_sequences(adjacent_bases) + self._find_consecutive_TC_sequences(adjacent_bases) + self._find_consecutive_AC_sequences(adjacent_bases) + self._find_consecutive_AG_sequences(adjacent_bases) + self._find_consecutive_CA_sequences(adjacent_bases) + self._find_consecutive_CT_sequences(adjacent_bases) + self._find_consecutive_GT_sequences(adjacent_bases) + self._find_consecutive_TG_sequences(adjacent_bases))
            #logger.info("The 60 adjacent bases to my target base are:", marked_adjacent)
            left_adjacent_bases = adjacent_bases[:30] #(base 0 to 29)
            right_adjacent_bases = adjacent_bases[31:] #(base 31 to 59)
            logger.info("The left and right adjacent bases are: %s and %s", left_adjacent_bases, right_adjacent_bases)
            FLAG=None
            #checking the position of the adjacent base --> if T present on either side --> then off-target
            if right_adjacent_bases[0] == 'T':
                dummy = 1
                dum.append(pos+1)
                FLAG=True
            if left_adjacent_bases[-1] =='T':
                dummy += 1
                dum.append(pos-1)
                FLAG=True
            if right_adjacent_bases[0] == 'A':
                FLAG=True
            if left_adjacent_bases[-1] =='A':
                FLAG=True
            for window_source in ["sTALED with AD on the right_TALE", "sTALED with AD on the left_TALE"]:
                for window_size in range(14, 19): #for window sizes of 14-18bp long --> MAJOR ASSUMPTION!!
                    if window_source == "sTALED with AD on the left_TALE":
                        sTALED_left_T_windows = self._generate_windows_from_right_sTALED(circular_seq, pos, window_size) #generating the TC windows(14bp-18bp)
                        TALES=False
                        for num, window in enumerate(sTALED_left_T_windows, start=5):
                            window_desc = f"Position {num} from the 5' end"
                            ws = f"{window_size}bp"
                            #off_target_sites_set = count_TC_sequences(window[4:13]) | count_GA_sequences(window[3:12]) | count_AC_sequences(window[4:13]) | count_AG_sequences(window[4:13]) | count_CA_sequences(window[3:12]) | count_CT_sequences(window[3:12]) | count_GT_sequences(window[3:12]) | count_TG_sequences(window[4:13]) #finds the off-target sites in each window
                            #if dummy:
                            #    off_target_sites_set.add(dummy)
                            marked_window = self._mark_bases(window, num, [(x + 3) for x in self._find_consecutive_GA_sequences(window[3:12])] + [(x + 4) for x in self._find_consecutive_AC_sequences(window[4:13])] + [(x + 4) for x in self._find_consecutive_AG_sequences(window[4:13])] + [(x + 3) for x in self._find_consecutive_CA_sequences(window[3:12])] + [(x + 3) for x in self._find_consecutive_CT_sequences(window[3:12])] + [(x + 3) for x in self._find_consecutive_GT_sequences(window[3:12])] + [(x + 4) for x in self._find_consecutive_TG_sequences(window[4:13])] + [(x + 4) for x in self._find_consecutive_TC_sequences(window[4:13])])
                            off_target_sites = marked_window.count('{') + dummy #because the original counting method was including the duplicate values in case of multiple contextss --> so now im counting on basis of number of '{' after marking the windows
                            int_pos_end = marked_window.find(']')
                            int_pos_ini = marked_window.find('[')
                            #to deal with T on either sides of the --> so it can be present in any context
                            if marked_window[int_pos_end + 1] == 'T':
                                final_window = self._mark_base_at_position(marked_window, int_pos_end + 1)
                            elif marked_window[int_pos_ini - 1] == 'T':
                                final_window = self._mark_base_at_position(marked_window, int_pos_ini - 1)
                            else:
                                final_window = marked_window
                            start_position = pos -num + 1
                            ac_positions = self._find_N_positions(window[4:13], start_position + 3 , 'AC')
                            tc_positions = self._find_N_positions(window[4:13], start_position + 3 , 'TC')
                            ag_positions = self._find_N_positions(window[4:13], start_position + 3 , 'AG')
                            tg_positions = self._find_N_positions(window[4:13], start_position + 3 , 'TG')
                            ca_positions = self._find_N_positions(window[3:12], start_position + 3, 'CA')
                            ct_positions = self._find_N_positions(window[3:12], start_position + 3, 'CT')
                            ga_positions = self._find_N_positions(window[3:12], start_position + 3, 'GA')
                            gt_positions = self._find_N_positions(window[3:12], start_position + 3, 'GT')
                            ftc = tc_positions + tg_positions + ct_positions + gt_positions
                            ftg = ac_positions + ag_positions + ca_positions + ga_positions
                            ftg = [x for x in ftg if x != pos]
                            ftc = [x for x in ftc if x != pos]
                            combined_set = set(ftc + ftg)
                            combined_set = set(ftc + ftg)  # Combine into a set
                            sorted_combined = sorted(combined_set) if combined_set else []  # Sort if not empty, else return [0]
                            all_windows.append((self.pipeline_name, window_source, pos, ref, mut, ws, final_window, window_desc, off_target_sites, sorted_combined, TALES, FLAG))
    
                    elif window_source == "sTALED with AD on the right_TALE":
                        sTALED_right_T_windows = self._generate_windows_from_left_sTALED(circular_seq, pos, window_size) #generating the TC windows(14bp-18bp)
                        TALES=False
                        for num, window in enumerate(sTALED_right_T_windows, start=5):
                            window_desc = f"Position {num} from the 3' end"
                            ws = f"{window_size}bp"
                            #off_target_sites_set = count_TC_sequences(window[-13:-4]) | count_GA_sequences(window[-13:-4]) | count_AC_sequences(window[-13:-4]) | count_AG_sequences(window[-13:-4]) | count_CA_sequences(window[-13:-4]) | count_CT_sequences(window[-13:-4]) | count_GT_sequences(window[-13:-4]) | count_TG_sequences(window[-13:-4]) #finds the off-target sites in each window
                            #if dummy:
                            #    off_target_sites_set.add(dummy)
                            #off_target_sites = len(off_target_sites_set)
                            marked_window = self._mark_bases(window, window_size - num + 1, [(x + window_size - 12) for x in self._find_consecutive_AG_sequences(window[-12:-3])] + [(x + window_size - 12) for x in self._find_consecutive_AC_sequences(window[-12:-3])] + [(x + window_size - 13) for x in self._find_consecutive_GA_sequences(window[-13:-4])] + [(x + window_size - 13) for x in self._find_consecutive_CA_sequences(window[-13:-4])] + [(x + window_size - 12) for x in self._find_consecutive_TC_sequences(window[-12:-3])] + [(x + window_size - 12) for x in self._find_consecutive_TG_sequences(window[-12:-3])] + [(x + window_size - 13) for x in self._find_consecutive_GT_sequences(window[-13:-4])] + [(x + window_size - 13) for x in self._find_consecutive_CT_sequences(window[-13:-4])])
                            off_target_sites = marked_window.count('{') + dummy
                            int_pos_end = marked_window.find(']')
                            int_pos_ini = marked_window.find('[')
                            #to deal with T on either sides of the --> so it can be present in any context
                            if marked_window[int_pos_end + 1] == 'T':
                                final_window = self._mark_base_at_position(marked_window, int_pos_end + 1)
                            elif marked_window[int_pos_ini - 1] == 'T':
                                final_window = self._mark_base_at_position(marked_window, int_pos_ini - 1)
                            else:
                                final_window = marked_window
                            start_position = pos - (window_size-num)
                            ac_positions = self._find_N_positions(window[-12:-3], start_position +window_size- 12 - 1, 'AC')
                            tc_positions = self._find_N_positions(window[-12:-3], start_position +window_size- 12 - 1, 'TC')
                            ag_positions = self._find_N_positions(window[-12:-3], start_position +window_size -12 - 1, 'AG')
                            tg_positions = self._find_N_positions(window[-12:-3], start_position +window_size -12 - 1, 'TG')
                            ca_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'CA')
                            ct_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'CT')
                            ga_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'GA')
                            gt_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'GT')
                            ftc = tc_positions + tg_positions + ct_positions + gt_positions
                            ftg = ac_positions + ag_positions + ca_positions + ga_positions
                            ftg = [x for x in ftg if x != pos]
                            ftc = [x for x in ftc if x != pos]
                            combined_set = set(ftc + ftg)
                            combined_set = set(ftc + ftg)  # Combine into a set
                            sorted_combined = sorted(combined_set) if combined_set else []  # Sort if not empty, else return [0]
                            all_windows.append((self.pipeline_name, window_source, pos, ref, mut, ws, final_window, window_desc, off_target_sites, sorted_combined, TALES, FLAG))

        elif pos in A_positions:
            logger.info("Base at position %d is in a editable context.", pos)
            ref, mut, all_windows, dum = 'A', 'G', [], []
            circular_seq = nospace_mtDNA + nospace_mtDNA
            start_index = pos - (16 + 15)
            end_index = pos + (15 + 15)
            adjacent_bases = circular_seq[start_index:end_index]
            marked_adjacent = self._mark_bases(adjacent_bases, 31, self._find_consecutive_GA_sequences(adjacent_bases) + self._find_consecutive_TC_sequences(adjacent_bases) + self._find_consecutive_AC_sequences(adjacent_bases) + self._find_consecutive_AG_sequences(adjacent_bases) + self._find_consecutive_CA_sequences(adjacent_bases) + self._find_consecutive_CT_sequences(adjacent_bases) + self._find_consecutive_GT_sequences(adjacent_bases) + self._find_consecutive_TG_sequences(adjacent_bases))
            left_adjacent_bases = adjacent_bases[:30] #(base 0 to 29)
            right_adjacent_bases = adjacent_bases[31:] #(base 31 to 59)
            logger.info("The left and right adjacent bases are: %s and %s", left_adjacent_bases, right_adjacent_bases)
            FLAG=None
            #checking the position of the adjacent base --> if T present on either side --> then off-target
            if right_adjacent_bases[0] == 'A':
                dummy = 1
                dum.append(pos+1)
                FLAG=True
            if left_adjacent_bases[-1] =='A':
                dummy += 1
                dum.append(pos-1)
                FLAG=True
            if right_adjacent_bases[0] == 'T':
                FLAG=True
            if left_adjacent_bases[-1] =='T':
                FLAG=True
            for window_source in ["sTALED with AD on the right_TALE", "sTALED with AD on the left_TALE"]: 
                for window_size in range(14, 19): #for window sizes of 14-18bp long
                    if window_source == "sTALED with AD on the right_TALE":
                        sTALED_left_A_windows = self._generate_windows_from_left_sTALED(circular_seq, pos, window_size) #generating the A windows(14bp-18bp)
                        TALES=False
                        for num, window in enumerate(sTALED_left_A_windows, start=5):
                            window_desc = f"Position {num} from the 3' end"
                            ws = f"{window_size}bp"
                            #off_target_sites = count_T_sequences(window[-13:-4]) + count_A_sequences(window[-13:-4])  #+dummy #finds the off-target sites in each window
                            marked_window = self._mark_bases(window, window_size - num + 1, [(x + window_size - 12) for x in self._find_consecutive_AG_sequences(window[-12:-3])] + [(x + window_size - 12) for x in self._find_consecutive_AC_sequences(window[-12:-3])] + [(x + window_size - 13) for x in self._find_consecutive_GA_sequences(window[-13:-4])] + [(x + window_size - 13) for x in self._find_consecutive_CA_sequences(window[-13:-4])] + [(x + window_size - 12) for x in self._find_consecutive_TC_sequences(window[-12:-3])] + [(x + window_size - 12) for x in self._find_consecutive_TG_sequences(window[-12:-3])] + [(x + window_size - 13) for x in self._find_consecutive_GT_sequences(window[-13:-4])] + [(x + window_size - 13) for x in self._find_consecutive_CT_sequences(window[-13:-4])])
                            off_target_sites = marked_window.count('{') + dummy
                            int_pos_end = marked_window.find(']')
                            int_pos_ini = marked_window.find('[')
                            #to deal with T on either sides of the --> so it can be present in any context
                            if marked_window[int_pos_end + 1] == 'A':
                                final_window = self._mark_base_at_position(marked_window, int_pos_end + 1)
                            elif marked_window[int_pos_ini - 1] == 'A':
                                final_window = self._mark_base_at_position(marked_window, int_pos_ini - 1)
                            else:
                                final_window = marked_window
                            start_position = pos - (window_size-num)
                            ac_positions = self._find_N_positions(window[-12:-3], start_position +window_size- 12 - 1, 'AC')
                            tc_positions = self._find_N_positions(window[-12:-3], start_position +window_size- 12 - 1, 'TC')
                            ag_positions = self._find_N_positions(window[-12:-3], start_position +window_size -12 - 1, 'AG')
                            tg_positions = self._find_N_positions(window[-12:-3], start_position +window_size -12 - 1, 'TG')
                            ca_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'CA')
                            ct_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'CT')
                            ga_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'GA')
                            gt_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'GT')
                            ftc = tc_positions + tg_positions + ct_positions + gt_positions
                            ftg = ac_positions + ag_positions + ca_positions + ga_positions
                            ftg = [x for x in ftg if x != pos]
                            ftc = [x for x in ftc if x != pos]
                            combined_set = set(ftc + ftg)  # Combine into a set
                            sorted_combined = sorted(combined_set) if combined_set else []  # Sort if not empty, else return [0]
                            all_windows.append((self.pipeline_name, window_source, pos, ref, mut, ws, final_window, window_desc, off_target_sites, sorted_combined, TALES, FLAG))

                    elif window_source == "sTALED with AD on the left_TALE":
                        sTALED_right_A_windows = self._generate_windows_from_right_sTALED(circular_seq, pos, window_size) #generating the A windows(14bp-18bp)
                        TALES=False
                        for num, window in enumerate(sTALED_right_A_windows, start=5):
                            window_desc = f"Position {num} from the 5' end"
                            ws = f"{window_size}bp"

                            marked_window = self._mark_bases(window, num, [(x + 3) for x in self._find_consecutive_GA_sequences(window[3:12])] + [(x + 4) for x in self._find_consecutive_AC_sequences(window[4:13])] + [(x + 4) for x in self._find_consecutive_AG_sequences(window[4:13])] + [(x + 3) for x in self._find_consecutive_CA_sequences(window[3:12])] + [(x + 3) for x in self._find_consecutive_CT_sequences(window[3:12])] + [(x + 3) for x in self._find_consecutive_GT_sequences(window[3:12])] + [(x + 4) for x in self._find_consecutive_TG_sequences(window[4:13])] + [(x + 4) for x in self._find_consecutive_TC_sequences(window[4:13])])
                            off_target_sites = marked_window.count('{') + dummy
                            int_pos_end = marked_window.find(']')
                            int_pos_ini = marked_window.find('[')
                            #to deal with T on either sides of the --> so it can be present in any context
                            if marked_window[int_pos_end + 1] == 'A':
                                final_window = self._mark_base_at_position(marked_window, int_pos_end + 1)
                            elif marked_window[int_pos_ini - 1] == 'A':
                                final_window = self._mark_base_at_position(marked_window, int_pos_ini - 1)
                            else:
                                final_window = marked_window
                            start_position = pos -num + 1
                            ac_positions = self._find_N_positions(window[4:13], start_position + 3 , 'AC')
                            tc_positions = self._find_N_positions(window[4:13], start_position + 3 , 'TC')
                            ag_positions = self._find_N_positions(window[4:13], start_position + 3 , 'AG')
                            tg_positions = self._find_N_positions(window[4:13], start_position + 3 , 'TG')
                            ca_positions = self._find_N_positions(window[3:12], start_position + 3, 'CA')
                            ct_positions = self._find_N_positions(window[3:12], start_position + 3, 'CT')
                            ga_positions = self._find_N_positions(window[3:12], start_position + 3, 'GA')
                            gt_positions = self._find_N_positions(window[3:12], start_position + 3, 'GT')
                            ftc = tc_positions + tg_positions + ct_positions + gt_positions
                            ftg = ac_positions + ag_positions + ca_positions + ga_positions
                            ftg = [x for x in ftg if x != pos]
                            ftc = [x for x in ftc if x != pos]
                            combined_set = set(ftc + ftg)  # Combine into a set
                            sorted_combined = sorted(combined_set) if combined_set else []  # Sort if not empty, else return [0]
                            all_windows.append((self.pipeline_name, window_source, pos, ref, mut, ws, final_window, window_desc, off_target_sites, sorted_combined, TALES, FLAG))


        else:
            logger.warning("Base at position %d is not in a editable context and cannot be edited by the %s pipeline.", pos, self.pipeline_name)
            print(f"Position {pos} is not in a editable context and cannot be edited by the {self.pipeline_name}.")
            return [], []  # Return empty lists to indicate failure

        return all_windows, adjacent_bases

    def process_bystander_data(self, all_windows, additional_file):
        """Process bystander information from additional file and return DataFrames."""
        logger.info("Appending additional bystanders information to the Excel file.")
        
        # Create a DataFrame from all_windows
        all_windows_df = pd.DataFrame(all_windows, columns=[
            'Pipeline', 'sTALED Type', 'Position', 'Reference Base', 'Mutant Base', 'Window Size',
            'Window Sequence', 'Target Location', 'Number of Bystanders',
            'Position of Bystanders', 'Optimal Flanking TALEs', 'Flag (CheckBystanderEffect)'
        ])

        # Only read and process the additional file if it is provided
        if additional_file:
            if not os.path.isfile(additional_file):
                logger.warning("The additional bystander file - %s does not exist. Skipping appending bystander information.", additional_file)
                new_data = pd.DataFrame()  # Create an empty DataFrame
            else:
                ftc_fga_positions = set(pos for _, _, _, _, _, _, _, _, _, positions, _, _ in all_windows for pos in positions)
                additional_df = pd.read_excel(additional_file)
                filtered_df = additional_df[additional_df['mtDNA_pos'].isin(ftc_fga_positions)]
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
    parser.add_argument('additional_file', type=str, help='Excel file containing additional bystander information') #this can be hardcoded ig
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

    pipeline = ChosTALEDsPipeline()

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

