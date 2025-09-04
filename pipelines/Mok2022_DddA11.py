from ..logging import logger
from .base_pipeline import BasePipeline


class Mok2022DddA11Pipeline(BasePipeline):
    """
    Mok2022 DddA11 Base Editing Pipeline (Multiple Context Support)
    
    Implements the Mok et al. 2022 DddA11 base editing approach for mitochondrial DNA editing.
    This advanced pipeline supports multiple sequence contexts for both cytosine-to-thymine (C→T)
    and guanine-to-adenine (G→A) conversions, providing the broadest editing capability.
    
    Supported Editing Contexts:
    - C→T edits in: TC, AC, CC contexts (5'-XC-3' where X = T, A, or C)
    - G→A edits in: GA, GT, GG contexts (5'-GX-3' where X = A, T, or G)
    
    Key Features:
    - Supports 6 different sequence contexts (most comprehensive)
    - Generates 14-20bp editing windows with flexible target positioning
    - Advanced bystander detection across all supported contexts
    - Optimized window positioning (positions 4 to window_size-4)
    - Comprehensive context analysis and reporting
    
    Unique Capabilities:
    - Multi-context editing in a single pipeline
    - Context-specific editing type reporting
    - Enhanced bystander effect analysis across all contexts
    - Most versatile pipeline for diverse editing requirements
    """
    
    def __init__(self):
        super().__init__()
        self.pipeline_name = "Mok2022_G1397_DddA11"

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
            return [], []


