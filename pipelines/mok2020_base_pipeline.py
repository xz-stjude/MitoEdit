import logging
from multiprocessing import Value
logger = logging.getLogger(__name__)
from .base_pipeline import BasePipeline


class Mok2020BasePipeline(BasePipeline):
    """Base class for Mok2020 pipelines (G1397 and G1333) that handle C→T and G→A editing"""
    
    def __init__(self):
        super().__init__()
    
    def _find_C_positions(self, sequence):
        """Find all positions where 'C' occurs in TC contexts"""
        logger.debug("Finding C positions in TC contexts.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'TC':
                positions.append(i + 2)  # 1-based indexing for C position
        return positions

    def _find_G_positions(self, sequence):
        """Find all positions where 'G' occurs in GA contexts"""
        logger.debug("Finding G positions in GA contexts.")
        positions = []
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == 'GA':
                positions.append(i + 1)  # 1-based indexing for G position
        return positions

    def _find_GA_positions(self, window, start_position):
        """Find positions of GA contexts in the window relative to the mtDNA sequence."""
        logger.debug("Finding GA contexts in the window.")
        ga_positions = []
        start_index = 0
        while True:
            index_ga = window.find('GA', start_index)
            if index_ga == -1:
                break
            ga_positions.append(start_position + index_ga + 1)  # +1 for 1-indexed
            start_index = index_ga + 1
        return ga_positions

    def _find_TC_positions(self, window, start_position):
        """Find positions of TC contexts in the window relative to the mtDNA sequence."""
        logger.debug("Finding TC contexts in the window.")
        tc_positions = []
        start_index = 0
        while True:
            index_tc = window.find('TC', start_index)
            if index_tc == -1:
                break
            tc_positions.append(start_position + index_tc + 1+1)  # +1 for 1-indexed and additional +1 to get the position of C
            start_index = index_tc + 1
        return tc_positions

    def _get_position_range(self, window_size):
        """Get the position range for window generation. Override in subclasses."""
        raise NotImplementedError("Subclasses must implement _get_position_range")

    def _process_editing_context(self, nospace_mtDNA, pos, ref, mut, editing_type):
        """Process a specific editing context (C→T or G→A)"""
        all_windows = []
        
        # Create 60bp window around the target position
        start_index = pos - 31
        end_index = pos + 29
        adjacent_bases = self._create_window(nospace_mtDNA, pos, start_index, end_index)
        
        # Generate editing windows of different sizes
        for window_size in range(14, 21):  # 14-20bp windows
            position_range = self._get_position_range(window_size)
            for position_in_window in position_range:
                window_start = pos - position_in_window
                window_end = pos + (window_size - position_in_window)
                window = self._create_window(nospace_mtDNA, pos, window_start, window_end)
                
                # Find bystander positions
                ga_positions = self._find_GA_positions(window, window_start)
                tc_positions = self._find_TC_positions(window, window_start)
                bystander_positions = [p for p in ga_positions + tc_positions if p != pos]
                
                # Mark the window
                marked_window = self._mark_bases(window, position_in_window + 1, 
                                               [p - window_start + 1 for p in bystander_positions])
                
                all_windows.append((
                    self.pipeline_name, editing_type, pos, ref, mut, f"{window_size}bp",
                    marked_window, f"Position {position_in_window + 1}", 
                    len(bystander_positions), sorted(bystander_positions), True, None
                ))
        
        return all_windows, adjacent_bases

    def process_mtDNA(self, mtDNA_seq, pos):
        """Main function which processes the DNA sequence."""
        logger.info("Processing mtDNA sequence for position %d.", pos)
        
        nospace_mtDNA = self._capitalize(self._remove_whitespace(mtDNA_seq))
        C_positions = self._find_C_positions(nospace_mtDNA)
        logger.debug("C positions in TC contexts: %s", C_positions)
        G_positions = self._find_G_positions(nospace_mtDNA)
        logger.debug("G positions in GA contexts: %s", G_positions)
        
        if pos in C_positions:
            logger.info("Base at position %d is C in TC context and can be edited to T.", pos)
            return self._process_editing_context(nospace_mtDNA, pos, 'C', 'T', "C→T editing")
        
        elif pos in G_positions:
            logger.info("Base at position %d is G in GA context and can be edited to A.", pos)
            return self._process_editing_context(nospace_mtDNA, pos, 'G', 'A', "G→A editing")
        
        else:
            raise ValueError("Base at position %d is not in an editable context and cannot be edited by the %s pipeline.", pos, self.pipeline_name)