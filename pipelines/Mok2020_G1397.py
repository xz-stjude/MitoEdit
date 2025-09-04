from .mok2020_base_pipeline import Mok2020BasePipeline


class Mok2020G1397Pipeline(Mok2020BasePipeline):
    """
    Mok2020 G1397 Base Editing Pipeline (C→T and G→A edits)
    
    Implements the Mok et al. 2020 G1397 base editing approach for mitochondrial DNA editing.
    This pipeline targets cytosine-to-thymine (C→T) and guanine-to-adenine (G→A) conversions
    in specific sequence contexts (TC and GA respectively).
    
    Key Features:
    - Targets C→T edits in TC contexts (5'-TC-3')
    - Targets G→A edits in GA contexts (5'-GA-3')
    - Uses G1397-specific window positioning (positions 4 to window_size-4)
    - Generates 14-20bp editing windows with flexible target positioning
    - Identifies and tracks bystander editing positions
    
    Positioning Strategy:
    - More conservative positioning compared to G1333
    - Target can be positioned at positions 4 to window_size-4 within the editing window
    - Optimized for reduced off-target effects
    """
    
    def __init__(self):
        super().__init__()
        self.pipeline_name = "Mok2020_G1397"

    def _get_position_range(self, window_size):
        """Get the position range for G1397 window generation."""
        return range(4, window_size - 3)  # Target can be at positions 4 to window_size-4
