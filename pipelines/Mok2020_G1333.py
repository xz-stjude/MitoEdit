from .mok2020_base_pipeline import Mok2020BasePipeline


class Mok2020G1333Pipeline(Mok2020BasePipeline):
    """
    Mok2020 G1333 Base Editing Pipeline (C→T and G→A edits)
    
    Implements the Mok et al. 2020 G1333 base editing approach for mitochondrial DNA editing.
    This pipeline targets cytosine-to-thymine (C→T) and guanine-to-adenine (G→A) conversions
    in specific sequence contexts (TC and GA respectively), with G1333-specific positioning.
    
    Key Features:
    - Targets C→T edits in TC contexts (5'-TC-3')
    - Targets G→A edits in GA contexts (5'-GA-3')
    - Uses G1333-specific window positioning (positions 3 to window_size-3)
    - Generates 14-20bp editing windows with flexible target positioning
    - Identifies and tracks bystander editing positions
    
    Positioning Strategy:
    - More aggressive positioning compared to G1397
    - Target can be positioned at positions 3 to window_size-3 within the editing window
    - Allows for broader target positioning flexibility
    - May have higher editing efficiency but potentially more off-target effects
    """
    
    def __init__(self):
        super().__init__()
        self.pipeline_name = "Mok2020_G1333"

    def _get_position_range(self, window_size):
        """Get the position range for G1333 window generation."""
        return range(3, window_size - 2)  # Different positioning logic for G1333
