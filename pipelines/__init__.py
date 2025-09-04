from .base_pipeline import BasePipeline
from .mok2020_base_pipeline import Mok2020BasePipeline
from .Cho_sTALEDs import ChosTALEDsPipeline
from .Mok2020_G1397 import Mok2020G1397Pipeline
from .Mok2020_G1333 import Mok2020G1333Pipeline
from .Mok2022_DddA11 import Mok2022DddA11Pipeline

PIPELINE_CATALOG = {
    "Cho_sTALEDs": ChosTALEDsPipeline,
    "Mok2020_G1397": Mok2020G1397Pipeline,
    "Mok2020_G1333": Mok2020G1333Pipeline,
    "Mok2022_DddA11": Mok2022DddA11Pipeline
}

__all__ = [
    "BasePipeline",
    "Mok2020BasePipeline", 
    "ChosTALEDsPipeline",
    "Mok2020G1397Pipeline",
    "Mok2020G1333Pipeline",
    "Mok2022DddA11Pipeline",
    "PIPELINE_CATALOG"
]