from .base_pipeline import BasePipeline
from .Cho_sTALEDs import ChosTALEDsPipeline
from .Mok2020_unified import Mok2020UnifiedPipeline

PIPELINE_CATALOG = {
    "Cho_sTALEDs": ChosTALEDsPipeline,
    "Mok2020_Unified": Mok2020UnifiedPipeline,
}

__all__ = [
    "BasePipeline",
    "ChosTALEDsPipeline",
    "Mok2020UnifiedPipeline",
    "PIPELINE_CATALOG",
]
