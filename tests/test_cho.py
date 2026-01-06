import pytest
import tempfile
import os
import sys
from pathlib import Path
import pandas as pd


test_dir = Path(__file__).parent
mitoedit_dir = test_dir.parent
if str(mitoedit_dir) not in sys.path:
    sys.path.insert(0, str(mitoedit_dir))

from mitoedit.pipelines.Cho_sTALEDs import ChosTALEDsPipeline


class TestChosTALEDsPipeline:
    def test_run_find_bystanders_T(self):
        mtdna_seq = (
            "GGGGGTGGGGGGGGGGGGGGGGGGGGGGGGAGGTGGGGTGGGGGGGGGGGGGGGGGGGGGGAGGGGG"
        )
        position = 34
        all_windows, adjacent_bases = ChosTALEDsPipeline().process_mtDNA(
            mtdna_seq, position
        )
        print(all_windows)

    def test_run_find_bystanders_A(self):
        mtdna_seq = (
            "GGGGGTGGGGGGGGGGGGGGGGGGGGGGGTGGGAGGGGAGGGGGGGGGGGGGGGGGGGGGGAGGGGG"
        )
        position = 34
        all_windows, adjacent_bases = ChosTALEDsPipeline().process_mtDNA(
            mtdna_seq, position
        )
