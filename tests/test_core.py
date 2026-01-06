import pytest
import tempfile
import os
import sys
from pathlib import Path
import pandas as pd
from mitoedit import process_mitoedit


def test_process_mitoedit_default():
    results = process_mitoedit(  # Your DNA sequence
        position=11696,  # Target position
        mutant_base="A",  # Desired base
    )

    # Verify the function returns the expected structure
    print(
        results["talen_output_df"][
            ["TAL1 start", "TAL1 length", "TAL2 start", "TAL2 length", "Spacer range"]
        ]
    )
    assert isinstance(results, dict)
    assert "windows_df" in results
    assert "bystanders_df" in results
    assert "adjacent_bases" in results
    assert "fasta_content" in results
    assert "talen_output_df" in results
