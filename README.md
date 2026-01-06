# MitoEdit: a tool for targeted mitochondrial base editing

## About

**MitoEdit** is a novel Python workflow designed for targeted base editing of
the mitochondrial DNA (mtDNA). This tool streamlines the selection of optimal
targeting windows for base editing, helping researchers efficiently target
mtDNA mutations that are linked to various diseases. By enabling targeted base
editing of the mitochondrial genome, MitoEdit will speed up the study of
mtDNA-related diseases, help with preclinical drug testing, and enable
therapeutic approaches to correct pathogenic mutations.

![Rough workflow](imgs/MitoEdit_Pipeline.png)

## Overview

MitoEdit lets users input DNA sequences in text format, specify the target base
position, and indicate the desired modification. The tool processes this
information to identify candidate target windows, list the number and position
of potential bystander edits, and find optimal flanking TALE sequences, where
applicable. All the results are provided in a structured CSV format including
detailed logs to track progress.

### Pipelines

The current MitoEdit workflow includes three different pipelines based on the
editing patterns observed in the following base editor systems:

1. [Mok2020 Unified](https://rdcu.be/dXUIG) - Combines G1333, G1397, and DddA11 variants with all positioning strategies
2. [Cho_sTALEDs](https://www.sciencedirect.com/science/article/pii/S0092867422003890)

A detailed description for each pipeline can be found in the
[README](mitoedit/pipelines/README.md) file in the _pipelines_ folder.

**Note**: The unified Mok2020 pipeline automatically includes all variants (G1333, G1397, DddA11) and positioning strategies, eliminating the need for separate pipeline selection.

### TALE-NT Tool

The workflow uses the [TALE-NT tool](https://tale-nt.cac.cornell.edu/) to
identify optimal flanking TALE arrays around each candidate target window.

## Quick Start

If you have all the necessary tools installed, you can use these commands to
access MitoEdit:

#### For targeting the human mitochondrial DNA:

```
mitoedit <position> <mutant_base>
```

**Note:** The `position` is based on the [human mitochondrial
genome](https://www.ncbi.nlm.nih.gov/nuccore/251831106) sequence (NC_012920.1)

#### For targeting any other DNA sequence:

```
mitoedit --mtdna_seq_path <input_DNA_file> <position> <mutant_base>
```

## Prerequisites

- Python 3.8 or newer
- Required Python packages (automatically installed with pip):
  - `pandas>=2.0.0`
- Optional dependencies for web interface:
  - `fastapi>=0.100.0`
  - `uvicorn>=0.20.0`
  - `python-multipart>=0.0.5`
  - `jinja2>=3.0.0`

## Installation

#### 1. Clone the repository:

```
git clone https://github.com/Kundu-Lab/mitoedit.git
cd mitoedit
```

#### 2. Install the package:

```
pip install .
```

#### 3. For web interface support (optional):

```
pip install .[web]
```

#### 4. For development (optional):

```
pip install .[dev]
```

## Library Usage

MitoEdit can be used as a Python library in addition to the command-line interface:

```python
from mitoedit import process_mitoedit
import pandas as pd

# Process mitochondrial DNA editing
results = process_mitoedit(
    mtdna_seq="ATCGATCG...",  # Your DNA sequence
    position=100,             # Target position
    mutant_base="A",         # Desired base
    bystander_df=None,       # Optional bystander data
    tale_nt_params={         # TALE-NT parameters
        'min_spacer': 14,
        'max_spacer': 18,
        'array_min': 14,
        'array_max': 18,
    }
)

# Access results
windows_df = results['windows_df']
bystanders_df = results['bystanders_df']
fasta_content = results['fasta_content']
talen_output_df = results['talen_output_df']
```

## Web Interface

MitoEdit includes a web interface that makes it easy to analyze DNA sequences
through your browser. Here's how to set it up:

### System Requirements

- Python 3.8 or newer
- Web interface dependencies installed (`pip install .[web]`)

### Setup Instructions

1. **Install Dependencies**

```bash
pip install .[web]
```

2. **Run the Web Server**

```bash
# Set the required password for web interface access
export MITOEDIT_PASSWORD=your_secure_password

# Start the web server
python -m mitoedit.web.main
```

The web interface will be available at http://localhost:8000, where you can:

- Input DNA sequence position and bases
- Upload custom DNA sequence files
- View and download analysis results in CSV format

### Troubleshooting

If you encounter issues:

- Check the terminal output for error messages
- Ensure the `mitoedit` package is properly installed
- Verify all required files and directories exist

### Docker Installation

MitoEdit is also available as a Docker container, which provides an isolated environment with all dependencies pre-installed.

#### Prerequisites

- Docker installed on your system ([Install Docker](https://docs.docker.com/get-docker/))

#### Building and Running with Docker

1. **Build the Docker image:**

```bash
sudo docker build -t mitoedit .
```

2. **Run the container:**

For HTTP (port 80):

```bash
sudo docker run -d -p 80:80 -e MITOEDIT_PASSWORD=your_secure_password mitoedit
```

For HTTPS (port 443):

```bash
sudo docker run -d -p 443:443 -e PORT=443 -e MITOEDIT_PASSWORD=your_secure_password mitoedit
```

The web interface will be available at:

- HTTP: http://localhost:80 (when using port 80)
- HTTPS: https://localhost:443 (when using port 443)

#### Stopping/Terminating the Docker Container

To stop or terminate the running MitoEdit Docker container:

1. **List all running containers to find the container ID:**

```bash
sudo docker ps
```

2. **Stop the container:**

```bash
sudo docker stop <container_id>
```

For example: `sudo docker stop f11e4cc1615e`

3. **Remove the container (optional, if you want to completely remove it):**

```bash
sudo docker rm <container_id>
```

#### Docker Troubleshooting

- Ensure Docker is running on your system
- Check if port 80 is available (if not, you can map to a different port using `-p 8000:80` for example)
- For HTTPS, check if port 443 is available (if not, you can map to a different port using `-p 8443:443` for example)
- If you encounter permission issues, ensure you have the necessary privileges to run Docker commands

### Running as a Systemd Service

You can set up MitoEdit as a systemd service to run automatically on system boot:

1. **Create the service file:**

The repository includes a `mitoedit.service` file that defines how to run the Docker container as a systemd service.

2. **Install the service:**

Run the provided installation script with sudo:

```bash
sudo ./scripts/install_mitoedit_service.sh
```

3. **Managing the service:**

Once installed, you can manage the service using standard systemd commands:

```bash
sudo systemctl start mitoedit.service
sudo systemctl stop mitoedit.service
sudo systemctl restart mitoedit.service
sudo systemctl status mitoedit.service
```

4. **Accessing the web interface:**

When the service is running, access the MitoEdit web interface at:

```
https://localhost:443
```

## What input parameters does MitoEdit require?

MitoEdit requires the following parameters:

### Required Data Parameters:

#### 1. **Position**: The position of the base to target (1-based index).

#### 2. **Mutant Base**: The desired base conversion (A, T, C, or G).

### Optional Parameters:

#### Input File:
- `--mtdna_seq_path, -i`: The path to a file containing the DNA sequence as plain text.
- If not provided, MitoEdit will use the [human mtDNA sequence](https://www.ncbi.nlm.nih.gov/nuccore/251831106) from NCBI by default.

#### Bystander Analysis:
- `--bystander_file`: Excel file containing bystander effect annotations (optional, for human mtDNA analysis).

#### Output Configuration:
- `--output_prefix, -o`: Prefix for output CSV files (default: output).

#### TALE-NT Parameters:
- `--min_spacer`: Minimum spacer length for TALE-NT (default: 14).
- `--max_spacer`: Maximum spacer length for TALE-NT (default: 18).
- `--array_min`: Minimum array length for TALE-NT (default: 14).
- `--array_max`: Maximum array length for TALE-NT (default: 18).
- `--filter`: TALE-NT filter setting (default: 1).
- `--cut_pos`: TALE-NT cut position (default: 31).

## What does MitoEdit output?

MitoEdit generates the following outputs in the specified output directory:

### Output Files

#### CSV Files:

- `pipeline_windows.csv`: Lists the target windows generated from the pipeline.
- `pipeline_bystanders.csv`: Contains bystander effect information (if available).
- `all_windows.csv`: Contains all target windows from the pipeline.
- `all_bystanders.csv`: Contains all bystander effect information (if available).
- `talen_output.txt`: Contains the output from TALE-NT Tool describing the optimal flanking TALE sequences possible.

#### FASTA File:

- `adjacent_bases.fasta`: Contains the sequence adjacent to the target base, extending 30bp on each side.

## Usage

- For a full list of parameters, use the --help flag from the command line.

```
mitoedit --help
mitoedit -h
```

- To run the tool, use the following commands:

```
mitoedit <position> <mutant_base>
mitoedit --mtdna_seq_path <DNA.txt> <position> <mutant_base>
```

### Advanced Usage Examples

#### Basic usage with default settings:
```
mitoedit 11696 A
```

#### Custom output prefix:
```
mitoedit 11696 A --output_prefix my_results
```

#### Custom input file with bystander analysis:
```
mitoedit --mtdna_seq_path custom_dna.txt --bystander_file bystanders.xlsx 100 T
```

#### Custom TALE-NT parameters:
```
mitoedit 11696 A --min_spacer 12 --max_spacer 20 --filter 2
```

### Examples

#### To target the human mitochondrial DNA:

```
mitoedit 11696 A
```

**Expected Output:**
When you run this command, MitoEdit generates CSV files in the `output` directory. The main results are in `all_windows.csv` and `all_bystanders.csv` files.

**Note**: The [ ] represents the target base and { } represent bystander edits.

**1. all_windows.csv**
| Pipeline| Position |Reference Base | Mutant Base | Window Size | Window Sequence | Target Location | Number of bystanders | Position of Bystanders | Optimal Flanking TALEs | Flag (CheckBystanderEffect) |
|--------------|----------------------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|
|Mok2020_unified| 11696 |G| A |14bp| GCA[G]TCATT{C}TCAT| Position 4 from the 5' end |1| [11702] |FALSE | -
|Mok2020_unified | 11696 |G |A| 14bp| CGCA[G]TCATT{C}TCA| Position 5 from the 5' end| 1| [11702]| FALSE| -
|Mok2020_unified |11696|G |A| 14bp| GCGCA[G]T{C}ATTCTC| Position 6 from the 5' end| 1 |[11698] |FALSE | -
|Mok2020_unified | 11696 |G| A| 14bp |GGCGCA[G]T{C}ATTCT| Position 7 from the 5' end| 1| [11698]| FALSE | -

**Note:** If the column `Flag (CheckBystanderEffect)=TRUE`, you should manually check the results for potential amino acid changes caused by neighbouring bystanders on the same codon.

**2. all_bystanders.csv**
|Bystander Position| Reference Base| Mutant Base| Location On Genome| Predicted Mutation Impact |SNV Type| AA Variant| Functional Impact| MutationAssessor Score|
|---------|---------|---------| ---------| --------- |---------|---------|---------|---------|
|11698 |C |T |Complex 1| Predicted Benign |synonymous SNV|V313V| ||
|11702 |C |T |Complex 1| Predicted Pathogenic |nonsynonymous SNV|V315F|medium|3.44| ||
|11704 |C |T |Complex 1| Predicted Benign |synonymous SNV|V315L| ||

#### To target any other DNA sequence:

```
mitoedit --mtdna_seq_path test.txt 33 A
```

**Expected Output:**
When using an input file, the generated CSV files will contain results similar to the following: (example taken from the file provided in the [test file](test/input/test.txt) folder)

**1. all_windows.csv**
| Pipeline| Position |Reference_Base | Mutant Base | Window Size | Window Sequence | Target Location| Number of bystanders | Position of Bystanders | Optimal Flanking TALEs | Flag_CheckBystanderEffect |LeftTALE1 | RightTALE1|LeftTALE2|RightTALE2 |
|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------|---------|-------------|
|Mok2020_unified| 33| G| A| 14bp| TG{G}[G]A{G}AACT{C}TCT| Position 4 from the 5' end |3| [32, 35, 40] |FALSE| TRUE|
|Mok2020_unified| 33| G| A| 14bp| CTG{G}[G]A{G}AACTCTC| Position 5 from the 5' end| 2| [32, 35]| FALSE| TRUE |
|Mok2020_unified| 33| G| A| 14bp| ACTG{G}[G]AGAACTCT| Position 6 from the 5' end |1| [32] |FALSE| TRUE|
|Mok2020_unified| 33| G| A| 14bp| TACTG{G}[G]AGAACTC| Position 7 from the 5' end |1| [32] |TRUE |TRUE |T TACCCCCCACTATTAACC |TCTGTGCTAGTAACC A |T ACCCCCCACTATTAACC |TCTGTGCTAGTAACC A|

**Note**: When optimal flanking TALE sequences are found, the sequence is added to the `LeftTALE` and `RightTALE` columns respectively. The impact of bystander edits is not provided when using an input DNA file.

## Notes

- **Python 3 Compatibility**: MitoEdit has been fully updated to support Python 3.8+ for improved performance and modern library support.
- **Simplified Dependencies**: The tool now uses minimal dependencies with pandas as the core requirement.
- **CSV Output Format**: All results are now provided in CSV format for better compatibility and easier data processing.
- **Unified Pipeline**: The Mok2020 pipeline now combines all variants (G1333, G1397, DddA11) and positioning strategies in a single unified approach.
- **Package Installation**: Install using `pip install .` for proper package management and console script registration.
- **Library Usage**: MitoEdit can now be imported and used as a Python library in addition to the command-line interface.
- **Input File Formatting**: Ensure that your input file is correctly formatted as plain text, with the reference base matching the base at the specified position.
- **Supported File Formats**: MitoEdit can read plain text files for DNA sequences.
- **Minimum Upload Sequence Length**: The input file must contain at least 35 bases, covering the target base on either side, for accurate processing.
- **Output File Organization**: All output files are organized in the specified output directory with clear naming conventions.
- **Logging:** MitoEdit logs its progress and any issues encountered during execution to the console.
- **Species Support**: While the tool is designed primarily for human mtDNA, other DNA sequences can also be uploaded and used.
- **Modifying TALE-NT Workflow**: If no matching flanking TALE sequences are identified, consider modifying the TALE-NT parameter by setting `--filter 2`. This will identify all TALE pairs targeting any base in the target window, not just those for the target base. For further information, refer to the [TALE-NT FAQs](https://tale-nt.cac.cornell.edu/faqs).

## How to Cite?

If you use MitoEdit in your research, please cite the following paper:

- [MitoEdit: a pipeline for optimizing mtDNA base editing and predicting bystander effects](https://doi.org/10.1101/2025.01.22.634390)

## Contact

If you have any questions, please do not hesitate to contact me at:

- Xun Zhu: xun.zhu@stjude.org
- Devansh Shah: Devansh.Shah@stjude.org, shah.16@iitj.ac.in

## License

This project is licensed under the MIT License. See the `LICENSE` file for details. You are free to use and modify it, but please give credit to the original authors.

## Contributing

We welcome contributions! If you want to help improve MitoEdit, please fork the repository and submit a pull request.

## TALE-NT License

Copyright (c) 2011-2015, Nick Booher njbooher@gmail.com and Erin Doyle edoyle@iastate.edu.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE
