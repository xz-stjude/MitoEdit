# MitoEdit: targeted mitochondrial base editing

## About
**MitoEdit** is a novel Python workflow designed for the targeted editing of mitochondrial DNA (mtDNA). This tool streamlines the selection of optimal spacing regions for base editing, enabling researchers to efficiently identify and correct pathogenic mtDNA mutations implicated in various diseases.

![Rough workflow](/imgs/rough_pipeline.png)

## Overview
MitoEdit allows users to input mtDNA sequences, specify target base positions, and indicate the desired base changes. The tool processes this information to generate a list of potential spacing regions, evaluates bystander edits, and outputs the results in a structured format including detailed logging for tracking progress.

### Pipelines
The current MitoEdit workflow uses four pipelines to predict the editing potential of a particular base. Please find more information related to each [here](paper link) (should I list pipelines?)

### TALE-NT Tool
The workflow utilizes the [TALE-NT tool](https://academic.oup.com/nar/article/40/W1/W117/1752530) to identify TALE sequences flanking the specified target position. This tool aids in predicting the effectiveness of the designed TALE proteins for the desired edits.

## Quick Start
- If you have all the dependencies and pre-requisites already installed, you can use the following commands to access the workflow - 

#### For targeting the human mitochondrial DNA:
```
python main.py <position> <reference_base> <mutant_base>
```
**Note:** The position is based on the [NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/251831106) sequence.

#### For targeting any other DNA sequence:
```
python main.py --input_file <input_DNA> <position> <reference_base> <mutant_base>
```

## Prerequisites
- Python 3.x
- Required Python packages:
  - `pandas` 
  - `openpyxl`
  - `argparse`
- Ensure that [Conda](https://docs.conda.io/en/latest/) is installed to manage environments and dependencies.
- The `findTAL.py` script from the TALE-NT tool should be available in the specified directory (software/talent_tools_master)

## Installation

#### 1. Clone the repository:
```
git clone <repository-url.git>
cd <repository-directory>
```
#### 2. Set Up the Conda Environment:
* To create and activate the conda environment, you can use the provided `environment.yml` file.
```
conda env create -f environment.yml
conda activate run_talen_env
```
- You can also make your conda environment using the following step.
```
conda create -n run_talen_env python=2.7.18 biopython=1.70
conda activate run_talen_env
```
**Note:** The conda environment should named as `run_talen_env` for the pipeline to correctly reference the TALE-NT tool.

#### 3. Install the additional required packages using pip:

```
pip install pandas openpyxl argparse
```

## What data parameters does MitoEdit require?
MitoEdit requires the following input parameters:
### Required data and data parameters: 

#### 1. Position: 
- The position of the base you want to change (1-based index). For example, if you want to modify the first base of the sequence, you would enter 1.

#### 2. Reference Base: 
- The original base at the specified position. This should be one of the following: A, T, C, or G.

#### 3. Mutant Base: 
- The base you want to change to. Again, this should be one of A, T, C, or G.

### Optional data and data parameters:
#### Input File:
- The path to a file (.txt or .fasta) containing the mtDNA sequence. 
- If not provided, MitoEdit will use a default file of the human mtDNA based on [NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/251831106)

## What MitoEdit outputs are included?
MitoEdit generates the following outputs to assist users in analyzing the base editing process:

### Output Directories
MitoEdit creates the following directories to organize the output files:
- fasta: Contains the FASTA files.
- pipeline_windows: Stores results from individual base-editing pipelines.
- all_windows: Contains the combined results of windows from all the pipelines.
- talen: Stores outputs from the TALE-NT tool.
- matching_output: Includes windows with matching TALE sequences.
- **final_output: Contains the final consolidated output files.**

**Note:** The user should look at the `final_output` directory to see the final output. The other directories are stored under the `running` directory.

### Output Files

#### Excel Files:

- `{pipeline}_{args.position}.xlsx`: Contains the spacing windows generated from each pipeline applicable.
- `all_windows_{position}.xlsx`: Contains all the potential spacing windows combined from all the pipelines.
- `matching_tales_{args.position}.xlsx`: Summary of whether matching TALEs are present for each spacing window. 
- **`final_{position}.xlsx`: Consolidated results including the matching TALE sequences.**

#### FASTA File:
- `adjacent_bases_{position}.fasta`: Contains sequences adjacent to the specified target base.

#### TALE-NT File:
- `TALENT_{args.position}.txt`: Contains the output from TALE-NT Tool.

## Usage
To execute the tool from the command line, use the following command structure:
```
python main.py <position> <reference_base> <mutant_base>
```
### Examples
#### To target a specific position in the human mitochondrial DNA:

```
python main.py 3243 A G
```
**Expected Output:**
When you run the command, the tool generates an Excel file named `final_{position}.xlsx`, which includes multiple sheets. Below are examples of the first five rows of two sheets: **All_Windows** and **Bystander_Effect**.

**All_Windows Sheet**
| Pipeline| sTALED types | Position |Reference_Base | Mutant Base | Window Size | Window Sequence | Description of window| Number of bystanders | Position of Bystanders | Matching TALEs | Flag_CheckBystanderEffect |
|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|
|Cho_G1397_sTALEDs|	sTALED with AD on the right_TALE|	3243	|A|	G	|14bp|	AG{A}{T}GGC{A}G[A]GCCC|	Position 5 from the 3' end	|3|	[3236, 3237, 3241]	|FALSE	|
|Cho_G1397_sTALEDs	|sTALED with AD on the right_TALE|	3243	|A	|G|	14bp|	GA{T}GGC{A}G[A]GCCCG|	Position 6 from the 3' end|	2|	[3237, 3241]|	FALSE|	
|Cho_G1397_sTALEDs|sTALED with AD on the right_TALE	|3243	|A	|G|	14bp|	ATGGC{A}G[A]GCCCGG|	Position 7 from the 3' end|	1	|[3241]	|FALSE	|
|Cho_G1397_sTALEDs	|sTALED with AD on the right_TALE|	3243	|A|	G|	14bp	|TGGC{A}G[A]GCCCGGT|	Position 8 from the 3' end|	1|	[3241]|	FALSE	|

**Bystanders_Effects Sheet**
|Bystander Position|	Reference Base|	Mutant Base|	Location On Genome|	Predicted Mutation Impact!	|SNV_Type|	AA_Variant|	Functional Impact|	MutationAssessor Score|
|---------|---------|---------|	---------|	---------	|---------|---------|---------|---------|
|3236	|A	|C	|tRNA|	Predicted Benign		|||	benign	||
|3236	|A	|G	|tRNA|	Predicted Benign		|||	benign	||
|3236	|A	|T	|tRNA|	Predicted Benign		|||	benign	||
|3237	|T	|A	|tRNA|	Predicted Benign		|||	benign	||

**Note:** The column `Flag_CheckBystanderEffect=TRUE` means you need to check the bystander change manually since a neighbouring base to the target base is a bystander base.

#### To target any other DNA sequence using an input file:
```
python main.py --input_file test.txt 33 G A
```

**Expected Output:**
When using an input file, the Excel file will contain an **All_Windows** sheet similar to the following:

**All_Windows Sheet**
| Pipeline| Position |Reference_Base | Mutant Base | Window Size | Window Sequence | Description of window| Number of bystanders | Position of Bystanders | Matching TALEs | Flag_CheckBystanderEffect |LeftTALE1 |	RightTALE1|LeftTALE2|RightTALE2 |
|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------|---------|-------------|
|Mok2020_G1397|	33|	G|	A|	14bp|	TG{G}[G]A{G}AACT{C}TCT|	Position 4 from the 5' end	|3|	[32, 35, 40]	|FALSE|	TRUE|				
|Mok2020_G1397|	33|	G|	A|	14bp|	CTG{G}[G]A{G}AACTCTC|	Position 5 from the 5' end|	2|	[32, 35]|	FALSE|	TRUE	|			
|Mok2020_G1397|	33|	G|	A|	14bp|	ACTG{G}[G]AGAACTCT|	Position 6 from the 5' end	|1|	[32]	|FALSE|	TRUE|				
|Mok2020_G1397|	33|	G|	A|	14bp|	TACTG{G}[G]AGAACTC|	Position 7 from the 5' end	|1|	[32]	|TRUE	|TRUE	|T TACCCCCCACTATTAACC	|TCTGTGCTAGTAACC A	|T ACCCCCCACTATTAACC	|TCTGTGCTAGTAACC A|

**Note:** 
- When a matching flanking TALE sequence is found, the sequence is appended to columns `LeftTALE` and `RightTALE`, respectively. Since you have uploaded an input file, impact of the bystander edits is not provided.

## Notes

- **Input File Formatting**: Ensure that your input file is correctly formatted, with the reference base matching the base at the specified position in the input file.

- **Minimum Upload Sequence Length**: Ensure that the input file contains at least 35 adjacent bases covering the target base on both sides for accurate processing.

- **Output File Name Conflicts**: Before running the tool, ensure there are no existing output files with the same name to prevent overwriting.

- **Logging:** The tool logs its progress and any issues encountered during execution in `logging_main.log.`

- **Species Support**: While the tool is designed for human mtDNA for bystander effect predictions, custom mtDNA sequences from other species can be uploaded and used. In fact even nuclear DNA can be used in some cases.

- **Modifying the Workflow**: If no matching flanking TALE sequences are identified, consider modifying the workflow by adjusting `FILTER = 2` This allows for identification of all TALE pairs targeting any base in the spacing region, rather than only those targeting a specific base. For further information, refer to the [TALE-NT FAQs](https://tale-nt.cac.cornell.edu/faqs).


## How to Cite?
If you use MitoEdit in your research, please cite it as follows:

- {paper link}


## Contact
If you have any questions, please do not hesitate to contact me at:
- Devansh Shah: Devansh.Shah@stjude.org, shah.16@iitj.ac.in


## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contributing
Contributions are welcome! Please open an issue or submit a pull request for any changes or improvements.

## TALE-NT License
Copyright (c) 2011-2015, Nick Booher njbooher@gmail.com and Erin Doyle edoyle@iastate.edu.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE
