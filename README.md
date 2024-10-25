# MitoEdit: targeted mitochondrial base editing

## About
**MitoEdit** is a novel Python workflow designed for the targeted editing of mitochondrial DNA (mtDNA). This tool streamlines the selection of optimal spacing regions for base editing, helping researchers efficiently identify and correct pathogenic mtDNA mutations that are linked to various diseases. By enabling targeted base editing of the mitochondrial genome, MitoEdit will speed up the study of mtDNA-related diseases, help with preclinical drug testing, and allow for therapies that directly fix these mutations.

![Rough workflow](/imgs/rough_pipeline2.png)

## Overview
MitoEdit lets users input mtDNA sequences, specify target base positions, and indicate the desired base changes. The tool processes this information to create a list of potential spacing regions, evaluates bystander edits, and provides the results in a structured format including detailed logs to track progress.

### Pipelines
The current MitoEdit workflow uses four pipelines to predict how well a specific base can be edited. For more information about each pipeline click on the pipelines:
1. [Mok2020_G1333](https://rdcu.be/dXUIG)
2. [Mok2020_G1397](https://rdcu.be/dXUIG)
3. [Mok2022_G1397_DddA11](https://rdcu.be/dXUIm)
4. [Cho_sTALEDs](https://www.sciencedirect.com/science/article/pii/S0092867422003890)

**Note**: To use the DddA6 variant from the [Mok 2022 paper](https://rdcu.be/dXUIm) you can use the output from the `Mok2020_G1397` pipeline.

### TALE-NT Tool
The workflow uses the [TALE-NT tool](https://academic.oup.com/nar/article/40/W1/W117/1752530) to identify TALE sequences around the specified target position. This helps to design effective TALE proteins for the desired edits.

## Quick Start
If you have all the necessary tools installed, you can use these commands to access the workflow:

#### For targeting the human mitochondrial DNA:
```
python mitocraft.py <position> <reference_base> <mutant_base>
```
**Note:** The `position` is based on the [human mitochondrial genome](https://www.ncbi.nlm.nih.gov/nuccore/251831106) sequence (NC_012920.1)

#### For targeting any other DNA sequence:
```
python mitocraft.py --input_file <input_DNA> <position> <reference_base> <mutant_base>
```

## Prerequisites
- Python 3.x
- Required Python packages:
  - `pandas` 
  - `openpyxl`
  - `argparse`
- Download and install [Conda](https://docs.anaconda.com/anaconda/install/) if not already installed. Follow the prompts to complete installation.

## Installation

#### 1. Clone the repository:
```
git clone <repository-url.git>
cd mitoedit
```
#### 2. Set Up And Activate the Conda Environment:
* Use the provided `environment.yml` file To create and activate the conda environment:
```
conda env create -f environment.yml
```
- Alternatively, you can create the conda environment manually:
```
conda create -n run_talen_env python=2.7.18 biopython=1.70
```
**Note:** The conda environment should be named `run_talen_env` for the pipeline to correctly reference the TALE-NT tool. There is no need to activate the conda environment because the pipeline will automatically use `run_talen_env`. This environment is purely for the TALE-NT Tool.

## What data parameters does MitoEdit require?
MitoEdit needs the following input parameters:
### Required Data: 

#### 1. Position: 
- The position of the base you want to edit (1-based index). For example, if you want to modify the first base of the sequence, you would enter 1.

#### 2. Reference Base: 
- The original base at the specified position (A, T, C, or G).

#### 3. Mutant Base: 
- The desired base edit (A, T, C, or G).

### Optional Data Parameter:
#### Input File:
- The path to a file (.txt or .fasta) containing any DNA sequence. 
- If not provided, MitoEdit will use the [human mtDNA sequence](https://www.ncbi.nlm.nih.gov/nuccore/251831106) from NCBI by default.

## What MitoEdit outputs are included?
MitoEdit creates the following outputs to help you analyze the base editing process:

### Output Directories
MitoEdit organizes the output files in the following directories:
- fasta: Contains the FASTA files.
- pipeline_windows: Stores results from individual base-editing pipelines.
- all_windows: Contains combined results from all pipelines.
- talen: Stores outputs from the TALE-NT tool.
- matching_output: Includes windows with matching TALE sequences.
- **final_output: Contains the final consolidated output files.**

**Note:** Check the `final_output` directory to see the final results. The other directories are stored under the `running` directory.

### Output Files

#### Excel Files:

- `{pipeline}_{position}.xlsx`: Contains the spacing windows from each pipeline.
- `all_windows_{position}.xlsx`: Contains all potential spacing windows combined from all pipelines.
- `matching_tales_{position}.xlsx`: Summary of matching TALEs for each spacing window. 
- **`final_{position}.xlsx`: Consolidated results including the matching TALE sequences.**

#### FASTA File:
- `adjacent_bases_{position}.fasta`: Contains sequences adjacent to the target base.

#### TALE-NT File:
- `TALENT_{position}.txt`: Contains the output from TALE-NT Tool.

## Usage
- For a full list of parameters, use the --help flag.
```
python mitocraft.py --help
```
- To run the tool from the command line, use:
```
python mitocraft.py <position> <reference_base> <mutant_base>
```
### Examples
#### To target a specific position in the human mitochondrial DNA:

```
python mitocraft.py 3243 A G
```
**Expected Output:**
When you run this command, the tool generates an Excel file named `final_{position}.xlsx`, which includes two sheets. Below are examples of the first five rows from the sheets: **All_Windows** and **Bystander_Effect**.

**Note**: The [ ] represents the target base and { } represent bystander edits.

**1. All_Windows Sheet**
| Pipeline| sTALED types | Position |Reference_Base | Mutant Base | Window Size | Window Sequence | Target Location | Number of bystanders | Position of Bystanders | Matching TALEs | Flag_CheckBystanderEffect |
|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|
|Cho_G1397_sTALEDs|	sTALED with AD on the right_TALE|	3243	|A|	G	|14bp|	AG{A}{T}GGC{A}G[A]GCCC|	Position 5 from the 3' end	|3|	[3236, 3237, 3241]	|FALSE	|
|Cho_G1397_sTALEDs	|sTALED with AD on the right_TALE|	3243	|A	|G|	14bp|	GA{T}GGC{A}G[A]GCCCG|	Position 6 from the 3' end|	2|	[3237, 3241]|	FALSE|	
|Cho_G1397_sTALEDs|sTALED with AD on the right_TALE	|3243	|A	|G|	14bp|	ATGGC{A}G[A]GCCCGG|	Position 7 from the 3' end|	1	|[3241]	|FALSE	|
|Cho_G1397_sTALEDs	|sTALED with AD on the right_TALE|	3243	|A|	G|	14bp	|TGGC{A}G[A]GCCCGGT|	Position 8 from the 3' end|	1|	[3241]|	FALSE	|

**2. Bystanders_Effects Sheet**
|Bystander Position|	Reference Base|	Mutant Base|	Location On Genome|	Predicted Mutation Impact	|SNV_Type|	AA_Variant|	Functional Impact|	MutationAssessor Score|
|---------|---------|---------|	---------|	---------	|---------|---------|---------|---------|
|3236	|A	|C	|tRNA|	Predicted Benign		|||	benign	||
|3236	|A	|G	|tRNA|	Predicted Benign		|||	benign	||
|3236	|A	|T	|tRNA|	Predicted Benign		|||	benign	||
|3237	|T	|A	|tRNA|	Predicted Benign		|||	benign	||

**Note:** If the column `Flag_CheckBystanderEffect=TRUE` then you should manually check results for potential amino acid changes caused by neighbouring bystanders.

#### To target any other DNA sequence using an input file:
```
python mitocraft.py --input_file test.txt 33 G A
```

**Expected Output:**
When using an input file, the Excel file will contain an **All_Windows** sheet similar to the following:

**All_Windows Sheet**
| Pipeline| Position |Reference_Base | Mutant Base | Window Size | Window Sequence | Target Location| Number of bystanders | Position of Bystanders | Matching TALEs | Flag_CheckBystanderEffect |LeftTALE1 |	RightTALE1|LeftTALE2|RightTALE2 |
|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------|---------|-------------|
|Mok2020_G1397|	33|	G|	A|	14bp|	TG{G}[G]A{G}AACT{C}TCT|	Position 4 from the 5' end	|3|	[32, 35, 40]	|FALSE|	TRUE|				
|Mok2020_G1397|	33|	G|	A|	14bp|	CTG{G}[G]A{G}AACTCTC|	Position 5 from the 5' end|	2|	[32, 35]|	FALSE|	TRUE	|			
|Mok2020_G1397|	33|	G|	A|	14bp|	ACTG{G}[G]AGAACTCT|	Position 6 from the 5' end	|1|	[32]	|FALSE|	TRUE|				
|Mok2020_G1397|	33|	G|	A|	14bp|	TACTG{G}[G]AGAACTC|	Position 7 from the 5' end	|1|	[32]	|TRUE	|TRUE	|T TACCCCCCACTATTAACC	|TCTGTGCTAGTAACC A	|T ACCCCCCACTATTAACC	|TCTGTGCTAGTAACC A|

**Note**: When a matching flanking TALE sequence is found, the sequence is added to the `LeftTALE` and `RightTALE` columns. The impact of bystander edits is not provided when using an input file.

## Notes

- **Input File Formatting**: Ensure that your input file is correctly formatted, with the reference base matching the base at the specified position.
- **Supported File Formats**: MitoEdit can read .fasta and .txt formats DNA sequences.
- **Minimum Upload Sequence Length**: The input file must contain at least 35 bases covering the target base on both sides for accurate processing.
- **Output File Name Conflicts**: Check for existing output files with the same name before running the tool to prevent overwriting.
- **Logging:** The tool logs its progress and any issues encountered during execution in `logging_main.log.`
- **Species Support**: While the tool is designed for human mtDNA for bystander effect predictions, custom mtDNA sequences from other species can be uploaded and used. Nuclear DNA can also be used in some cases.
- **Modifying TALE-NT Workflow**: If no matching flanking TALE sequences are identified, consider modifying the workflow by setting `FILTER = 2`. This will identify all TALE pairs targeting any base in the spacing region, not just those targeting a specific base. For further information, refer to the [TALE-NT FAQs](https://tale-nt.cac.cornell.edu/faqs).


## How to Cite?
If you use MitoEdit in your research, please cite it as follows:

- {paper link}


## Contact
If you have any questions, please do not hesitate to contact me at:
- Devansh Shah: Devansh.Shah@stjude.org, shah.16@iitj.ac.in


## License
This project is licensed under the MIT License. See the `LICENSE` file for details. You are free to use and modify it, but please give credit to the original authors.

## Contributing
We welcome contributions! If you want to help improve MitoEdit, please fork the repository and submit a pull request.

## TALE-NT License
Copyright (c) 2011-2015, Nick Booher njbooher@gmail.com and Erin Doyle edoyle@iastate.edu.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE
