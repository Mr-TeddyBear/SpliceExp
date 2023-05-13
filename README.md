# SpliceExp
SpliceExp is a small colletion of scripts that tries to identify differently expressed exons that have splice region mutations. This is done by correlation mutations in splicing regions, with RNAseq data from the same tumor samples.
This software is developed using data on prostate adenocarcinoma patients, but could be employed on RNA count data for other types of cancer.


## Pre-requsites
This software is ment to use mutation annotation files (MAF), with one for each sample. SGSeq was used to process RNAseq data, create counts for all exons and junctions.
## Installation
Clone this repository and create a new Python virtual enviorment. <br>
This software was developed with Python 3.10, and all dependencies are listed in the `requerments.txt`. <br>
Install all dependencies using `python -m pip -r requierments.txt`. <br>

## Running

The src folder contains 2 scripts ment to be run from the commandline
### main.py
|Command | Description |
|----------------|--------------------|
|  -h, --help     |       show this help message and exit |
|  --MAF MAF       |      path folder containing MAF file for all samples and a map.txt file containing sample name and file name for sample. Sample name must correlate with sample name in count file. |
|  --features FEATURES |  path to csv file containing all feature counts as produced på SGSeq. |
|  --sep SEP          |  If a different seperation is used in features csv. Deafult \t |
|  --out OUT          |   Path to folder where outputs will be placed. Will create a new folder output at path given. |
|  --save             |   Arguemnt to write tables and figures to file. Deafult true. Will overwrite old fiures in outfolder. |
|  --plot_filetype PLOT_FILETYPE | Sets the filetype of exported plots. All formats accepted by matplotlib as accepted. |

<details>

<summary>Create heatmap</summary>

### create_heatmap.py
|Command | Description |
|----------------|--------------------|
|  -h, --help     |       show this help message and exit |
|  --MAF MAF       |      path folder containing MAF file for all samples and a map.txt file containing sample name and file name for sample. Sample name must correlate with sample name in count file. |
|  --features FEATURES |  path to csv file containing all feature counts as produced på SGSeq. |
|  --out OUT          |   Path to folder where outputs will be placed. Will create a new folder output at path given. |
|  --save             |   Arguemnt to write tables and figures to file. Deafult true. Will overwrite old fiures in outfolder. |
|  --plot_filetype PLOT_FILETYPE | Sets the filetype of exported plots. All formats accepted by matplotlib as accepted. |
| --geneID GENEID | The hugo symbol for the gene to create a heatmap. If none is provided the program will go in to a loop and await user input for a new Hugo symbol |
 | --sub_selection START END | Optional argument for defiing a genomic region inside the gene. The plot will only plot features inside this range.
  
  
</details>


## Author
This software was developed as a part of a master poject at the Oslo University hospital - Institute for Cancer Research - Dept. of Molecular Oncolgy and Oslo University - Institue of Informatics.
The software was written by Amund Isaksen.
