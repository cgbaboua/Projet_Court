# BQ4CY010 Programmation 3 - Projet court | CLUSTAL
**Author** : GBABOUA Cassandra

## Summary 
This program performs multiple sequence alignment of protein sequences using a Python-based implementation of the CLUSTAL algorithm. It utilizes the Needleman-Wunsch algorithm for pairwise alignment and supports the BLOSUM62 substitution matrix. The project is part of the "BQ4CY010 Programmation 3" course at the Paris Cité University.

## System requirement
This program is compatible with macOS and UNIX-based operating systems. Windows users can run this project using Windows Subsystem for Linux (WSL).
This program require a Python 3.11.5 version.

## Installation

### Clone Repository

Clone the workflow github repository :

`git clone https://github.com/cgbaboua/Projet_Court.git`

### Repository contents 
[Projet_Court/](https://github.com/cgbaboua/Projet_Court)
  - [README.md](https://github.com/cgbaboua/Projet_Court/blob/main/README.md)
  - [blosum62.txt](https://github.com/cgbaboua/Projet_Court/blob/main/blosum62.txt)
  - [clustal.py](https://github.com/cgbaboua/Projet_Court/blob/main/clustal.py)
  - [clustal_env.yml](https://github.com/cgbaboua/Projet_Court/blob/main/clustal_env.yml)
  - [Examples/](https://github.com/cgbaboua/Projet_Court/tree/main/Examples)
    - [CoreG-Dendrogram_2.png](https://github.com/cgbaboua/Projet_Court/blob/main/Examples/CoreG-Dendrogram_2.png)
    - [CoreG-MultiAlignmt_2.txt](https://github.com/cgbaboua/Projet_Court/blob/main/Examples/CoreG-MultiAlignmt_2.txt)
    - [H-Globin-Dendrogram_1.png](https://github.com/cgbaboua/Projet_Court/blob/main/Examples/H-Globin-Dendrogram_1.png)
    - [HGlobin-MultiAlignmt_1.txt](https://github.com/cgbaboua/Projet_Court/blob/main/Examples/HGlobin-MultiAlignmt_1.txt)
    - [HGlobin.fasta](https://github.com/cgbaboua/Projet_Court/blob/main/Examples/HGlobin.fasta)
    - [SeqFastaCoreG.fasta](https://github.com/cgbaboua/Projet_Court/blob/main/Examples/SeqFastaCoreG.fasta)
  
### Environment Setup

To set up the environment run :

`conda env create -f clustal_env.yml`

`conda activate clustal_env`


## Getting started

### Basic Usage
To launch the program, run : 

`python3 clustal.py --seq_file your_sequence_file.fasta --blosum_file blosum62.txt [--display_dendrogram]`

The --display_dendrogram argument is optional, and will generate a dendrogram if specified.

### Inputs 

- Seq_file: A FASTA file containing sequences (example files are provided, but you can use your own).
- Blosum_file: A .txt file containing the BLOSUM62 substitution matrix (included with the program).

### Features and Functions

`read_sequences(fasta_file)`

  Reads all the sequences from a FASTA file and stores them in a dictionary.

`matrix(matrix_file)`

  Creates the BLOSUM62 matrix, which is stored in a dictionary.

`needleman_wunsch(seq1, seq2, blosum62, gap_penalty, return_score)`

  Implements the Needleman-Wunsch algorithm for pairwise sequence alignment.

`pairwise_alignment_scores(seq1_id, seq2_id, score, dict_alignment_scores)`

  Initializes a dictionary to store all the alignment scores.

`convert_score_to_distance(alignment_scores_dict, list_keys_sequences_dict)`

  Converts the alignment scores to a distance matrix.

`upgma(distance_matrix, list_keys_sequences_dict, display_dendrogram, dendrogram_file_path)`

  Performs UPGMA (Unweighted Pair Group Method with Arithmetic Mean) clustering.

`multiple_alignment(sequences, blosum62, Z_matrix, list_keys_sequences_dict)`

  Performs multiple sequence alignment using hierarchical clustering.

`needleman_wunsch_ali_to_seq(current_alignment, new_sequence, blosum62, gap_penalty)`

  Performs sequence alignment of a new sequence with an existing multiple sequence alignment.

`needleman_wunsch_ali_to_ali(current_alignment, new_alignment, blosum62, gap_penalty)`
  
  Performs alignment of a new alignment with an existing multiple sequence alignment.


### Packages Used 
- `os`: Used for interacting with the operating system to create directories, read files, etc.
- `sys`: Provides access to some variables used or maintained by the interpreter and to functions that interact strongly with the interpreter.
- `numpy (aliased as np)`: Utilized for numerical operations, especially for handling arrays and matrices.
- `matplotlib.pyplot (aliased as plt)`: Used for plotting and generating any graphs, like dendrograms.
- `scipy.cluster.hierarchy`: Provides functions for hierarchical clustering, including `dendrogram`, `linkage`, and `leaves_list`, which are used to create and display dendrograms.
- `argparse`: Command-line parsing library to identify what arguments need to be passed to the script, making it more user-friendly.

## Outputs and Examples

### Example

To test the program, use the example [HGlobin.fasta](https://github.com/cgbaboua/Projet_Court/blob/main/Examples/HGlobin.fasta) file :

`python3 clustal.py --seq_file Examples/HGlobin.fasta --blosum_file blosum62.txt --display_dendrogram`

*Note: Example FASTA sequences are sourced from Uniprot and NCBI. The program can handle FASTA formats from both Uniprot and NCBI.*
*Note : In the [Examples/](https://github.com/cgbaboua/Projet_Court/tree/main/Examples) repository, you'll find example fasta files and the results obtained with these files.*
### Outputs 

At the end of execution, a new folder named `Outputs/` will be created to contain the alignment results. This folder will include the multiple sequence alignment in a `.txt` format, and optional dendrogram images if requested.






