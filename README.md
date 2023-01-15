# Creation of an objective function for the RNA folding problem
## Aim
The aim of this repository is to infer the Gibbs free energy of a pure (no interaction with another type macromolecule)
 RNA structure.

To do such a task, an objective function will be statistically trained on frequencies of distances of the C3 atom of nucleotides on a training set (n=10).

we consider the 10 different RNA pairings (AA, AC, AG, AU, CC, CG, CU, GG, GU, UU), whose distances have been classified into 
20 bins of 1 Angstrom (Å), from 0Å to 20Å (20Å excluded).
The rationale behind this upper boundary is that weak interatomic interactions fade away rapidly with distance and become negligible.

## Installation

Open a terminal and run the following commands:
```
git clone https://github.com/OGrondin/RNA_Folding
cd RNA_Folding
pip install -r requirements.txt 
```

## Scripts
This repository contains four scripts:
- resources.py, to be sourced by other scripts
- training.py, to compute frequencies
- plotting.py, to plot the results from training.py with score as a function of distance
- scoring.py, to score untrained structures using the RNA Puzzles Structures (Link below)

(https://github.com/RNA-Puzzles/raw_dataset_and_for_assessment)

## Training structures
10 training structures were chosen from the RNA-only PDB files.
Link: https://www.rcsb.org/stats/growth/growth-rna
Diversity and low number of chains (<3, for greater number of intrachain pairs) were criteria to choose these structures.

The chosen training structures have the following PDB IDs:
"2fqn", "4gxy", "5kpy", "5l4o", "5lyu", "5t83", "5u3g", "6wtl", "6ymc" and "7eem".

## References
### RNA Puzzles
Miao, Zhichao, et al. 
"RNA-Puzzles Round III: 3D RNA structure prediction of five riboswitches and one ribozyme." 
RNA 23.5 (2017): 655-672.

Miao, Zhichao, et al. 
"RNA-Puzzles Round II: assessment of RNA structure prediction programs applied to three large RNA structures." 
Rna (2015).

Cruz, José Almeida, et al. 
"RNA-Puzzles: a CASP-like evaluation of RNA three-dimensional structure prediction." 
Rna 18.4 (2012): 610-625.