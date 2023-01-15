'''
This second script takes as input the Gibbs energy pseudoscore for each couple (pair, bin), where pair
is a pairing between two nucleotides (A, C, T or U) and bin is a distance category of 1 Angstrom range between 0 and 20.
'''

# Variables
PDB_IDs_Puzzles = ["3_solution", "3_solution_0", "3_solution_1", "PZ3_solution_0", "PZ3_solution_0", "PZ3_Bujnicki_1", \
                   "PZ3_Bujnicki_2", "PZ3_Chen_1", "PZ3_Das_1", "PZ3_Das_2", "PZ3_Das_3", "PZ3_Das_4", "PZ3_Das_5", \
                   "PZ3_Dokholyan_1", "PZ3_Dokholyan_1", "PZ3_Major_1", "PZ3_Major_1"]

Parsed_PDBs = list()
MY_PDBS = resources.PDB_list_path(PDB_IDs_training)

from matplotlib.pyplot as plt

