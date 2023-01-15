'''
This third script takes as input the list of structures from the RNA_Puzzles3 folder and gives them a pseudoenergy score.
It is based on linear interpolation of the scores from the objective function trained by training.py.
These linearly interpolated values are then summed to amount to the Gibbs pseudoenergy score.
'''
import json
import resources

# Variables
PDB_IDs_Puzzles = ["3_solution", "3_solution_0", "3_solution_1", "PZ3_solution_0", "PZ3_solution_0", "PZ3_Bujnicki_1", \
                   "PZ3_Bujnicki_2", "PZ3_Chen_1", "PZ3_Das_1", "PZ3_Das_2", "PZ3_Das_3", "PZ3_Das_4", "PZ3_Das_5", \
                   "PZ3_Dokholyan_1", "PZ3_Dokholyan_1", "PZ3_Major_1", "PZ3_Major_1"]

scoring_function = dict()

Parsed_PDBs = list()
MY_PDBS = resources.PDB_list_path(PDB_IDs_Puzzles)

def interpolate(available, value):
    with open("Results/Score_Output.json", 'r') as json_file:
        scoring_function = json.load(json_file)

if __name__ == '__main__':

