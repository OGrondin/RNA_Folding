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

scoring_dict = dict()
with open("Results/Score_Output.json", 'r') as json_file:
    scoring_dict = json.load(json_file)

Parsed_PDBs = list()
MY_PDBS = resources.PDB_list_path(PDB_IDs_Puzzles, "RNA_Puzzles3")
print(MY_PDBS)

def scoring(scoring_dict, list_PDB):
    with open("Results/Scoring.tsv", 'w') as _:
        pass
    id_index = 0
    for Parsed_PDB in list_PDB:
        distances = dict()
        for i in range(len(Parsed_PDB)):
            for j in range(i + 3, len(Parsed_PDB)):
                if Parsed_PDB[i][5] == Parsed_PDB[j][5]:  # Testing to ensure intra-chain comparison
                    pair, vector_dist = resources.prepare_distance(Parsed_PDB[i],Parsed_PDB[j])
                    eucl_dist = resources.dist_3D(vector_dist)
                    # print(pair, vector_dist, eucl_dist) # Debugging only
                    if eucl_dist > 20:
                        continue
                    if pair in distances:
                        new_dist = distances[pair]
                        new_dist.append(eucl_dist)
                        distances[pair] = new_dist
                    else:
                        new_dist = [eucl_dist]
                        distances[pair] = new_dist
        score = interpolate(scoring_dict, distances)
        output = "RNA Puzzles 3 ID: {}\t".format(PDB_IDs_Puzzles[id_index]) + str(round(score,4)) + "\n"
        with open("Results/Scoring.tsv", 'a+') as s:
            s.write(output)
        id_index += 1



def interpolate(scoring_dict, distances_puzzle):
    pair_index = 0
    interpol = [0 for _ in range(10)]
    for key in scoring_dict:
        dist_list = scoring_dict[key]
        if key in distances_puzzle:
            for distance in distances_puzzle[key]:
                lower_bound = dist_list[int(distance) - 1]
                upper_bound = dist_list[int(distance)]
                interpol[pair_index] = interpol[pair_index] + (distance - int(distance)) * (upper_bound - lower_bound)
            pair_index += 1
    interpol = sum(interpol)
    return interpol


if __name__ == '__main__':
    for PDB_File in MY_PDBS:
        Parsed_PDBs.append(resources.PDB_Parser(PDB_File).exec())

    scoring(scoring_dict, Parsed_PDBs)
