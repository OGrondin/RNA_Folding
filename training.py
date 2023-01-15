'''
This first script computes distance (0-20 Angstroms, bins of 1 Angstrom range) distributions for the 10 usual RNA pairs.
It takes into account distances between the Carbon 3 atom of every two nucleotides i and j where j > i+3 inside a chain.
It also computes the distribution of 'XX' where X is any nucleotide and the log-ratio of each nucleotide to 'XX'.
'''

import sys
import resources
import time
import numpy as np
import pandas as pd


# Variables
PDB_IDs_training = ["2fqn", "4gxy", "5kpy", "5l4o", "5lyu", "5t83", "5u3g", "6wtl", "6ymc", "7eem"]
original_stdout = sys.stdout

# Input preparation
Parsed_PDBs = list()
MY_PDBS = resources.PDB_list_path(PDB_IDs_training)

def distances_computation(list_PDB):
    '''
    Computes 3D euclidean distances between C3 atoms of nucleotides
    :param Parsed_PDB: A list of C3 atoms, output of the PDB_Parser class
    :return: Dictionnary of distances between i and i+4, i+5 ...
    '''
    distances = dict()
    for Parsed_PDB in list_PDB:
        for i in range(len(Parsed_PDB)):
            for j in range(i+3, len(Parsed_PDB)):
                if Parsed_PDB[i][5] == Parsed_PDB[j][5]: # Testing to ensure intra-chain comparison
                    pair = sorted((Parsed_PDB[i][4], Parsed_PDB[j][4])) # Sorted ensures there are 10 pairings, not 16
                    pair = str(pair[0]) + str(pair[1])
                    vector_dist = Parsed_PDB[i][8:11] + Parsed_PDB[j][8:11]
                    vector_dist = [float(elem) for elem in vector_dist]
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

    # Processing of distances to counts in bins
    for key in distances:
        dist_list = distances[key]
        list_dist = [0 for _ in range(20)]
        for elem in dist_list:
            index = int(elem)
            old = list_dist[index]
            list_dist[index] = old + 1
        distances[key] = list_dist

        # with open("Training_distances.tsv",'w+') as g: !!!!!!!!!!!

    return distances

def frequencies_score(distance_dict):
    '''
    Takes the distance dictionnary as input, computes the frequencies for each pairing, the overall distances frequency,
    as well as the log ratios of observed pairing over
    :param distance_dict: Dictionnary with training distances as counts in bins
    :return: Two dictionnaries with frequency and score values as well as the frequencies of any pairing (list_XX)
    '''
    list_XX = [0 for _ in range(20)]
    frequency_dict = dict()
    # Computing frequency of each pair as well as overall distance frequency
    for key in distance_dict:
        dist_list = distance_dict[key]
        print(list_XX, dist_list)
        list_XX = [i+j for i, j in zip(list_XX, dist_list)]
        freq_list = [elem / sum(dist_list) for elem in dist_list]
        frequency_dict[key] = freq_list
    list_XX = [elem / sum(list_XX) for elem in dist_list]

    score_dict = dict()
    for key in frequency_dict:
        freq_list = frequency_dict[key]
        score_list = [(float(x),float(y)) for x, y in zip(freq_list, list_XX)]
        score_list = [resources.pseudo_score(x,y) for x, y in score_list]
        score_dict[key] = score_list
    return(frequency_dict, list_XX, score_dict)





if __name__ == '__main__':

    start = time.time()

    # Counting distances in bins
    for PDB_File in MY_PDBS:
        Parsed_PDBs.append(resources.PDB_Parser(PDB_File).exec())
    pairwise_distances = distances_computation(Parsed_PDBs)


    # Computing frequencies for (i,j) and (X,X) pairs, where i, j are A, C, G or U and X is any of the above
    freq_pairwise, freq_overall, score_pairwise = frequencies_score(pairwise_distances)

    # Pretty_Output
    with open("Results/Pretty_Output.txt", 'w+') as p:
        sys.stdout = p  # Change the standard output to the file we created.

        print("Pairwise distances matrix", "\n")
        resources.pretty(pairwise_distances)
        print("\n", "Pairwise Frequency matrix: Fobs_ij = N_ij(r) / N_ij, where r is the distance bin", "\n")
        resources.pretty(freq_pairwise)
        print("\n", "Overall Frequency matrix: Fref_XX = N_XX(r) / N_XX, where X is any base (A, C, G or U)", "\n")
        print([round(i, 5) for i in freq_overall])

        print("\n", "Pairwise Score matrix: u_ij = - log (Fobs_ij(r) / Fref_ij(r)", "\n")
        resources.pretty(score_pairwise)

        print("\n","Time spent:", round(time.time() - start, 4), "seconds")

        sys.stdout = original_stdout  # Reset the standard output to its original value