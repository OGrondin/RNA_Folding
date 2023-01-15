'''
This third script takes as input the list of structures from the RNA_Puzzles3 folder and gives them a pseudoenergy score.
It is based on linear interpolation of the scores from the objective function trained by training.py.
These linearly interpolated values are then summed to amount to the Gibbs pseudoenergy score.
'''



if __name__ == '__main__':
