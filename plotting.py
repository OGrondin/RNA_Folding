'''
This second script takes as input the Gibbs energy pseudoscore for each couple (pair, bin), where pair
is a pairing between two nucleotides (A, C, T or U) and bin is a distance category of 1 Angstrom range between 0 and 20.
'''

import json
import matplotlib.pyplot as plt

import resources

to_plot = dict()
with open("Results/Score_Output.json", 'r') as json_file:
    to_plot = json.load(json_file)

def plot(show_graph = False):
    for key in to_plot:
        data = to_plot[key]
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel('Distance in Angstroms')
        ax.set_ylabel('Pseudoenergy score')
        plt.plot(range(20), data)
        ax.set_title('Pseudoenergy per distance to each other: Pair {}'.format(key))
        if show_graph:
            plt.show()
        else:
            plt.savefig("Results/Plots/" + key + ".png")

if __name__ == '__main__':
    plot()