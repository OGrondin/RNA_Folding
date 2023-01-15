import numpy as np


class PDB_Parser:
    '''
    A class common to both the training.py and scoring.py, to process PDB files and get desired output format (C3 atoms)
    '''
    def __init__(self, PDB_file):
        self.path = PDB_file
        self.output = list()

    def split_pdb_line(self, line):
        '''
        Splitting each line by fixed length to retrieve PDB Data in a list
        :param line: One line of a PDB file
        :return: A list representing all the fields of line
        '''
        #PDB_HEADER = ("Atom_type", "Atom_#", "Atom_name", "Alt_loc", \
        # "Residue", "Chain_ID", "Residue_#", "Insertion", \
        # "X_coords", "Y_coords", "Z_coords", "Occupancy", \
        # "Temp_factor", "Element", "Charge")

        line = [line[0:6], line[6:11], line[12:16], line[16:17], \
                line[17:20], line[21:22], line[22:26], line[26:27], \
                line[30:38], line[38:46], line[46:54], line[54:60], \
                line[60:66], line[76:78], line[78:80]]
        line = [field.replace(' ', '') for field in line]
        return line

    def pdb_ATOM_df(self):
        '''
        Main PDB parsing function.
        :return: List of lines containing C3 atoms
        '''
        with open(self.path, 'r') as f:
            self.output = [k.strip("\n") for k in f.readlines()]
            self.output = [l for l in self.output if (l[0:4] == "ATOM")]
        # Parsing as a list of lists, one list per line
        self.output = [self.split_pdb_line(l) for l in self.output]
        # Selecting only rows corresponding to C3 atoms of the first chain "A"
        self.output = [l for l in self.output if l[2] == "C3'"]


    def exec(self):
        print("Log: Parsing data from:", self.path)
        self.pdb_ATOM_df()
        print("Number of nucleotides", len(self.output))
        return self.output

## Math functions

def dist_3D(coords_series):
    '''
    Given a list of six coordinates, outputs their 3D euclidean distance.
    :param coords_series: List following "x1, y1, z1, x2, y2, z2" prototype for two points 1 and 2.
    :return: Euclidean distance of coordinates
    '''
    x1, y1, z1, x2, y2, z2 = coords_series
    return ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ** .5

## Print functions

# From sth on https://stackoverflow.com/questions/3229419/how-to-pretty-print-nested-dictionaries
def pretty(d, indent=0):
   for key, value in d.items():
      print('\t' * indent + str(key))
      if isinstance(value, dict):
         pretty(value, indent+1)
      else:
         print('\t' * (indent+1) + str([round(i,5) for i in value]))

## File management

def PDB_list_path(list_PDB_IDs):
    for i, x in enumerate(list_PDB_IDs):
        list_PDB_IDs[i] = "PDB_Training/" + str(x) + ".pdb"
    return list_PDB_IDs

def pseudo_score(observed, reference):
    '''
    Computes the minimum between -log(observed/ref) and 10.
    :param observed: Observed frequency for known (pair, distance bin)
    :param reference: Reference frequency for known: a) distance bin (training.py) b) (pair, distance bin) (scoring.py)
    :return: 
    '''
    if reference != 0:
        return (min((- np.log(observed/reference)), 10))
    else:
        return(10)
