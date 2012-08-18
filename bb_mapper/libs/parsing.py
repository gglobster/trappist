import re, numpy as np
from loaders import td_txt_file_load
from config import mtype
from array_tetris import extract_nonzero, clump_rows

def mauver_load2_k0(file, threshold):
    """Parse Mauve coordinates file to extract segment coordinates.

    This loads the coordinates data into a Numpy array. All rows that contain
    a zero value are deleted from the array. Rows are then collapsed if the
    coordinates are close (under a user-specified threshold) and in the same
    orientation. The resulting array is then returned.

    Important note: this function only works for pairwise alignments!

    """
    # load file data into numpy array
    raw_array = np.loadtxt(file, skiprows=1, dtype=mtype)
    # set up an array stub to receive non-zero rows (leaves a (1,1) pair)
    stub_array = np.ones(1, dtype=mtype)
    # eliminate rows containing elements of value 0
    try:
        nz_array = extract_nonzero(raw_array, stub_array)
    except TypeError:
        nz_array = np.append(stub_array, raw_array)
    # collapse rows
    cl_array = clump_rows(nz_array, threshold)
    return cl_array

def parse_clustal_idstars(filename):
    """Parse ClustalW output file to estimate identity percentage."""
    #from analysis.text_manipulation import td_txt_file_load
    raw_lines = td_txt_file_load(filename, 3)
    single_line = ''.join(raw_lines)
    idnstar = re.compile(r"\*")
    idnalls = re.findall(idnstar, single_line)
    idntot = len(idnalls)
    return idntot # total number of identical nucleotide positions
