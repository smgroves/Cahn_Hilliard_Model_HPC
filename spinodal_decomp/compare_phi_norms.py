import numpy as np
import pandas as pd
import itertools
import os
import re

def read_specific_lines(file_path, line_numbers):
    result = []
    with open(file_path, "r") as file:
        for current_line_number, line in enumerate(file):
            if current_line_number in line_numbers:
                digits = [float(i) for i in line.strip("\n").split(" ")]
                result.append(digits)
            if current_line_number > max(line_numbers):
                break
    result = np.array(result)
    return result

# NMG_MATLAB_2000_dt_5.50e-06_Nx_512_neumann_n_relax_4_75p_phi

def compare_phis_norm_at_timepoint(indir, GridSize, timepoint=None):
    first_line = timepoint * GridSize
    last_line = first_line + GridSize
    line_list = range(first_line, last_line)

    list_of_n = []
    list_of_filenames = []
    pattern = re.compile(fr'Nx_(\d+)_([a-zA-Z]+).*phi.csv$')
    for root, dirs, files in os.walk(indir):
        # grab only files that match the structure
        for fname in files:
            match = pattern.search(fname)
            if match:
                nx_val = int(match.group(1))
            if nx_val == GridSize:
                print(fname)
                snapshot = read_specific_lines(fname, line_list)
                n = np.linalg.norm(snapshot)
                list_of_filenames.append(fname)
                list_of_n.append(n)
    return list_of_filenames, list_of_n

indir = ""
GridSize = 64
print(compare_phis_norm_at_timepoint(indir, GridSize))


