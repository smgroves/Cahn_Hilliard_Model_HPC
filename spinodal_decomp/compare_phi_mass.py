import numpy as np
import pandas as pd
import itertools
import os
import re
import sys
import traceback

def read_specific_lines(file_path, line_numbers):
    result = []
    with open(file_path, "r") as file:
        for current_line_number, line in enumerate(file):
            if current_line_number in line_numbers:
                mass = line.strip("\n")
                result.append(mass)
            if current_line_number > max(line_numbers):
                break
    result = np.array(result)
    return result

# NMG_MATLAB_2000_dt_5.50e-06_Nx_512_neumann_n_relax_4_75p_phi

def compare_phi_masses(indir,out_file, timepoints):


    # list_of_n = []
    # list_of_filenames = []

    # pattern = re.compile(fr'Nx_(\d+)_([a-zA-Z]+).*phi.csv$')
    pattern = re.compile(fr'^([A-Za-z0-9]+)_([A-Za-z0-9]+)_(\d+).*?_Nx_(\d+).*?_(neumann|periodic).*?(?:_dtout_(\d+))?.*?_?(\d+p)?_?mass\.csv$')

    for root, dirs, files in os.walk(indir):
        # grab only files that match the structure
        for fname in files:
            try:
                clean_fname = re.sub(r'^.*?([A-Z]{2,4})', r'\1', fname)

                match = re.match(pattern, clean_fname)

                if match:
                    method = match.group(1)
                    language = match.group(2)
                    total_timepoints = match.group(3)
                    nx_val = match.group(4)
                    bc = match.group(5)
                    dt_out = match.group(6) or "None (10)"
                    ic = match.group(7) or "50p"

                    # print(f"{method}, {language},{total_timepoints},{nx_val},{bc},{dt_out},{ic}")

                    if dt_out == "None (10)":
                        dt_out_practical = 10
                    else:
                        dt_out_practical = int(dt_out)

                        lines = [int(i/int(dt_out_practical)) for i in timepoints]
                        masses = read_specific_lines(f"{root}/{fname}",lines)
                    for t, timepoint in enumerate(timepoints):
                        T = {}
                        T['language'] = language
                        T['method'] = method
                        T['GridSize'] = nx_val
                        T['boundary'] = bc
                        T['ic'] = ic
                        T['dt_out'] = dt_out
                        T["timepoint"] = t
                        T['mass'] = masses[t]
                        T['pathname'] = f"{root}/{fname}"

                        T = pd.DataFrame([T])
                        if not os.path.isfile(out_file):
                            T.to_csv(out_file, index=False)
                        else:
                            with open(out_file, "a") as f:
                                T.to_csv(f, header=False, index=False)
            except BaseException as ex:
                # Get current system exception
                ex_type, ex_value, ex_traceback = sys.exc_info()

                # Extract unformatter stack traces as tuples
                trace_back = traceback.extract_tb(ex_traceback)

                # Format stacktrace
                stack_trace = list()

                for trace in trace_back:
                    stack_trace.append("File : %s , Line : %d, Func.Name : %s, Message : %s" % (trace[0], trace[1], trace[2], trace[3]))
                print(fname)
                print("Exception type : %s " % ex_type.__name__)
                print("Exception message : %s" %ex_value)
                print("Stack trace : %s" %stack_trace)
                # return list_of_filenames, list_of_n

#reset file each time if rerunning the same timepoints
out_file = f"/home/xpz5km/Cahn_Hilliard_Model_HPC/spinodal_decomp/compare_masses.csv"


indir = "/project/g_bme-janeslab/SarahG/spinodal_decomp_04_2025/out_julia"

compare_phi_masses(indir, out_file, timepoints = [0,1000,2000])



