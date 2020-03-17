import scipy.optimize
import argparse
import functools
import subprocess
from collections import Counter
from typing import List, Tuple, Set, Dict
import random
import time
import re

#Filenames

DEUTERON_EXECUTABLE = "xdqq"
PHASE_SHIFT_EXECUTABLE = "xphqq"

DEUTERON_INPUT_FILE = "ddqq.d"
PHASE_SHIFT_INPUT_FILE = "dphqq.d"

DEUTERTON_OUTPUT_FILE = "deuqq.d"
PHASE_SHIFT_OUTPUT_FILE = "phqq.d"


def main():
    ##Need to write functions for pieces in target_function
    #target_function = functools.partial(
    #    compute_phase_shifts,
    #    input_writier=input_writer,
    #    output_reader=output_reader,
    #    states_to_check=states_to_check,
    #    command_line_args=args,
    #)

    #TODO: Find coefficients in PHASE_SHIFT_INPUT_FILE= "dphqq.d", store them in an array, and print them

    with open(PHASE_SHIFT_INPUT_FILE, "r") as file:
        all_content = file.readlines()

    
    initial_coefficients = []

    _1s0string = "c01\s*1s0\s*(-*\d+\.\d*\s*)[0,1]."
    regex1s0 = re.compile(_1s0string)

    # Written to get coefficient floats
    for i in all_content:
        if regex1s0.findall(i) != []:
            x = regex1s0.findall(i)[0]
            initial_coefficients.append(float( [j for j in filter(None,x.split(' '))][0]))
        
    print(initial_coefficients)





"""
    least_squares_output = scipy.optimize.leastsq(
        func = target_function
        x0 = initial_coefficients,
        ftol=1e-11,  # max relative error of squares sum
        xtol=1e-13,  # max rel error in approximate solution
        epsfcn=1e-10,  # parameter step size for jacobian approximation.
        # must be >> 1e-sigfig = (rounding treshold) for plugging
        # data into input file. See modify_input() in interface_functions
        maxfev=10000,  # max number of function evaluations
        factor=0.1,
        full_output=True,  # show all auxiliary info in the output.
    )



def compute_phase_shifts(
    coefficients, input_writer, output_reader, states_to_check, command_line_args
):
    input_writer.modify_coefficients(coefficients)
    input_writer.write_lines()
    run_phase_shift_executable()
    results = output_reader.get_results_from_output()

    # Compute chi values
    chi = [
        (ps.theoretical - ps.experimental) / (ps.upper - ps.experimental)
        for state in states_to_check
        for ps in results[state].phase_shifts
        if ps.upper != ps.lower
    ]

    # TODO: Add if's for 1S0 case: add computation for a and r values, np vs pp cases



    # TODO: Add ifs for 3S1, 3D1 case: add computation for a_t, r_t, B_D, P_d values


    args = command_line_args
    if not args.skip_verify_lsj:
        verify_lsj()

    chi2 = sum(x ** 2 for x in chi)
    print(f"Sum chi**2 = {chi2:0.4f}")
    print(f"    {coefficients}")
    return chi
"""

if __name__ == "__main__":
    main()
