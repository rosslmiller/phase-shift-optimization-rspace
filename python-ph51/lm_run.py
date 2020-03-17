import argparse
import functools
import scipy.optimize
import subprocess

# local imports
from output_reader import OutputReader
from input_writer import InputWriter
from filenames import (
    PHASE_SHIFT_EXECUTABLE,
    PHASE_SHIFT_INPUT_FILE,
    PHASE_SHIFT_OUTPUT_FILE,
)
from verify_lsj import verify_lsj

from literature_values import _1S0, _3S1


def main():
    args = parse_args()
    input_writer = InputWriter(PHASE_SHIFT_INPUT_FILE)

    if not args.test_case:
        modify_input_file(input_writer, args)

    run_phase_shift_executable()

    output_reader = OutputReader(PHASE_SHIFT_OUTPUT_FILE)
    results = output_reader.get_results_from_output()
    input_writer.determine_modifiable_states_from_output(results)

    states_to_check = list(results)
    print(f"Checking states: {states_to_check!r}")

    if args.init_random:
        insert_random_coefficients(input_writer)

    target_function = functools.partial(
        compute_phase_shifts,
        input_writer=input_writer,
        output_reader=output_reader,
        states_to_check=states_to_check,
        command_line_args=args,
    )

    initial_coefficients = input_writer.initial_coefficients

    least_squares_output = scipy.optimize.leastsq(
        func=target_function,
        x0=initial_coefficients,
        ftol=1e-11,  # max relative error of squares sum
        xtol=1e-13,  # max rel error in approximate solution
        epsfcn=1e-10,  # parameter step size for jacobian approximation.
        # must be >> 1e-sigfig = (rounding treshold) for plugging
        # data into input file. See modify_input() in interface_functions
        maxfev=10000,  # max number of function evaluations
        factor=0.1,
        full_output=True,  # show all auxiliary info in the output.
    )

    format_and_print(initial_coefficients, least_squares_output)
    #TODO(ross): Make the following assignment and print statments a function.
    #TODO: Add Deuteron calculations to script
    if "1S0" in results:
        a_theory=results["1S0"].low_energy_params.a
        delta_a=results["1S0"].low_energy_params.a - _1S0["a_np"]
        r_theory=results["1S0"].low_energy_params.r
        delta_r=results["1S0"].low_energy_params.r - _1S0["r_np"]
        print("For 1S0: \n")
        print("a is calculated to be {a_th} and the difference from experiment is: {da}.\n".format(a_th = a_theory, da = delta_a ) )
        print("r is calculated to be {r_th} and the difference from experiment is: {dr}.\n".format( r_th = r_theory, dr = delta_r) )
    if "3S1" in results:
        a_theory=results["3S1"].low_energy_params.a
        delta_a=results["3S1"].low_energy_params.a - _3S1["a_t"]
        r_theory=results["3S1"].low_energy_params.r
        delta_r=results["3S1"].low_energy_params.r - _3S1["r_t"]
        print("For 3S1: \n")
        print("a is calculated to be {a_th} and the difference from experiment is: {da}.\n".format(a_th = a_theory, da = delta_a ) )
        print("r is calculated to be {r_th} and the difference from experiment is: {dr}.\n".format( r_th = r_theory, dr = delta_r) )
        #print("The deuteron binding energy is calculated to be {B_d} and the difference from experiment is: {dB_d}.\n ".format(B_d= , dB_d=))
        #print("The d-state probability is calcuated to be {P_D} and the difference from experiment is: {dP_D}.\n".format(P_D= , dP_D=))
    return initial_coefficients, least_squares_output


def parse_args():
    parser = argparse.ArgumentParser(description="")  # TODO
    parser.add_argument("--test-case", action="store_true")
    parser.add_argument("--init-random", action="store_true")
    parser.add_argument("--skip-verify-lsj", action="store_true")

    # The following default to None if not specified
    parser.add_argument("--J-start", type=int)
    parser.add_argument("--J-end", type=int)
    parser.add_argument("--singlet", type=int, choices=(0, 1))
    parser.add_argument("--triplet", type=int, choices=(0, 1))
    parser.add_argument("--coupled", type=int, choices=(0, 1))
    parser.add_argument("--ke", dest="kinetic_energies_in_MeV", nargs="+", type=float)
    parser.add_argument("--nuclear-states", nargs="+", type=str.upper)
    parser.add_argument("--coefficients", nargs="+", type=float)
    args = parser.parse_args()
    return args


def modify_input_file(input_writer, args):
    if args.J_start is not None:
        input_writer.J_start = args.J_start

    if args.J_end is not None:
        input_writer.J_end = args.J_end

    if args.singlet is not None:
        input_writer.singlet_states_enabled = bool(args.singlet)

    if args.triplet is not None:
        input_writer.triplet_states_enabled = bool(args.triplet)

    if args.coupled is not None:
        input_writer.coupled_states_enabled = bool(args.coupled)

    # Nuclear states can only be checked after J_start, J_end, and singlet/
    # triplet/coupled in order to determine if the provided nuclear states
    # are valid.
    if args.nuclear_states is not None:
        input_writer.nuclear_states_to_modify = args.nuclear_states

    if args.kinetic_energies_in_MeV is not None:
        input_writer.kinetic_energies_in_MeV = args.kinetic_energies_in_MeV

    if args.coefficients is not None:
        if args.nuclear_states is None:
            raise CoefficientsWithoutNuclearStates(
                "--nuclear-states must be specified if coefficients are provided"
            )

        input_writer.modify_coefficients(args.coefficients)
        # TODO(ben): Should we create a .coefficients property?

    input_writer.write_lines()


def run_phase_shift_executable():
    process = subprocess.run(["./" + PHASE_SHIFT_EXECUTABLE])
    assert process.returncode == 0


def insert_random_coefficients(input_writer):
    rand_coeff = input_writer.create_random_coefficients()
    input_writer.modify_coefficients(rand_coeff)
    input_writer.write_lines()


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


def format_and_print(initial_coefficients, least_squares_output):
    print("Output from scipy.optimize.leastsq:")
    print(f"\n{least_squares_output}\n")
    print(f"Initial parameters: {initial_coefficients!r}")
    final_coefficients = least_squares_output[0]
    print(f"Final parameters:   {final_coefficients.tolist()!r}")


class LevenbergMarquardtException(Exception):
    """Base exception for this module."""


class CoefficientsWithoutNuclearStates(LevenbergMarquardtException):
    """The user supplied coefficients but did not supply nuclear states."""


if __name__ == "__main__":
    main()
