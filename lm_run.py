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

from literature_values import _1S0, _3S1


TARGET_STATES = ["1S0"]


def main():
    input_writer = InputWriter(PHASE_SHIFT_INPUT_FILE)

    run_phase_shift_executable()

    output_reader = OutputReader(PHASE_SHIFT_OUTPUT_FILE)
    results = output_reader.get_results_from_output()
    input_writer.target_states = TARGET_STATES

    states_to_check = list(results)
    print(f"Checking states: {states_to_check!r}")

    target_function = functools.partial(
        compute_phase_shifts,
        input_writer=input_writer,
        output_reader=output_reader,
        states_to_check=states_to_check,
    )

    initial_coefficients = input_writer.target_coefficients

    least_squares_output = scipy.optimize.leastsq(
        func=target_function,
        x0=initial_coefficients,
        ftol=1e-11,  # max relative error of squares sum
        xtol=1e-13,  # max rel error in approximate solution
        epsfcn=1e-10,  # parameter step size for Jacobian approximation.
        # must be >> 1e-sigfig = (rounding threshold) for plugging
        # data into input file.
        maxfev=10000,  # max number of function evaluations
        factor=0.1,
        full_output=True,  # show all auxiliary info in the output.
    )

    format_and_print(initial_coefficients, least_squares_output)

    # Fetch results after optimization
    results = output_reader.get_results_from_output()
    show_error_in_a_and_r(results)

    return initial_coefficients, least_squares_output


def run_phase_shift_executable():
    process = subprocess.run(["./" + PHASE_SHIFT_EXECUTABLE])
    assert process.returncode == 0


def compute_phase_shifts(coefficients, input_writer, output_reader, states_to_check):
    input_writer.modify_coefficients(coefficients)
    input_writer.write_lines()
    run_phase_shift_executable()
    results = output_reader.get_results_from_output()

    # Compute chi values
    chi = [
        (ps.theoretical - ps.experimental) / (ps.upper - ps.experimental)
        for state in states_to_check
        for ps in results[state].phase_shifts
        if ps.upper != ps.experimental
    ]

    # TODO: Add if's for 1S0 case: add computation for a and r values, np vs pp cases
    # TODO: Add ifs for 3S1, 3D1 case: add computation for a_t, r_t, B_D, P_d values

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


def show_error_in_a_and_r(results):
    # TODO: Add Deuteron calculations to script
    if "1S0" in results:
        at = results["1S0"].low_energy_params.a
        da = results["1S0"].low_energy_params.a - _1S0["a_np"]
        rt = results["1S0"].low_energy_params.r
        dr = results["1S0"].low_energy_params.r - _1S0["r_np"]
        print("\nFor 1S0:")
        print(f"a is calculated to be {at}; difference from experiment is: {da}")
        print(f"r is calculated to be {rt}; difference from experiment is: {dr}")
    if "3S1" in results:
        at = results["3S1"].low_energy_params.a
        da = results["3S1"].low_energy_params.a - _3S1["a_t"]
        rt = results["3S1"].low_energy_params.r
        dr = results["3S1"].low_energy_params.r - _3S1["r_t"]
        print("\nFor 3S1:")
        print(f"a is calculated to be {at}; difference from experiment is: {da}")
        print(f"r is calculated to be {rt}; difference from experiment is: {dr}")
        # print("The deuteron binding energy is calculated to be {B_d} and the difference from experiment is: {dB_d}.\n ".format(B_d= , dB_d=))
        # print("The d-state probability is calcuated to be {P_D} and the difference from experiment is: {dP_D}.\n".format(P_D= , dP_D=))


class LevenbergMarquardtException(Exception):
    """Base exception for this module."""


class CoefficientsWithoutNuclearStates(LevenbergMarquardtException):
    """The user supplied coefficients but did not supply nuclear states."""


if __name__ == "__main__":
    main()
