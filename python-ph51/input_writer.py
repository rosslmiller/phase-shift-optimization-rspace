import random
import time

from typing import List, Tuple, Dict

# Local imports
from output_reader import StateResults


class NuclearState:
    orbitals = "SPDFGHIJKLMNOPQRSTUVWXYZ"


class InputWriter:

    _set_of_nuclear_states = "1S0 1P1 3S1 3PJ".split()
    # TODO(ben): should we determine these states by reading the input file
    # instead of hard-coding them here?

    def __init__(self, filename: str):
        pass

    def modify_coefficients(self, coefficients: List[float]):
        pass

    def _translate_to_input_states(self, output_states: List[str]) -> List[str]:
        input_states = []
        for state in output_states:
            if state in self._set_of_nuclear_states:
                input_states.append(state)
            elif state == "E1":
                input_states.append("3SD1")
            elif state == "E2":
                input_states.append("3PF2")
            elif state == "E3":
                input_states.append("3D-G3")
            elif state == "E4":
                input_states.append("3F-H4")
            elif state == "E5":
                input_states.append("3G-I5")
            else:
                raise NotImplementedError(f"output state: {state}")
        return input_states

    def write_lines(self):
        assert self.J_start <= self.J_end
        # TODO(ben): consider J_range property
        with open(self._output_filename, "w") as f:
            f.writelines(self.lines)


def truncate_coefficient_to_max_width(coefficient: float) -> str:
    coefficient = f"{coefficient: < .15f}"
    return coefficient[:MAX_COEFFICIENT_WIDTH]


class InputWriterException(Exception):
    """Base exception for this module."""
