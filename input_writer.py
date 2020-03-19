from typing import List

# Local imports
from input_reader import InputReader


class InputWriter(InputReader):
    @property
    def target_states(self):
        return tuple(self._target_states)

    @target_states.setter
    def target_states(self, states: List[str]):
        # TODO(ben): validate the list of states
        self._target_states = tuple(states)
        self.compute_coefficient_slices()

    @property
    def target_coefficients(self):
        if self._target_states is None:
            raise TargetStatesNotSet(
                "please specify which partial wave coefficients are being "
                "optimized by assigning a list of strings to the `target_states`"
                'attribute (e.g., `x.target_states = ["1S0", "1P1"]`)'
            )
        coefficients = []
        for state in self.target_states:
            coefficients.extend(self.coefficients_for(state))
        return coefficients

    def __init__(self, filename: str):
        super().__init__(filename)
        self._target_states = None

    def compute_coefficient_slices(self):
        start = 0
        self.coefficient_slices = {}
        for state in self.target_states:
            end = start + self.multiplicity[state]
            self.coefficient_slices[state] = slice(start, end)
            start += self.multiplicity[state]

    def modify_coefficients(self, coefficients: List[float]):
        for state in self.target_states:
            s = self.coefficient_slices[state]
            self.set_state_coefficients(state, coefficients[s])

    def set_state_coefficients(self, state: str, coefficients: List[float]):
        assert len(coefficients) == self.multiplicity[state]
        for i, coeff in zip(self.line_numbers[state], coefficients):
            self._lines[i] = inject_coefficient_into(self._lines[i], coeff)

    def write_lines(self):
        with open(self.filename, "w") as f:
            f.writelines(self.lines)


COEFFICIENT_START_COLUMN = 10
MAX_COEFFICIENT_WIDTH = 10


def inject_coefficient_into(line: str, coeff: float):
    coeff = format_coefficient_as_string(coeff)
    character_list = list(line)
    start = COEFFICIENT_START_COLUMN
    end = start + MAX_COEFFICIENT_WIDTH
    character_list[start:end] = list(coeff)
    new_line = "".join(character_list)
    return new_line


def format_coefficient_as_string(coeff: float):
    coeff = f"{coeff: < 10.6f}"
    coeff = truncate_coefficient_to_max_width(coeff)
    return coeff


def truncate_coefficient_to_max_width(coeff: str):
    return coeff[:MAX_COEFFICIENT_WIDTH]


class InputWriterException(Exception):
    """Base exception for this module."""


class TargetStatesNotSet(InputWriterException):
    """The InputWriter.target_states attribute has not yet been set;
    the instance does not know which nuclear state's coefficients are
    being optimized."""
