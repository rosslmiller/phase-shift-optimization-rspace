from typing import List

# Local imports
from input_reader import InputReader
from output_reader import StateResults


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

        # Sanjoy wants the first 1P1 coefficient to remain fixed
        if state.upper() == "1P1":
            for i, coeff in zip(self.line_numbers[state], coefficients):
                if i == 0:
                    continue
                self._lines[i] = inject_coefficient_into(self._lines[i], coeff)
            return

        for i, coeff in zip(self.line_numbers[state], coefficients):
            self._lines[i] = inject_coefficient_into(self._lines[i], coeff)

    def determine_target_states_from(self, results: StateResults):
        # The following dependencies were shown to exist by direct computation.
        # See `docs/determine_dependent_coefficients.ipynb`

        target_states = []

        # Check for 1S0 dependencies
        if any(state in results for state in ["1S0", "1D2", "1G4"]):
            target_states.append("1S0")

        # Check for 1P1 dependencies
        if any(state in results for state in ["1P1", "1F3"]):
            target_states.append("1P1")

        # Check for 3S1 dependencies
        if any(
            state in results
            for state in ["3S1", "3D1", "E1", "3D2", "3D3", "3G3", "E3", "3G4"]
        ):
            target_states.append("3S1")

        # Check for 3PJ dependencies
        if any(
            state in results
            for state in ["3P0", "3P1", "3P2", "3F2", "E2", "3F3", "3F4", "3H4", "E4"]
        ):
            target_states.append("3PJ")

        self.target_states = target_states

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
