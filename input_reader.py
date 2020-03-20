from collections import Counter
from typing import Tuple


class InputReader:
    @property
    def filename(self):
        return self._filename

    @property
    def lines(self):
        return list(self._lines)

    @property
    def J_range(self) -> Tuple[int, int]:
        i = self.J_range_line_index
        label, lower, upper = self._lines[i].split()
        return (int(lower), int(upper))

    @property
    def singlet_states_enabled(self) -> bool:
        return self.check_multiplicity_selector("singlet")

    @property
    def triplet_states_enabled(self) -> bool:
        return self.check_multiplicity_selector("triplet")

    @property
    def coupled_states_enabled(self) -> bool:
        return self.check_multiplicity_selector("coupled")

    @property
    def states(self):
        return list(self.multiplicity.keys())

    @property
    def coefficients(self):
        return {state: self.coefficients_for(state) for state in self.states}

    def __init__(self, filename):
        self._filename = filename
        self.read_lines()
        self.locate_J_range()
        self.locate_multiplicity_selector()

        # The following order of operations is intentional
        self.compute_state_multiplicities()
        self.locate_coefficients()

    def read_lines(self):
        with open(self.filename, "r") as f:
            self._lines = f.readlines()

    def locate_J_range(self):
        self.J_range_line_index = self.index_of_line_starting_with("jb,je")

    def locate_multiplicity_selector(self):
        self.multiplicity_selector_line_index = self.index_of_line_starting_with("sing,tr,cp")

    def locate_coefficients(self):
        self.line_numbers = {
            state: self.line_numbers_for(state) for state in self.states
        }

    def line_numbers_for(self, state: str):
        # NOTE: we assume state names, such as "1S0", do not appear
        # in /any/ other lines of the input file.
        return tuple(i for i, line in enumerate(self.lines) if state in line.upper())

    def compute_state_multiplicities(self):
        self.multiplicity = Counter(
            get_state_from(line)
            for line in self.lines
            if contains_coefficient_label(line)
        )

    def coefficients_for(self, state: str):
        return tuple(
            get_coefficient_from(self._lines[i]) for i in self.line_numbers[state]
        )

    def index_of_line_starting_with(self, target: str):
        for i, line in enumerate(self._lines):
            if target in line:
                return i
        raise LineNotFound("no line containing {target!r} was found.")

    def check_multiplicity_selector(self, target: str):
        i = self.multiplicity_selector_line_index
        label, singlet, triplet, coupled = self._lines[i].split()

        if target.lower() in ["singlet", "sing", "s"]:
            return int(singlet) != 0
        if target.lower() in ["triplet", "tr", "t"]:
            return int(triplet) != 0
        if target.lower() in ["coupled", "cp", "c"]:
            return int(coupled) != 0

        raise MultiplicityTypeNotFound(target)


def get_state_from(line: str):
    label, state, *rest = line.split()
    return state.upper()


COEFFICIENT_START_INDEX = 10
MAX_COEFFICIENT_WIDTH = 10


def get_coefficient_from(line: str):
    start = COEFFICIENT_START_INDEX
    end = start + MAX_COEFFICIENT_WIDTH
    coefficient = line[start:end]
    return float(coefficient)


COEFFICIENT_LINE_LABELS = [
    "c01",
    "c00",
    "c10",
    "s120",
    "ls0",
    "c11",
    "ls1",
    "s121",
]


def contains_coefficient_label(line: str):
    return any(line.startswith(label) for label in COEFFICIENT_LINE_LABELS)


class InputReaderException(Exception):
    """Base exception for this module."""


class LineNotFound(InputReaderException):
    """Could not locate line within input file."""


class MultiplicityTypeNotFound(InputReaderException):
    """When attempting the read the multiplicity selector, an unknown
    target string was supplied."""
