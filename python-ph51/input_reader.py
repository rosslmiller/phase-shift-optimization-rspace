from collections import Counter


class InputReader:
    @property
    def filename(self):
        return self._filename

    @property
    def lines(self):
        return list(self._lines)

    @property
    def states(self):
        return list(self.multiplicity.keys())

    @property
    def coefficients(self):
        return {state: self.coefficients_for(state) for state in self.states}

    def __init__(self, filename):
        self._filename = filename
        self.read_lines()
        self.compute_state_multiplicities()
        self.locate_coefficients()

    def read_lines(self):
        with open(self.filename, "r") as f:
            self._lines = f.readlines()

    def compute_state_multiplicities(self):
        self.multiplicity = Counter(
            get_state_from(line)
            for line in self.lines
            if contains_coefficient_label(line)
        )

    def locate_coefficients(self):
        self.line_numbers = {
            state: self.line_numbers_for(state) for state in self.states
        }

    def line_numbers_for(self, state: str):
        # NOTE: we assume state names, such as "1S0", do not appear
        # in /any/ other lines of the input file.
        return tuple(i for i, line in enumerate(self.lines) if state in line.upper())

    def coefficients_for(self, state: str):
        return tuple(
            get_coefficient_from(self._lines[i]) for i in self.line_numbers[state]
        )


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
