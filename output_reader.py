import re

from collections import namedtuple

# TODO(ben): What are the technical names for
#     1) terms like "1s0", "3p0", "1p1", etc.; and
#     2) terms like "e1", "e2", "e3"
# Answer: I've settled on "states" (as in nuclear energy state)
type_1_state_np = re.compile(r"^\s*(\d [spdfghijk] \d)\s+$")
type_1_state_pp = re.compile(r"^\s*(\d [spdfghijk] \d)\s+p-p\s*$")
type_2_state_np = re.compile(r"^\s*(e \d)\s+$")
type_2_state_pp = re.compile(r"^\s*(e \d)\s+p-p\s*$")

# TODO(ben): Will we ever see phases such as '3 S D 1' or '3 F - H 4' in the output?
# Answer: No, these only appear in the input.

phase_shift_result = re.compile(
    r"""
    ^                                           # beginning of line
    \s+(?P<elab>\d+\.\d+)                       # kinetic energy, I think...
    \s+(?P<theoretical>-?\d\.\d+[eE][\+-]\d+)   # theoretical value
    \s+(?P<experimental>-?\d\.\d+[eE][\+-]\d+)  # experimental value
    \s+(?P<upper>-?\d\.\d+[eE][\+-]\d+)         # experimental upper bound
    \s+(?P<lower>-?\d\.\d+[eE][\+-]\d+)         # experimental lower bound
    \s(?:n93n)?                                 # optional string; unsure of meaning
    \s+(?P<chi2>\d+\.\d+)                       # chi squared
    $                                           # end of line
""",
    re.VERBOSE,
)

low_energy_params = re.compile(
    r"""
    ^                                          # beginning of line
    \s+(?:low\ energy\ parameters)             # literal string: low energy parameters
    \s+a\s+=\s+(?P<a>-?\d\.\d+[eE][\+-]\d+)    # a = ...
    \s+r\s+=\s+(?P<r>-?\d\.\d+[eE][\+-]\d+)    # r = ...
    $                                          # end of line
""",
    re.VERBOSE,
)


PhaseShift = namedtuple("PhaseShift", "elab theoretical experimental upper lower chi2")
LowEnergyParams = namedtuple("LowEnergyParams", "a r")


class StateResults:
    def __init__(self, phase_shifts, low_energy_params=None):
        self._phase_shifts = phase_shifts
        self._low_energy_params = low_energy_params

    @property
    def phase_shifts(self):
        return self._phase_shifts

    @property
    def low_energy_params(self):
        return self._low_energy_params


class OutputReader:

    # TODO: Rewrite this entire class, applying principles from Clean Code
    order_of_nuclear_states = [
        "1S0",
        "3P0",
        "1P1",
        "3P1",
        "3S1",
        "3D1",
        "E1",
        "1D2",
        "3D2",
        "3P2",
        "3F2",
        "E2",
        "1F3",
        "3F3",
        "3D3",
        "3G3",
        "E3",
        "1G4",
        "3G4",
        "3F4",
        "3H4",
        "E4",
        "1H5",
        "3H5",
        "3G5",
        "3I5",
        "E5",
        "1I6",
        "3I6",
        "3H6",
        "3K6",
        "E6",
        "1K7",
        "3K7",
        "3I7",
        "E7",
    ]

    def __init__(self, filename: str):
        self._filename = filename
        self.count_elabs()
        self.create_ordinals_dict()

    @property
    def filename(self):
        return self._filename

    def get_results_from_output(self):
        with open(self.filename, "r") as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if "p h a s e - s h i f t s" in line:
                break

        line_iter = iter(lines[i:])
        results = dict()

        def is_blank(line):
            return line.strip() == ""

        def parse_phase_shift(line):
            elab, th, exp, up, low, *maybe, chi2 = line.split()
            return PhaseShift(
                float(elab), float(th), float(exp), float(up), float(low), float(chi2)
            )

            match = re.match(phase_shift_result, line)
            assert match
            return PhaseShift(**{k: float(v) for k, v in match.groupdict().items()})

        for line in line_iter:

            match = (
                re.match(type_1_state_np, line)
                or re.match(type_2_state_np, line)
                or re.match(type_1_state_pp, line)
                or re.match(type_2_state_pp, line)
            )

            if match:
                state = match.group(1)
                state = (
                    state.strip().replace(" ", "").upper()
                )  # Remove spaces + make uppercase
                assert is_blank(next(line_iter))

                phase_shifts = tuple(
                    parse_phase_shift(next(line_iter)) for i in range(self.num_elabs)
                )

                # Check for low energy parameters
                next_line = next(line_iter)
                if "total chi**2" in next_line:
                    # No low energy parameters displayed
                    le_params = None
                else:
                    assert is_blank(next_line)
                    next_line = next(line_iter)
                    lep_match = re.match(low_energy_params, next_line)

                    le_params = (
                        LowEnergyParams(
                            **{k: float(v) for k, v in lep_match.groupdict().items()}
                        )
                        if lep_match
                        else None
                    )

                # Store results in dictionary
                results[state] = StateResults(phase_shifts, le_params)

        return results

    def elab_generator(self):
        with open(self.filename, "r") as f:
            # Skip over first few lines of input parameters
            for line in f:
                if "wn" in line:
                    break

            # Expect a blank line just before elab values
            assert f.readline().strip() == ""

            # Return floats until we hit "elab end"
            for line in f:
                if "elab end" in line:
                    return
                yield float(line.strip())

    def count_elabs(self):
        self.num_elabs = len(list(self.elab_generator()))

    def create_ordinals_dict(self):
        self._ordinal_for_nuclear_state = {
            s: i for i, s in enumerate(self.order_of_nuclear_states)
        }

    def ordinal_for_nuclear_state(self, state: str):
        return self._ordinal_for_nuclear_state[state]
