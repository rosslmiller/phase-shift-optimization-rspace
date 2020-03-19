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


low_energy_params = re.compile(
    r"""
    ^                               # beginning of line
    \s+(?:low\ energy\ parameters)  # literal string: low energy parameters
    \s+a\s+=\s+(?P<a>-?\d+\.\d+)    # a = ...
    \s+r\s+=\s+(?P<r>-?\d+\.\d+)    # r = ...
    $                               # end of line
""",
    re.VERBOSE,
)


PhaseShift = namedtuple("PhaseShift", "elab theoretical experimental upper lower chi2")
LowEnergyParams = namedtuple("LowEnergyParams", "a r")


class StateResults:
    @property
    def state(self):
        return self._state

    @property
    def phase_shifts(self):
        return self._phase_shifts

    @property
    def low_energy_params(self):
        return self._low_energy_params

    def __init__(self, state, phase_shifts, low_energy_params=None):
        self._state = state
        self._phase_shifts = phase_shifts
        self._low_energy_params = low_energy_params

    def __sub__(self, other):
        if self.state != other.state:
            raise ValueError(
                "cannot subtract results of two different states "
                f"(self.state == {self.state!r}, other.state == {other.state!r})"
            )

        self_elabs, self_theoretical = zip(
            *[(ps.elab, ps.theoretical) for ps in self.phase_shifts]
        )

        other_elabs, other_theoretical = zip(
            *[(ps.elab, ps.theoretical) for ps in other.phase_shifts]
        )

        if self_elabs != other_elabs:
            raise ValueError("results do not have the same energies")

        return sum((x - y) ** 2 for x, y in zip(self_theoretical, other_theoretical))

    def __repr__(self):
        class_name = type(self).__name__

        string_representation = (
            f"{class_name}(\n" f"    state={self.state!r},\n" "    phase_shifts=(\n"
        )

        for ps in self.phase_shifts:
            string_representation += f"        {ps},\n"
        string_representation += "    ),\n"
        string_representation += (
            f"    low_energy_params={self.low_energy_params}\n" ")\n"
        )

        return string_representation

    def __str__(self):
        header = (
            f'Phase Shifts for "{self.state}":\n'
            "---------------------------------------------------\n"
            "      elab    theoretical    experimental    chi**2\n"
            "---------------------------------------------------\n"
        )

        body = "".join(
            [
                f"{ps.elab: > 10.2f}{ps.theoretical: > 15.5f}"
                f"{ps.experimental: > 16.3f}{ps.chi2: > 10.2f}\n"
                for ps in self.phase_shifts
            ]
        )

        footer = "---------------------------------------------------\n"

        if self.low_energy_params:
            a = self.low_energy_params.a
            r = self.low_energy_params.r
            footer += (
                "\nLow Energy Parameters:\n"
                f"    a = {a:10.4f}\n"
                f"    r = {r:10.4f}\n"
            )
        return header + body + footer


class OutputReader:

    # TODO(ben): Rewrite this entire class, applying principles from Clean Code

    @property
    def filename(self):
        return self._filename

    def __init__(self, filename: str):
        self._filename = filename
        self.count_elabs()

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
                results[state] = StateResults(state, phase_shifts, le_params)

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
