from collections import Counter
from typing import List, Tuple, Set


class InputReader:
    @property
    def filename(self) -> str:
        return self._filename

    @property
    def J_start(self) -> int:
        return self._J_start

    @property
    def J_end(self) -> int:
        return self._J_end

    @property
    def singlet_states_enabled(self) -> bool:
        return self._singlet_states_enabled

    @property
    def triplet_states_enabled(self) -> bool:
        return self._triplet_states_enabled

    @property
    def coupled_states_enabled(self) -> bool:
        return self._coupled_states_enabled

    @property
    def kinetic_energies_in_MeV(self) -> List[float]:
        return list(self._kinetic_energies_in_MeV)

    @property
    def set_of_nuclear_states(self) -> Set[str]:
        return set(self._set_of_nuclear_states)

    @property
    def coefficients(self) -> List[float]:
        return list(self._coefficients)

    def __init__(self, filename: str):
        self._filename = filename
        self.read_current_parameters()

    def read_current_parameters(self) -> None:
        self.read_lines()
        self._locate_all_parameters()
        self._read_parameters()
        # Coefficients must be handled separately because we must
        # determine which nuclear states are available before
        # we can compute coefficient slices / read coefficients.
        self._locate_coefficients()
        self._read_coefficients()

    def read_lines(self) -> None:
        with open(self._filename, "r") as f:
            self.lines = f.readlines()

    def _locate_all_parameters(self) -> None:
        self._locate_J_range_index()
        self._locate_multiplicity_selector_index()
        self._locate_kinetic_energies_slice()
        self._locate_lsj_slice()

    def _read_parameters(self) -> None:
        self._read_J_range()
        self._read_multiplicity_selector()
        self._read_kinetic_energies()
        self._read_nuclear_states()

    def _locate_coefficients(self) -> None:
        slices = {
            state: self._locate_slice_for_state(state)
            for state in self._set_of_nuclear_states
        }
        self.coefficient_slices = slices

    def _locate_J_range_index(self) -> None:
        self.J_range_index = self._get_index_of_line_starting_with("jb,je")

    def _locate_multiplicity_selector_index(self) -> None:
        self.multiplicity_selector_index = self._get_index_of_line_starting_with(
            "sing,tr,cp"
        )

    def _locate_kinetic_energies_slice(self) -> None:
        start = (
            self._get_index_of_line_starting_with("wn") + 1
        )  # "wn" is one line ahead
        end = self._get_index_of_line_starting_with("elab end")
        self.kinetic_energies_slice = slice(start, end)

    def _locate_lsj_slice(self) -> None:
        start = self._get_index_of_line_starting_with("lsj")
        end = self._get_index_of_line_starting_with("end param.")
        self.lsj_slice = slice(start, end)

    def _read_J_range(self) -> None:
        line = self.lines[self.J_range_index]
        label, J_start, J_end = line.split()
        self._J_start = int(J_start)
        self._J_end = int(J_end)

    def _read_multiplicity_selector(self) -> None:
        line = self.lines[self.multiplicity_selector_index]
        label, sing, tr, cp = line.split()
        self._singlet_states_enabled = bool(int(sing))
        self._triplet_states_enabled = bool(int(tr))
        self._coupled_states_enabled = bool(int(cp))

    def _read_kinetic_energies(self) -> None:
        lines = self.lines[self.kinetic_energies_slice]
        kinetic_energies = tuple(float(line.strip()) for line in lines)
        self._kinetic_energies_in_MeV = kinetic_energies

    def _read_nuclear_states(self) -> None:
        lsj_section = self.lines[self.lsj_slice]
        states = []
        for line in lsj_section:
            label, state, *rest = line.split()
            states.append(state)
        self._set_of_nuclear_states = set(states)
        self._nuclear_state_multiplicity = Counter(states)
        self._create_nuclear_state_ordinal_dict(states)

    def _read_coefficients(self) -> None:
        coefficients = []
        for state in self._set_of_nuclear_states:
            slicer = self.coefficient_slices[state]
            lines = self.lines[slicer]
            for line in lines:
                try:
                    _, _, coeff, _, _ = line.split()
                except ValueError as error:
                    msg = error.args[0]
                    assert msg.startswith("not enough values to unpack")
                else:
                    coefficients.append(float(coeff))
        self._coefficients = tuple(coefficients)

    def _create_nuclear_state_ordinal_dict(self, states: List[str]) -> None:
        ordinal_dict = {}
        seen_states = set()
        i = 0
        for state in states:
            if state in seen_states:
                continue
            seen_states.add(state)
            ordinal_dict[state] = i
            i += 1
        self._ordinal_for_nuclear_state = ordinal_dict

    def _get_index_of_line_starting_with(self, string: str) -> int:
        for i, line in enumerate(self.lines):
            if line.startswith(string):
                return i
        raise LineNotFound(f"No line starting with {string!r} was found.")

    def _locate_slice_for_state(self, state: str) -> slice:
        start, end = self._find_start_and_end_indexes_for_state(state)
        return slice(start, end)

    def _find_start_and_end_indexes_for_state(self, state: str) -> Tuple[int, int]:
        for i, line in enumerate(self.lines):
            if state in line:
                break
        for j, line in enumerate(self.lines[i:], i):
            if state not in line:
                break
        start, end = i, j
        return start, end

    def multiplicity_of_nuclear_state(self, state: str) -> int:
        return self._nuclear_state_multiplicity[state.upper()]

    def ordinal_of_nuclear_state(self, state: str) -> int:
        return self._ordinal_for_nuclear_state[state.upper()]

    def coefficients_for_state(self, state: str) -> List[float]:
        slicer = self.coefficient_slices[
            state
        ]  # TODO: adopt consistent naming for slicers / slices
        lines = self.lines[slicer]
        coefficients = []
        for line in lines:
            try:
                line_label, term_symbol, coeff, *rest = parse_lsj_line(line)
            except ValueError as error:
                msg, *rest = error.args
                assert msg.startswith("not enough values to unpack")
                # Indeed, some lines in the input may not have coefficients listed.
            else:
                coefficients.append(float(coeff))
        return list(coefficients)

    def cutoff_energy_for_state(self, state: str) -> float:
        slicer = self.coefficient_slices[state]
        lines = self.lines[slicer]
        cutoff_energies = set()
        for line in lines:
            label, state, coeff, decay, cutoff = parse_lsj_line(line)
            cutoff_energies.add(cutoff)
        if len(cutoff_energies) != 1:
            raise MultipleCutoffEnergiesForState(f"{state}, {list(cutoff_energies)}")
        return cutoff_energies.pop()


def parse_lsj_line(line: str):
    label, state, *rest = line.split()
    coefficient, decay_constant, cutoff_energy = [float(x) for x in rest]
    return label, state, coefficient, decay_constant, cutoff_energy


class InputReaderException(Exception):
    """Base exception for this module."""


class LineNotFound(InputReaderException):
    """Could not locate line within input file."""


class MultipleCutoffEnergiesForState(InputReaderException):
    """The given state has more than one value of cutoff energy."""
