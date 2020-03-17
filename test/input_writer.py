import random
import time

from typing import List, Tuple, Dict

# Local imports
from input_reader import InputReader, parse_lsj_line
from output_reader import StateResults


MAX_KINETIC_ENERGY_WIDTH = 4
MAX_COEFFICIENT_WIDTH = 13

RANDOM_UNIFORM_LOWER_BOUND = 0
RANDOM_UNIFORM_UPPER_BOUND = 1


class NuclearState:
    orbitals = "SPDFGHIJKLMNOPQRSTUVWXYZ"


class InputWriter(InputReader):

    _coupled_states = set(
        [
            "3P0",
            "3S1",
            "3D1",
            "3SD1",
            "3P2",
            "3F2",
            "3PF2",
            "3D3",
            "3G3",
            "3D-G3",
            "3F4",
            "3H4",
            "3F-H4",
        ]
    )

    @property
    def nuclear_states_to_modify(self) -> List[str]:
        return list(self._nuclear_states_to_modify)

    @nuclear_states_to_modify.setter
    def nuclear_states_to_modify(self, states: List[str]):
        nuclear_states_to_modify = list(s.upper() for s in states)
        self._validate_nuclear_states(nuclear_states_to_modify)
        self._nuclear_states_to_modify = nuclear_states_to_modify

    @property
    def initial_coefficients(self) -> List[float]:
        # TODO(ben): consider rename initial_coefficients_of_states_to_modify
        coefficients = []
        for state in self._nuclear_states_to_modify:
            coefficients.extend(self.coefficients_for_state(state))
        return coefficients

    @InputReader.J_start.setter
    def J_start(self, value: int):
        assert value >= 0
        self._modify_J_range(value, self.J_end)
        self._read_J_range()

    @InputReader.J_end.setter
    def J_end(self, value: int):
        assert value >= 0
        self._modify_J_range(self.J_start, value)
        self._read_J_range()

    # TODO(ben): Consider new property J_range
    # Allows us to set J_start and J_end simultaneously, which means
    # we can immediately verify that 0 <= J_start <= J_end

    @InputReader.singlet_states_enabled.setter
    def singlet_states_enabled(self, is_enabled: bool):
        self._modify_multiplicity_selector(
            singlet=is_enabled,
            triplet=self.triplet_states_enabled,
            coupled=self.coupled_states_enabled,
        )
        self._read_multiplicity_selector()

    @InputReader.triplet_states_enabled.setter
    def triplet_states_enabled(self, is_enabled: bool):
        self._modify_multiplicity_selector(
            singlet=self.singlet_states_enabled,
            triplet=is_enabled,
            coupled=self.coupled_states_enabled,
        )
        self._read_multiplicity_selector()

    @InputReader.coupled_states_enabled.setter
    def coupled_states_enabled(self, is_enabled: bool):
        self._modify_multiplicity_selector(
            singlet=self.singlet_states_enabled,
            triplet=self.triplet_states_enabled,
            coupled=is_enabled,
        )
        self._read_multiplicity_selector()

    @InputReader.kinetic_energies_in_MeV.setter
    def kinetic_energies_in_MeV(self, energies: List[float]):
        energies.sort()  # May not be required...
        self._modify_kinetic_energies(energies)
        # Since this modification may change the number of lines in self.lines,
        # we need to re-locate all parameters below this section in the input file.
        self._locate_kinetic_energies_slice()
        self._locate_lsj_slice()
        self._locate_coefficients()
        self._read_kinetic_energies()

    def __init__(self, filename: str, output_filename=None):
        # Note: output_filename is used only for testing
        super().__init__(filename)
        if output_filename is None:
            self._output_filename = self._filename
        else:
            self._output_filename = output_filename
        set_random_seed_to_current_time_in_ns()

    def _validate_nuclear_states(self, states: List[str]):
        for state in states:
            assert state in self._set_of_nuclear_states
            self._check_state_matches_input_parameters(state)

    def _check_state_matches_input_parameters(self, state: str):
        # We already checked that state is in self._set_of_nuclear_states
        multiplicity, *orbitals, J_value = list(state)
        multiplicity = int(multiplicity)
        J_value = int(J_value)
        self._validate_multiplicity(multiplicity, state)
        # TODO(ben): Why does 3P0 show up when requesting coupled parameters?
        self._validate_orbitals(orbitals, state)
        self._validate_J_value(J_value, state)

    def _validate_multiplicity(self, multiplicity: int, state: str):
        if state in self._coupled_states and self.coupled_states_enabled is True:
            return

        if multiplicity == 1 and self.singlet_states_enabled is False:
            raise InvalidMultiplicity(
                f"multiplicity == 1 for state {state} "
                "but self.singlet_states_enabled is False"
            )

        if multiplicity == 3 and self.triplet_states_enabled is False:
            raise InvalidMultiplicity(
                f"multiplicity == 3 for state {state} "
                "but self.triplet_states_enabled is False"
            )

    def _validate_orbitals(self, orbitals: List[int], state: str):
        if len(orbitals) > 1 and self.coupled_states_enabled is False:
            raise CoupledStateError(
                f"{state} is a coupled state "
                "but self.coupled_states_enabled is False"
            )

    def _validate_J_value(self, J_value: int, state: str):
        if not self.J_start <= J_value <= self.J_end:
            raise InvalidJValue(
                f"J_value must be in [{self.J_start}, {self.J_end}]; "
                f"J_value == {J_value} for state {state}"
            )

    def _modify_J_range(self, start: int, end: int):
        line = f"jb,je       {start}  {end}\n"
        self.lines[self.J_range_index] = line

    def _modify_multiplicity_selector(self, singlet, triplet, coupled):
        line = f"sing,tr,cp  {int(singlet)}  {int(triplet)}  {int(coupled)}\n"
        self.lines[self.multiplicity_selector_index] = line

    def _modify_kinetic_energies(self, energies: List[float]):
        new_lines = []
        for energy in energies:
            energy_string = f"{energy:.4f}"
            energy_string = energy_string[:MAX_KINETIC_ENERGY_WIDTH]
            assert "." in energy_string
            line = f"          {energy_string}\n"
            new_lines.append(line)
        self.lines[self.kinetic_energies_slice] = new_lines

    def multiplicity_selector(self, singlet, triplet, coupled):
        self._modify_multiplicity_selector(singlet, triplet, coupled)
        self._read_multiplicity_selector()

    def modify_coefficients(self, vector: List[float]):
        start = 0
        for state in self._nuclear_states_to_modify:
            multiplicity = self.multiplicity_of_nuclear_state(state)
            end = start + multiplicity
            coefficients = vector[start:end]
            self.modify_nuclear_state_coefficients(state, coefficients)
            start = end

    def modify_nuclear_state_coefficients(self, state: str, coefficients: List[float]):
        slicer = self.coefficient_slices[state]
        old_lines = self.lines[slicer]
        new_lines = []
        for coeff, old_line in zip(coefficients, old_lines):
            label, state, old_coeff, decay, cutoff = parse_lsj_line(old_line)
            line = create_lsj_line(label, state, coeff, decay, cutoff)
            new_lines.append(line)
        self.lines[slicer] = new_lines

    def determine_modifiable_states_from_output(self, results: Dict[str, StateResults]):
        output_states = list(results)
        modifiable_states = self._translate_to_input_states(output_states)
        self.nuclear_states_to_modify = modifiable_states

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

    def create_random_coefficients(
        self, lower=RANDOM_UNIFORM_LOWER_BOUND, upper=RANDOM_UNIFORM_UPPER_BOUND
    ):
        coefficients = []
        for state in self.nuclear_states_to_modify:
            m = self.multiplicity_of_nuclear_state(state)
            coefficients.extend([random.uniform(lower, upper) for _ in range(m)])
        return coefficients

    def modify_cutoff_energy(self, cutoff_energy: float):
        for state in self.nuclear_states_to_modify:
            self._modify_nuclear_state_cutoff_energies(state, cutoff_energy)

    def _modify_nuclear_state_cutoff_energies(self, state: str, cutoff_energy: float):
        slicer = self.coefficient_slices[state]
        old_lines = self.lines[slicer]
        new_lines = []
        for old_line in old_lines:
            # TODO(ben): It's hard to see exactly what changes between the following
            # two lines. Make the change more apparent.
            label, state, coeff, decay, old_cutoff = parse_lsj_line(old_line)
            new_line = create_lsj_line(label, state, coeff, decay, cutoff_energy)
            new_lines.append(new_line)
        self.lines[slicer] = new_lines

    def write_lines(self):
        assert self.J_start <= self.J_end
        # TODO(ben): consider J_range property
        with open(self._output_filename, "w") as f:
            f.writelines(self.lines)


def set_random_seed_to_current_time_in_ns():
    seed = time.time_ns()
    random.seed(seed)


def create_lsj_line(
    label: str,
    state: str,
    coefficient: float,
    decay_constant: float,
    cutoff_energy: float,
):
    coefficient = truncate_coefficient_to_max_width(coefficient)
    decay_constant = f"{decay_constant: < .1f}"
    cutoff_energy = f"{cutoff_energy: < .1f}"
    return f"{label:<5}{state:<6}{coefficient:<14}{decay_constant:<10}{cutoff_energy}\n"


def truncate_coefficient_to_max_width(coefficient: float) -> str:
    coefficient = f"{coefficient: < .15f}"
    return coefficient[:MAX_COEFFICIENT_WIDTH]


class InputWriterException(Exception):
    """Base exception for this module."""


# _validate_nuclear_states(self, state) exceptions
class InvalidMultiplicity(InputWriterException):
    """Multiplicity did not match multiplicity selectors."""


class CoupledStateError(InputWriterException):
    """A coupled state was provided but coupled states are not enabled."""


class InvalidJValue(InputWriterException):
    """The provided J_value is not within the required range."""
