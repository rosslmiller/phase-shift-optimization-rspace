
## What's the problem?

As we work towards adapting our previous project's code base to the new
input file format, it's imperative that we determine which API elements
should be kept and which should be discarded.

For the [momentum space project][momentum], classes such as `InputReader` and
`InputWriter` were developed to accommodate two different types of user
interaction, specifically
    1. direct modification of the input file `dph51.d` (i.e., with `vim`,
       `nano`, or sublime), and
    2. via the command line, e.g.
```
$ python lm_run.py --nuclear-states 1s0 --coefficients -0.15 2.45 3.41 -15.21
```

However, for users who never use the command line interface, these classes
contain a significant amount of extraneous features, making them harder to
understand. At this time, I would like to determine which aspects of these
classes should be kept if the intention is to accommodate only the first type
of user interaction.

[momentum]: https://github.com/pricebenjamin/phase-shift-analysis

## How do we interact with `InputWriter` and `OutputReader` instances?

```python
    # lm_run:main; lines 21-30
    input_writer = InputWriter(PHASE_SHIFT_INPUT_FILE)

    if not args.test_case:
        modify_input_file(input_writer, args)

    run_phase_shift_executable()

    output_reader = OutputReader(PHASE_SHIFT_OUTPUT_FILE)
    results = output_reader.get_results_from_output()
    input_writer.determine_modifiable_states_from_output(results)
```

```python
    # lm_run:main; line 46
    initial_coefficients = input_writer.initial_coefficients
```

```python
    # lm_run:compute_phase_shifts; lines 155-158
    input_writer.modify_coefficients(coefficients)
    input_writer.write_lines()
    run_phase_shift_executable()
    results = output_reader.get_results_from_output()
```

## How might you like to interact with these instances?

```python
input_writer = InputWriter("dphqq.d")

results = output_reader.get_results()
target_states = list(results)
input_writer.target_states = target_states
```

## `InputWriter` API

```python
class InputWriter:

    @property
    def states_to_modify(self):        # TODO(ben): consider better name
        return self._states_to_modify

    @states_to_modify.setter
    def states_to_modify(self, states: List[str]):
        # Verify and store list of states

    def __init__(self, filename):
        # Read the file
        # Locate and parse partial wave coefficients

    def modify_coefficients(self, coefficients: List[float]):
        # Responsible for interpreting the list of coefficients
        # and injecting them into the proper locations within
        # the input file

    def write_lines(self):
        # Overwrite the input file with new contents
```

## `OutputReader` API

```python
class OutputReader:

    def __init__(self, filename: str):
        # Not much to initialize

    def fetch_results(self) -> Dict[str, StateResults]:
        # Read and store results for all available states
```

## Further questions

#### Should we separate duties between an `InputReader` and an `InputWriter`?
