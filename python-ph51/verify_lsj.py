import numpy
import textwrap

# Local imports
from filenames import PHASE_SHIFT_INPUT_FILE, PHASE_SHIFT_OUTPUT_FILE


def verify_lsj():
    with open(PHASE_SHIFT_INPUT_FILE, "r") as inp, open(
        PHASE_SHIFT_OUTPUT_FILE, "r"
    ) as out:
        in_lines = inp.readlines()
        out_lines = out.readlines()

    # Jump to "lsj" section of both files; store indexes i, j
    i = index_of_line_starting_with("lsj", in_lines)
    j = index_of_line_starting_with(" lsj", out_lines)  # Note the leading space.

    # Parse each line
    zipped_lines = zip(in_lines[i:], out_lines[j:])
    for k, (in_line, out_line) in enumerate(zipped_lines):

        # Check if we've reached the end of the lsj section
        if in_line.startswith("end"):
            assert out_line.startswith(" end")
            break

        ilabel, istate, *iparams = in_line.split()
        olabel, ostate, *oparams = out_line.split()

        assert istate == ostate

        if iparams:
            # Verify first three columns match
            ix, iy, iz, *i_remainder = iparams
            ox, oy, oz, *o_remainder = oparams

            if not numpy.allclose(
                [float(x) for x in (ix, iy, iz)], [float(x) for x in (ox, oy, oz)]
            ):
                raise InputOutputMismatch(
                    in_line=in_line.strip(),
                    out_line=out_line.strip(),
                    locations=(i + k + 1, j + k + 1)  # List indexes start at 0,
                    # but line numbers start at 1
                )

            # Remaining columns should be zero, if they exist
            assert all(float(x) == 0 for x in i_remainder)
            assert all(float(x) == 0 for x in o_remainder)
        else:
            # If no parameters provided in input, output should be zero.
            assert all(float(x) == 0 for x in oparams)

    print(
        f"`lsj` parameters within {PHASE_SHIFT_INPUT_FILE} and {PHASE_SHIFT_OUTPUT_FILE} match."
    )


def index_of_line_starting_with(string, lines):
    for i, line in enumerate(lines):
        if line.startswith(string):
            break
    assert line.startswith(string)  # TODO(ben): Raise an exception if not found?
    return i


class InputOutputMismatch(Exception):
    """Some values in the input did not match the values in the output."""

    def __init__(self, in_line, out_line, locations):
        iloc, oloc = locations
        msg = textwrap.dedent(
            f"""\
            \n
            Error: input values do not match output values (see below)

            Input (line {iloc}):    {in_line}
            Output (line {oloc}):   {out_line}
        """
        )
        super().__init__(msg)


if __name__ == "__main__":
    verify_lsj()
