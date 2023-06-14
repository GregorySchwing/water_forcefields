"""Create atomistic representation of ethanol."""
import os

import mbuild as mb

from reproducibility_project.src import molecules


class a99SB_disp(mb.Compound):
    """Create a single particle water compound."""

    def __init__(self):
        super(a99SB_disp, self).__init__()
        abs_path = os.path.dirname(os.path.abspath(molecules.__file__))
        self.add(mb.load(f"{abs_path}/a99SB_disp.mol2"), label="WAT")


def main():
    """Create a a99SB_disp compound and print basic properties."""
    water = a99SB_disp()
    print(water)
    print(water.name)
    print(water.labels)
    print(water["WAT"])


if __name__ == "__main__":
    main()
