"""Create atomistic representation of ethanol."""
import os

import mbuild as mb

from reproducibility_project.src import molecules


class tip3p(mb.Compound):
    """Create a single particle water compound."""

    def __init__(self):
        super(tip3p, self).__init__()
        abs_path = os.path.dirname(os.path.abspath(molecules.__file__))
        self.add(mb.load(f"{abs_path}/tip3p.mol2"), label="WAT")


def main():
    """Create a tip3p compound and print basic properties."""
    water = tip3p()
    print(water)
    print(water.name)
    print(water.labels)
    print(water["WAT"])


if __name__ == "__main__":
    main()
