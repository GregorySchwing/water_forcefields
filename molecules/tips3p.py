"""Create atomistic representation of ethanol."""
import os

import mbuild as mb

from reproducibility_project.src import molecules


class tips3p(mb.Compound):
    """Create a single particle water compound."""

    def __init__(self):
        super(tips3p, self).__init__()
        abs_path = os.path.dirname(os.path.abspath(molecules.__file__))
        self.add(mb.load(f"{abs_path}/tips3p.mol2"), label="WAT")


def main():
    """Create a tisp3p compound and print basic properties."""
    water = tips3p()
    print(water)
    print(water.name)
    print(water.labels)
    print(water["WAT"])


if __name__ == "__main__":
    main()
