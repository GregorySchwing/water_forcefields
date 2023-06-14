"""Create atomistic representation of ethanol."""
import os

import mbuild as mb

from reproducibility_project.src import molecules


class opc(mb.Compound):
    """Create a single particle water compound."""

    def __init__(self):
        super(opc, self).__init__()
        abs_path = os.path.dirname(os.path.abspath(molecules.__file__))
        self.add(mb.load(f"{abs_path}/opc.mol2"), label="WAT")


def main():
    """Create a opc compound and print basic properties."""
    water = opc()
    print(water)
    print(water.name)
    print(water.labels)
    print(water["WAT"])


if __name__ == "__main__":
    main()
