from collections import Counter
import re
import pytest
import confi
from pathlib import Path
import subprocess
import numpy as np
from ase import Atoms
from ase.io import read

import confi.parameters
import confi.io


def renumber_residuenumbers(atoms: Atoms):
    # Extract the residue numbers which are not numbered consistently
    orig_residuenumbers = atoms.get_array("residuenumbers")
    residuenames = atoms.get_array("residuenames")

    residuenumbers = atoms.get_array("residuenumbers")
    counter = 1
    res_map = {}

    for i, (chain, old) in enumerate(zip(residuenames, orig_residuenumbers)):
        key = (chain, old)
        if key not in res_map:
            res_map[key] = counter
            counter += 1
        residuenumbers[i] = res_map[key]

    # 4. Overwrite the ASE atom‐array and write out a new PDB
    atoms.set_array("residuenumbers", residuenumbers)


def test_gromacs_input(test_packmol_input_gromacs_pdb):
    """Using the PACKMOL generated XYZ file as input, create a LAMMPS data file"""
    package_root = Path(__file__).parent.parent.resolve()  # top-level directory
    output_g96 = package_root / Path(
        "tests/output/gromacs/system_from_pdb.g96"
    )  # relative to the top-level directory

    output_g96.parent.mkdir(parents=True, exist_ok=True)

    packmol_params = test_packmol_input_gromacs_pdb

    # Dictionary from the PackmolParams class
    packmol_param_dict = packmol_params.model_dump()

    # Read in the XYZ file
    atoms = read(packmol_param_dict["system_file"], format="proteindatabank")
    # Update the box lengths
    atoms.set_cell(
        [
            packmol_param_dict["x_box_length"],
            packmol_param_dict["y_box_length"],
            packmol_param_dict["z_box_length"],
        ]
    )
    atoms.set_pbc([True, True, True])
    # # Renumber the molecule IDs or residuenumbers
    # renumber_residuenumbers(atoms)

    # Write out the gromacs files
    with open(output_g96, "w") as f:
        confi.io.write_gromos(f, atoms)
