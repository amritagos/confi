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
import confi.misc


def create_bonds_angles(atoms: Atoms):
    # Create the bonds and angles arrays
    # bonds array inner loop bond_type atom1 atom2
    # angles array inner loop angle_type atom1 atom2 atom3
    # bond type 1 is water
    # bond type 2 is for monomers and dimers
    # angle type 1 is water
    bonds_in = []
    angles_in = []
    ids = []  # Atom IDs; index + 1

    residuenames = atoms.get_array("residuenames")
    atomtypes = atoms.get_array("atomtypes")

    for i, (atom, atomtype, residuename) in enumerate(
        zip(atoms, atomtypes, residuenames)
    ):
        ids.append(i + 1)  # ID = index + 1

        if residuename == "MON" and atomtype == "FeM":
            bonds_in.append([2, i + 1, i + 2])
        elif residuename == "DIM" and atomtype == "FeD":
            bonds_in.append([2, i + 1, i + 2])
            bonds_in.append([2, i + 1, i + 3])
        elif residuename == "WAT" and atomtype == "OW1":
            bonds_in.append([1, i + 1, i + 2])
            bonds_in.append([1, i + 1, i + 3])
            angles_in.append([1, i + 2, i + 1, i + 3])

    return bonds_in, angles_in


def test_convert_pdb(test_packmol_input_gromacs_pdb):
    """Using the PACKMOL generated PDB file as input, create a LAMMPS data file"""
    package_root = Path(__file__).parent.parent.resolve()  # top-level directory
    output_data = package_root / Path(
        "tests/output/conversion/system_from_pdb.data"
    )  # relative to the top-level directory

    output_data.parent.mkdir(parents=True, exist_ok=True)

    packmol_params = test_packmol_input_gromacs_pdb

    # Dictionary from the PackmolParams class
    packmol_param_dict = packmol_params.model_dump()

    # Read in the PDB file
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

    # Delete the massless fourth site
    del atoms[
        [
            atom.index
            for i, atom in enumerate(atoms)
            if atoms.get_array("atomtypes")[i] == "MW4"
        ]
    ]

    # dictionary going from atom type (from PDB) to the ASE symbol
    Z_atomtype_to_symbol = {
        "FeM": "Fe",
        "ClM": "Cl",
        "FeD": "Fe",
        "ClD": "Cl",
        "Fe": "Fe",
        "Cl": "Cl",
        "HW2": "H",
        "HW3": "H",
        "OW1": "O",
    }

    bonds_in, angles_in = create_bonds_angles(atoms)
    bonds = confi.misc.convert_bond_array(bonds_in, len(atoms))
    angles = confi.misc.convert_angle_array(angles_in, len(atoms))
    # Set bond and angle arrays
    atoms.arrays["bonds"] = bonds
    atoms.arrays["angles"] = angles

    # Set the correct symbols in the atoms object, according to the atomtypes object
    for i, atom in enumerate(atoms):
        atomtype = atoms.get_array("atomtypes")[i]
        atom.symbol = Z_atomtype_to_symbol[atomtype]

    # The molecule IDs are residuenumbers
    atoms.arrays["mol-id"] = atoms.get_array("residuenumbers")

    # uniform random numbers between -0.1 and 0.1 Å/fs
    vels = np.random.uniform(-0.1, 0.1, size=(len(atoms), 3))
    atoms.set_velocities(vels)

    # Write out the LAMMPS data file
    with open(output_data, "w") as f:
        confi.io.write_lammps_data(
            f,
            atoms,
            specorder=["O", "H", "Fe", "Cl"],
            atom_style="full",
            velocities=True,
        )
