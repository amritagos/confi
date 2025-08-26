from collections import Counter
import re
import pytest
import confi
from pathlib import Path
import subprocess
import numpy as np
from ase import Atoms, Atom
from ase.io import read, lammpsdata

import confi.parameters
import confi.io
import confi.misc


def virtual_site_position(o_pos, h1_pos, h2_pos, om_distance=0.1546):
    # Order: OHHM,OHHM,...
    # DOI: 10.1002/(SICI)1096-987X(199906)20:8
    n = (h1_pos + h2_pos) / 2 - o_pos
    n /= np.linalg.norm(n)
    m_pos = o_pos + om_distance * n
    return m_pos


def test_convert_data_to_g96():
    """Using the PACKMOL generated PDB file as input, create a LAMMPS data file"""
    package_root = Path(__file__).parent.parent.resolve()  # top-level directory
    output_g96 = package_root / Path(
        "tests/output/conversion/system.g96"
    )  # relative to the top-level directory

    output_g96.parent.mkdir(parents=True, exist_ok=True)
    input_data = package_root / Path("tests/resources/system.data")

    Z_of_type = {1: 8, 2: 1, 3: 26, 4: 17}
    # Read in the LAMMPS data file
    with open(input_data, "r") as f:
        atoms_no_massless_sites = lammpsdata.read_lammps_data(
            f, Z_of_type=Z_of_type, units="real", atom_style="full"
        )

    molecule_ids_no_massless_sites = atoms_no_massless_sites.get_array("mol-id")

    atoms = Atoms(cell=atoms_no_massless_sites.cell, pbc=atoms_no_massless_sites.pbc)
    residuenumbers = []
    residuenames = []
    atomtypes = []

    for i, (atom_prev, mol_id_prev) in enumerate(
        zip(atoms_no_massless_sites, molecule_ids_no_massless_sites)
    ):
        if atom_prev.symbol == "H":
            continue
        # O, Fe, Cl
        atoms.append(atom_prev)
        residuenumbers.append(mol_id_prev)
        if atom_prev.symbol == "O":
            # also add the Hs
            atoms.append(atoms_no_massless_sites[i + 1])
            atoms.append(atoms_no_massless_sites[i + 2])
            residuenumbers.append(molecule_ids_no_massless_sites[i + 1])
            residuenumbers.append(molecule_ids_no_massless_sites[i + 2])
            # Get position of the massless site
            o_pos = atoms_no_massless_sites.positions[i]
            h1_pos = atoms_no_massless_sites.positions[i + 1]
            h2_pos = atoms_no_massless_sites.positions[i + 2]
            m_pos = virtual_site_position(o_pos, h1_pos, h2_pos, om_distance=0.1546)
            atoms.append(
                Atom(
                    symbol="X", position=m_pos, momentum=(0.0, 0.0, 0.0), charge=-1.1128
                )
            )
            residuenumbers.append(mol_id_prev)

    # g96 expects residuenames, atomtypes and molecule ids

    atomtypes = [atom.symbol for atom in atoms]
    residuenames = np.array([atom.symbol for atom in atoms])

    # Skip a particular atom wh"ose atomtype and residuename has been set
    skip = np.array([False for atom in atoms])

    for i, atom in enumerate(atoms):
        if skip[i]:
            continue
        if atom.symbol == "Fe":
            mol1 = residuenumbers[i]  # molecule ID
            mol2 = residuenumbers[i + 1]
            mol3 = residuenumbers[i + 2]
            if mol1 == mol2 == mol3:  # dimer
                atomtypes[i] = "FeD"
                atomtypes[i + 1] = "ClD"
                atomtypes[i + 2] = "ClD"
                residuenames[i : i + 3] = "DIM"
                skip[i : i + 3] = True
            elif mol1 == mol2:  # monomer
                atomtypes[i] = "FeM"
                atomtypes[i + 1] = "ClM"
                residuenames[i : i + 2] = "MON"
                skip[i : i + 2] = True
            else:
                atomtypes[i] = atom.symbol
        elif atom.symbol == "Cl":
            atomtypes[i] = atom.symbol
        if atom.symbol == "O":
            atomtypes[i] = "OW1"
            atomtypes[i + 1] = "HW2"
            atomtypes[i + 2] = "HW3"
            atomtypes[i + 3] = "MW4"
            residuenames[i : i + 4] = "WAT"
            skip[i : i + 4] = True

    atoms.arrays["atomtypes"] = atomtypes
    atoms.arrays["residuenames"] = residuenames
    atoms.arrays["residuenumbers"] = residuenumbers

    # Write this out to a g96 file
    with open(output_g96, "w") as f:
        confi.io.write_gromos(
            f,
            atoms,
            write_velocities=True,
        )
