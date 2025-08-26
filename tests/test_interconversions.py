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


def test_convert_data_to_g96():
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
    residuenames = np.array(["DUM"] * len(atoms), dtype=object)

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
                residuenames[i] = atom.symbol
        elif atom.symbol == "Cl":
            atomtypes[i] = atom.symbol
            residuenames[i] = atom.symbol
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
    velocities_prev = atoms.get_velocities()
    # -----------------------------------------------
    output_data = package_root / Path(
        "tests/output/conversion/system_from_g96.data"
    )  # relative to the top-level directory
    # Read the g96 file into an Atoms object
    with open(output_g96, "r") as f:
        atoms_gromos = confi.io.read_gromos(f, dummy_site_symbols=["MW4"])

    velocities_read_in = atoms_gromos.get_velocities()

    # Velocities should be the same
    # 1) sanity checks
    assert velocities_prev is not None
    assert velocities_read_in is not None
    assert velocities_prev.shape == velocities_read_in.shape

    # 2) numeric comparison with tolerance
    assert np.allclose(velocities_prev, velocities_read_in, rtol=1e-7, atol=1e-9)

    # Delete the massless fourth site
    del atoms_gromos[
        [
            atom.index
            for i, atom in enumerate(atoms_gromos)
            if atoms_gromos.get_array("atomtypes")[i] == "MW4"
        ]
    ]

    # dictionary going from atom type (from the g96 file) to the ASE symbol
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

    bonds_in, angles_in = create_bonds_angles(atoms_gromos)
    bonds = confi.misc.convert_bond_array(bonds_in, len(atoms_gromos))
    angles = confi.misc.convert_angle_array(angles_in, len(atoms_gromos))
    # Set bond and angle arrays
    atoms_gromos.arrays["bonds"] = bonds
    atoms_gromos.arrays["angles"] = angles

    # Set the correct symbols in the atoms object, according to the atomtypes object
    for i, atom in enumerate(atoms_gromos):
        atomtype = atoms_gromos.get_array("atomtypes")[i]
        atom.symbol = Z_atomtype_to_symbol[atomtype]

    # The molecule IDs are residuenumbers
    atoms_gromos.arrays["mol-id"] = atoms_gromos.get_array("residuenumbers").astype(int)

    # Write out the LAMMPS data file
    with open(output_data, "w") as f:
        confi.io.write_lammps_data(
            f,
            atoms_gromos,
            specorder=["O", "H", "Fe", "Cl"],
            atom_style="full",
            velocities=True,
            units="real",
        )
