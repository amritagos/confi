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


def get_neighbour_list(
    atoms: Atoms,
    cation_symbol: str,
    anion_symbol: str,
    cutoff: float,
    tolerance: float = 0.01,
):
    neigh_list = {}
    # cation-anion pairs
    for i in range(len(atoms) - 1):
        if atoms[i].symbol != cation_symbol:
            continue
        neigh_i = []
        for j in range(len(atoms)):
            if atoms[j].symbol == anion_symbol:
                # check cutoff
                r_ij = atoms.get_distance(i, j, mic=True)
                if r_ij <= cutoff + tolerance:
                    neigh_i.append(j)
        if len(neigh_i) > 0:
            neigh_list[i] = neigh_i
    return neigh_list


def find_ion_distribution(
    atoms: Atoms, cation_symbol: str, anion_symbol: str, cutoff: float
):
    neigh_list = get_neighbour_list(atoms, cation_symbol, anion_symbol, cutoff)

    free_cations = []
    free_anions = []
    monomers = (
        []
    )  # cation anion pairs corresponding to the species cation-anion (list of indices in the Atoms object)
    dimers = (
        []
    )  # cation-anion-anion groups corresponding to the species anion-cation-anion (list of lists of indices in the Atoms object)
    visited_anions = set()

    for i, atom in enumerate(atoms):
        if atom.symbol == cation_symbol:
            if i in neigh_list:
                # i has Cl nearby
                n_neigh_i = len(neigh_list[i])
                # for monomer pairs
                if n_neigh_i == 1:
                    anion_idx = neigh_list[i][0]
                    pair = [i, anion_idx]
                    monomers.append(pair)
                    visited_anions.add(anion_idx)
                # for dimers
                if n_neigh_i == 2:
                    anion_idx = neigh_list[i][0]
                    visited_anions.add(anion_idx)
                    anion_idx1 = neigh_list[i][1]
                    visited_anions.add(anion_idx1)
                    dimer_group = [i, anion_idx, anion_idx1]
                    dimers.append(dimer_group)
            else:
                free_cations.append(i)

    # Find the free anions
    for i, atom in enumerate(atoms):
        if atom.symbol == anion_symbol:
            if atom.index not in visited_anions:
                free_anions.append(atom.index)
                visited_anions.add(atom.index)

    return free_cations, free_anions, monomers, dimers


def set_atom_arrays(atoms: Atoms, atomtypes, residuenames, residuenumbers):
    residuenames = np.array(residuenames)
    atomtypes = np.array(atomtypes)
    residuenumbers = np.array(residuenumbers)
    atoms.set_array("residuenames", residuenames)
    atoms.set_array("atomtypes", atomtypes)
    atoms.set_array("residuenumbers", residuenumbers)


def get_water_molecules(
    atoms_unordered: Atoms,
    atoms: Atoms,
    ase_to_gro_dict: dict[str:str],
    atomtypes: list[str],
    residuenames: list[str],
    residuenumbers: list[int],
):
    """From an input Atoms object (atoms_unordered), add water molecules in a specific way into the output atoms object (atoms).
    The lists atomtypes, residuenames and residuenumbers are updated as molecules are added to the ASE atoms object.

    Args:
        atoms_unordered (Atoms): input Atoms object which may contain several types of molecules and species
        atoms (Atoms): Water molecules will be added into this ASE atoms object
        ase_to_gro_dict (dict): ASE atom symbols are keys and the GROMACS symbol is the value (does not have to correspond to an actual atom symbol)
        atomtypes (list[str]): GROMACS symbols for each atom
        residuenames (list[str]): Name of the residue in GROMACS
        residuenumbers (list[int]): Molecule ID of the atom (all atoms belonging to the same molecule have the same molecule ID)
    """
    atom_count = 1
    # Get the water molecules first
    for atom in atoms_unordered:
        if atom.symbol == "O":
            o_idx = atom.index
            atoms.append(atoms_unordered[o_idx])  # O
            atoms.append(atoms_unordered[o_idx + 1])  # H1
            atoms.append(atoms_unordered[o_idx + 2])  # H2
            atoms.append(atoms_unordered[o_idx + 3])  # M
            # if we get oxygen, write out O H H M
            atomtypes.append(ase_to_gro_dict[atom.symbol][0])  # O
            atomtypes.append(ase_to_gro_dict["H"][0])  # H1
            atomtypes.append(ase_to_gro_dict["H"][1])  # H2
            atomtypes.append(ase_to_gro_dict["X"][0])  # M
            for i in range(4):
                residuenames.append(f"water")
                residuenumbers.append(atom_count)
            atom_count += 1


def get_ions(
    atoms_unordered: Atoms,
    atoms: Atoms,
    indices: list[int],
    ase_to_gro_dict: dict[str:str],
    atomtypes: list[str],
    residuenames: list[str],
    residuenumbers: list[int],
):
    atom_count = len(atoms)
    for i in indices:
        atom = atoms_unordered[i]
        atoms.append(atom)
        atomtypes.append(ase_to_gro_dict[atom.symbol][0])
        residuenames.append(f"{ase_to_gro_dict[atom.symbol][0]}")
        residuenumbers.append(atom_count)
        atom_count += 1


def get_monomer(
    atoms_unordered: Atoms,
    atoms: Atoms,
    monomers: list[list[int]],
    residuename: str,
    ase_to_gro_dict: dict[str:str],
    atomtypes: list[str],
    residuenames: list[str],
    residuenumbers: list[int],
):
    """From an input Atoms object (atoms_unordered), add monomer units in a specific way into the output atoms object (atoms).
    The lists atomtypes, residuenames and residuenumbers are updated as "molecules" or units are added to the ASE atoms object.

    Args:
        atoms_unordered (Atoms): Atoms object which may contain many species
        atoms (Atoms): Atoms object into which monomer units are added to (output)
        monomers (list[list[int]]): list of list of indices corresponding to monomer units. The indices are wrt to atoms_unordered
        residuename (str): Name of the residue of each monomer unit
        ase_to_gro_dict (dict[str:str]): ASE atom symbols are keys and the GROMACS symbol is the value (does not have to correspond to an actual atom symbol)
        atomtypes (list[str]):  GROMACS symbols for each atom
        residuenames (list[str]): Name of the residue in GROMACS
        residuenumbers (list[int]): Molecule ID of each atom (all atoms belonging to the same molecule have the same molecule ID)
    """
    try:
        molecule_count = residuenumbers[-1] + 1
    except:
        molecule_count = len(atoms)

    for i_monomer in monomers:
        i_cat = i_monomer[0]
        i_anion = i_monomer[1]
        atoms.append(atoms_unordered[i_cat])  # FeM
        atoms.append(atoms_unordered[i_anion])  # ClM
        atomtypes.append(ase_to_gro_dict[atoms_unordered[i_cat].symbol][1])  # FeM
        atomtypes.append(ase_to_gro_dict[atoms_unordered[i_anion].symbol][1])  # ClM
        for i in range(2):
            residuenames.append(f"{residuename}")
            residuenumbers.append(molecule_count)
        molecule_count += 1


def get_dimer(
    atoms_unordered: Atoms,
    atoms: Atoms,
    dimers: list[list[int]],
    residuename: str,
    ase_to_gro_dict: dict[str:str],
    atomtypes: list[str],
    residuenames: list[str],
    residuenumbers: list[int],
):
    """From an input Atoms object (atoms_unordered), add dimer units in a specific way into the output atoms object (atoms).
    The lists atomtypes, residuenames and residuenumbers are updated as "molecules" or units are added to the ASE atoms object.

    Args:
        atoms_unordered (Atoms): Atoms object which may contain many species
        atoms (Atoms): Atoms object into which monomer units are added to (output)
        dimers (list[list[int]]): list of list of indices corresponding to dimer units. The indices are wrt to atoms_unordered
        residuename (str): Name of the residue of each dimer unit
        ase_to_gro_dict (dict[str:str]): ASE atom symbols are keys and the GROMACS symbol is the value (does not have to correspond to an actual atom symbol)
        atomtypes (list[str]):  GROMACS symbols for each atom
        residuenames (list[str]): Name of the residue in GROMACS
        residuenumbers (list[int]): Molecule ID of each atom (all atoms belonging to the same molecule have the same molecule ID)
    """
    try:
        molecule_count = residuenumbers[-1] + 1
    except:
        molecule_count = len(atoms)

    for i_dimer in dimers:
        i_cat = i_dimer[0]
        i_anion0 = i_dimer[1]
        i_anion1 = i_dimer[2]
        atoms.append(atoms_unordered[i_cat])  # FeD
        atoms.append(atoms_unordered[i_anion0])  # ClD
        atoms.append(atoms_unordered[i_anion1])  # ClD
        atomtypes.append(ase_to_gro_dict[atoms_unordered[i_cat].symbol][2])  # FeD
        atomtypes.append(ase_to_gro_dict[atoms_unordered[i_anion0].symbol][2])  # ClD
        atomtypes.append(ase_to_gro_dict[atoms_unordered[i_anion1].symbol][2])  # ClD
        for i in range(3):
            residuenames.append(f"{residuename}")
            residuenumbers.append(molecule_count)
        molecule_count += 1


def test_gromacs_input(test_packmol_input_gromacs):
    """Using the PACKMOL generated XYZ file as input, create a LAMMPS data file"""
    package_root = Path(__file__).parent.parent.resolve()  # top-level directory
    output_g96 = package_root / Path(
        "tests/output/gromacs/system.g96"
    )  # relative to the top-level directory
    cation_symbol_packmol = "Fe"  # symbol in the PACKMOL generated XYZ file
    anion_symbol_packmol = "Cl"
    cation_symbol_gromacs = "Fe"  # symbol in the output g96 file
    cutoff = 2.33  # distance between cation and anion

    output_g96.parent.mkdir(parents=True, exist_ok=True)

    packmol_params = test_packmol_input_gromacs

    # Dictionary from the PackmolParams class
    packmol_param_dict = packmol_params.model_dump()

    # Read in the XYZ file
    atoms_unordered = read(packmol_param_dict["system_file"])
    atoms = Atoms()
    atoms.set_cell(
        [
            packmol_param_dict["x_box_length"],
            packmol_param_dict["y_box_length"],
            packmol_param_dict["z_box_length"],
        ]
    )
    atoms.pbc = [True, True, True]
    ase_to_gro_dict = {
        "O": ["OW1"],
        "H": ["HW2", "HW3"],
        "X": ["MW4"],
        cation_symbol_packmol: [
            cation_symbol_gromacs,
            cation_symbol_gromacs + "M",
            cation_symbol_gromacs + "D",
        ],
        "Cl": ["Cl", "ClM", "ClD"],
    }

    # Things for writing out GROMOS files
    atomtypes = []
    residuenames = []
    residuenumbers = []

    get_water_molecules(
        atoms_unordered, atoms, ase_to_gro_dict, atomtypes, residuenames, residuenumbers
    )

    # Now get the free cations, free anions, bound cations and bound anions
    free_cations, free_anions, monomers, dimers = find_ion_distribution(
        atoms_unordered, cation_symbol_packmol, anion_symbol_packmol, cutoff
    )

    n_free_cations = len(free_cations)
    n_free_anions = len(free_anions)
    n_monomers = len(monomers)
    n_dimers = len(dimers)

    assert n_free_cations == packmol_param_dict["n_free_cations"]
    assert n_free_anions == packmol_param_dict["n_free_anions"]
    assert n_monomers == packmol_param_dict["n_monomer"]
    assert n_dimers == packmol_param_dict["n_dimer"]

    # Free cations
    get_ions(
        atoms_unordered,
        atoms,
        free_cations,
        ase_to_gro_dict,
        atomtypes,
        residuenames,
        residuenumbers,
    )
    # Free anions
    get_ions(
        atoms_unordered,
        atoms,
        free_anions,
        ase_to_gro_dict,
        atomtypes,
        residuenames,
        residuenumbers,
    )
    # Fe-Cl monomer
    name_monomer = "FeClM"
    get_monomer(
        atoms_unordered,
        atoms,
        monomers,
        name_monomer,
        ase_to_gro_dict,
        atomtypes,
        residuenames,
        residuenumbers,
    )

    # Cl-Fe-Cl dimer
    name_dimer = "FeClD"
    get_dimer(
        atoms_unordered,
        atoms,
        dimers,
        name_dimer,
        ase_to_gro_dict,
        atomtypes,
        residuenames,
        residuenumbers,
    )

    # Set the arrays
    set_atom_arrays(atoms, atomtypes, residuenames, residuenumbers)
    # Write out the gromacs files
    with open(output_g96, "w") as f:
        confi.io.write_gromos(f, atoms)
