from ase import Atoms
import numpy as np


def convert_bond_array(bonds_in: list[list], natoms: int):
    """Create an array similar to the kind of array created when bonds are read in by a LAMMPS data file,
    given a list of lists with the bond information

    Args:
        bonds_in (list[list]): List of lists such that each inner list corresponds to the bond_type, atom1 ID, and atom2 ID.
        Each inner list corresponds to the last three numbers in a line of the LAMMPS data file for Bonds.
        natoms (int) : Number of atoms
    """
    # Output array, which will have the information for type atom1 atom2
    bonds = [""] * natoms if len(bonds_in) > 0 else None

    # Convert to the format written by ASE's read_lammps_data
    if bonds is not None:
        for bond_type, at1, at2 in bonds_in:
            i_a1 = (
                at1 - 1
            )  # Assume that the index is one less than the ID in the bonds array
            i_a2 = at2 - 1
            if len(bonds[i_a1]) > 0:
                bonds[i_a1] += ","
            bonds[i_a1] += f"{i_a2:d}({bond_type:d})"
        for i, bond in enumerate(bonds):
            if len(bond) == 0:
                bonds[i] = "_"
        return np.array(bonds)


def convert_angle_array(angles_in: list[list], natoms: int):
    """Create an array similar to the kind of array created when angles are read in by a LAMMPS data file,
    given a list of lists with the angle information

    Args:
        angles_in (list[list]): List of lists such that each inner list corresponds to the angle_type atom1 atom2 atom3.
        Each inner list corresponds to the last four numbers in a line of the LAMMPS data file for Angles.
        natoms (int) : Number of atoms
    """
    # Output array, which will have the information for type atom1 atom2
    angles = [""] * natoms if len(angles_in) > 0 else None

    # Convert to the format written by ASE's read_lammps_data
    if angles is not None:
        for angle_type, at1, at2, at3 in angles_in:
            i_a1 = at1 - 1  # Assume the ID is just 1 greater than the index in atoms
            i_a2 = at2 - 1
            i_a3 = at3 - 1
            if len(angles[i_a2]) > 0:
                angles[i_a2] += ","
            angles[i_a2] += f"{i_a1:d}-{i_a3:d}({angle_type:d})"
        for i, angle in enumerate(angles):
            if len(angle) == 0:
                angles[i] = "_"
        return np.array(angles)
