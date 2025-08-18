from ase import Atoms, units

import re
import warnings

import numpy as np

from ase.atoms import Atoms
from ase.calculators.lammps import Prism, convert
from ase.data import atomic_masses, atomic_numbers


def _write_masses(fd, atoms: Atoms, species: list, units: str):
    symbols_indices = atoms.symbols.indices()
    fd.write("Masses\n\n")
    for i, s in enumerate(species):
        if s in symbols_indices:
            # Find the first atom of the element `s` and extract its mass
            # Cover by `float` to make a new object for safety
            mass = float(atoms[symbols_indices[s][0]].mass)
        else:
            # Fetch from ASE data if the element `s` is not in the system
            mass = atomic_masses[atomic_numbers[s]]
        # Convert mass from ASE units to LAMMPS units
        mass = convert(mass, "mass", "ASE", units)
        atom_type = i + 1
        fd.write(f"{atom_type:>6} {mass:23.17g} # {s}\n")
    fd.write("\n")


def write_lammps_data(
    fileobj,
    atoms: Atoms,
    *,
    specorder: list = None,
    reduce_cell: bool = False,
    force_skew: bool = False,
    prismobj: Prism = None,
    write_image_flags: bool = True,
    masses: bool = False,
    velocities: bool = False,
    units: str = "metal",
    bonds: bool = True,
    angles: bool = True,
    atom_style: str = "full",
):
    """Write atomic structure data to a LAMMPS data file.

    Parameters
    ----------
    fileobj : file|str
        File object which the output will be written.
    atoms : Atoms
        Atoms to be written.
    specorder : list[str], optional
        Chemical symbols in the order of LAMMPS atom types, by default None
    force_skew : bool, optional
        Force to write the cell as a
        `triclinic <https://docs.lammps.org/Howto_triclinic.html>`__ box,
        by default False
    reduce_cell : bool, optional
        Whether the cell shape is reduced or not, by default False
    prismobj : Prism|None, optional
        Prism, by default None
    write_image_flags : bool, default True
        If True, the image flags, i.e., in which images of the periodic
        simulation box the atoms are, are written and positions are wrapped into the box.
    masses : bool, optional
        Whether the atomic masses are written or not, by default False
    velocities : bool, optional
        Whether the atomic velocities are written or not, by default False
    units : str, optional
        `LAMMPS units <https://docs.lammps.org/units.html>`__,
        by default 'metal'
    bonds : bool, optional
        Whether the bonds are written or not. Bonds can only be written
        for atom_style='full', by default True
    atom_style : {'atomic', 'charge', 'full'}, optional
        `LAMMPS atom style <https://docs.lammps.org/atom_style.html>`__,
        by default 'atomic'
    """

    if isinstance(atoms, list):
        if len(atoms) > 1:
            raise ValueError("Can only write one configuration to a lammps data file!")
        atoms = atoms[0]

    fileobj.write("(written by ASE)\n\n")

    symbols = atoms.get_chemical_symbols()
    n_atoms = len(symbols)
    fileobj.write(f"{n_atoms} atoms\n")

    if specorder is None:
        # This way it is assured that LAMMPS atom types are always
        # assigned predictably according to the alphabetic order
        species = sorted(set(symbols))
    else:
        # To index elements in the LAMMPS data file
        # (indices must correspond to order in the potential file)
        species = specorder
    n_atom_types = len(species)

    bonds_in = []
    if bonds and (atom_style == "full") and (atoms.arrays.get("bonds") is not None):
        n_bonds = 0
        n_bond_types = 1
        for i, bondsi in enumerate(atoms.arrays["bonds"]):
            if bondsi != "_":
                for bond in bondsi.split(","):
                    dummy1, dummy2 = bond.split("(")
                    bond_type = int(dummy2.split(")")[0])
                    at1 = int(i) + 1
                    at2 = int(dummy1) + 1
                    bonds_in.append((bond_type, at1, at2))
                    n_bonds = n_bonds + 1
                    if bond_type > n_bond_types:
                        n_bond_types = bond_type
        fileobj.write(f"{n_bonds} bonds\n")

    angles_in = []
    if angles and (atom_style == "full") and (atoms.arrays.get("angles") is not None):
        n_angles = 0
        n_angle_types = 1
        for i, anglesi in enumerate(atoms.arrays["angles"]):
            if anglesi != "_":
                for angle in anglesi.split(","):
                    dummy1, dummy2 = angle.split("(")
                    angle_type = int(dummy2.split(")")[0])
                    at2 = int(i) + 1
                    at1 = int(dummy1.split("-")[0]) + 1
                    at3 = int(dummy1.split("-")[1]) + 1
                    angles_in.append((angle_type, at1, at2, at3))
                    n_angles = n_angles + 1
                    if angle_type > n_angle_types:
                        n_angle_types = angle_type
        fileobj.write(f"{n_angles} angles\n")

    # types 
    fileobj.write(f"\n{n_atom_types} atom types\n")
    if bonds and (atom_style == "full") and (atoms.arrays.get("bonds") is not None):
        fileobj.write(f"{n_bond_types} bond types\n")
    if angles and (atom_style == "full") and (atoms.arrays.get("angles") is not None):
        fileobj.write(f"{n_angle_types} angle types\n")

    if prismobj is None:
        prismobj = Prism(atoms.get_cell(), reduce_cell=reduce_cell)

    # Get cell parameters and convert from ASE units to LAMMPS units
    xhi, yhi, zhi, xy, xz, yz = convert(
        prismobj.get_lammps_prism(), "distance", "ASE", units
    )

    fileobj.write(f"\n0.0 {xhi:23.17g}  xlo xhi\n")
    fileobj.write(f"0.0 {yhi:23.17g}  ylo yhi\n")
    fileobj.write(f"0.0 {zhi:23.17g}  zlo zhi\n")

    if force_skew or prismobj.is_skewed():
        fileobj.write(f"{xy:23.17g} {xz:23.17g} {yz:23.17g}  xy xz yz\n")
    fileobj.write("\n")

    if masses:
        _write_masses(fileobj, atoms, species, units)

    # Write (unwrapped) atomic positions.  If wrapping of atoms back into the
    # cell along periodic directions is desired, this should be done manually
    # on the Atoms object itself beforehand.
    fileobj.write(f"Atoms # {atom_style}\n\n")

    if write_image_flags:
        scaled_positions = atoms.get_scaled_positions(wrap=False)
        image_flags = np.floor(scaled_positions).astype(int)

    # when `write_image_flags` is True, the positions are wrapped while the
    # unwrapped positions can be recovered from the image flags
    pos = prismobj.vector_to_lammps(
        atoms.get_positions(),
        wrap=write_image_flags,
    )

    if atom_style == "atomic":
        # Convert position from ASE units to LAMMPS units
        pos = convert(pos, "distance", "ASE", units)
        for i, r in enumerate(pos):
            s = species.index(symbols[i]) + 1
            line = f"{i + 1:>6} {s:>3}" f" {r[0]:23.17g} {r[1]:23.17g} {r[2]:23.17g}"
            if write_image_flags:
                img = image_flags[i]
                line += f" {img[0]:6d} {img[1]:6d} {img[2]:6d}"
            line += "\n"
            fileobj.write(line)
    elif atom_style == "charge":
        charges = atoms.get_initial_charges()
        # Convert position and charge from ASE units to LAMMPS units
        pos = convert(pos, "distance", "ASE", units)
        charges = convert(charges, "charge", "ASE", units)
        for i, (q, r) in enumerate(zip(charges, pos)):
            s = species.index(symbols[i]) + 1
            line = (
                f"{i + 1:>6} {s:>3} {q:>5}"
                f" {r[0]:23.17g} {r[1]:23.17g} {r[2]:23.17g}"
            )
            if write_image_flags:
                img = image_flags[i]
                line += f" {img[0]:6d} {img[1]:6d} {img[2]:6d}"
            line += "\n"
            fileobj.write(line)
    elif atom_style == "full":
        charges = atoms.get_initial_charges()
        # The label 'mol-id' has apparenlty been introduced in read earlier,
        # but so far not implemented here. Wouldn't a 'underscored' label
        # be better, i.e. 'mol_id' or 'molecule_id'?
        if atoms.has("mol-id"):
            molecules = atoms.get_array("mol-id")
            if not np.issubdtype(molecules.dtype, np.integer):
                raise TypeError(
                    f'If "atoms" object has "mol-id" array, then '
                    f"mol-id dtype must be subtype of np.integer, and "
                    f"not {molecules.dtype!s:s}."
                )
            if (len(molecules) != len(atoms)) or (molecules.ndim != 1):
                raise TypeError(
                    'If "atoms" object has "mol-id" array, then '
                    "each atom must have exactly one mol-id."
                )
        else:
            # Assigning each atom to a distinct molecule id would seem
            # preferableabove assigning all atoms to a single molecule
            # id per default, as done within ase <= v 3.19.1. I.e.,
            # molecules = np.arange(start=1, stop=len(atoms)+1,
            # step=1, dtype=int) However, according to LAMMPS default
            # behavior,
            molecules = np.zeros(len(atoms), dtype=int)
            # which is what happens if one creates new atoms within LAMMPS
            # without explicitly taking care of the molecule id.
            # Quote from docs at https://lammps.sandia.gov/doc/read_data.html:
            #    The molecule ID is a 2nd identifier attached to an atom.
            #    Normally, it is a number from 1 to N, identifying which
            #    molecule the atom belongs to. It can be 0 if it is a
            #    non-bonded atom or if you don't care to keep track of molecule
            #    assignments.

        # Convert position and charge from ASE units to LAMMPS units
        pos = convert(pos, "distance", "ASE", units)
        charges = convert(charges, "charge", "ASE", units)
        for i, (m, q, r) in enumerate(zip(molecules, charges, pos)):
            s = species.index(symbols[i]) + 1
            line = (
                f"{i + 1:>6} {m:>3} {s:>3} {q:>5}"
                f" {r[0]:23.17g} {r[1]:23.17g} {r[2]:23.17g}"
            )
            if write_image_flags:
                img = image_flags[i]
                line += f" {img[0]:6d} {img[1]:6d} {img[2]:6d}"
            line += "\n"
            fileobj.write(line)
        if bonds and (atoms.arrays.get("bonds") is not None):
            fileobj.write("\nBonds\n\n")
            for i in range(n_bonds):
                bond_type = bonds_in[i][0]
                at1 = bonds_in[i][1]
                at2 = bonds_in[i][2]
                fileobj.write(f"{i + 1:>3} {bond_type:>3} {at1:>3} {at2:>3}\n")
        if angles and (atoms.arrays.get("angles") is not None):
            fileobj.write("\nAngles\n\n")
            for i in range(n_angles):
                angle_type = angles_in[i][0]
                at1 = angles_in[i][1]
                at2 = angles_in[i][2]
                at3 = angles_in[i][3]
                fileobj.write(
                    f"{i + 1:>3} {angle_type:>3} {at1:>3} {at2:>3} {at3:>3}\n"
                )
    else:
        raise ValueError(atom_style)

    if velocities and atoms.get_velocities() is not None:
        fileobj.write("\n\nVelocities\n\n")
        vel = prismobj.vector_to_lammps(atoms.get_velocities())
        # Convert velocity from ASE units to LAMMPS units
        vel = convert(vel, "velocity", "ASE", units)
        for i, v in enumerate(vel):
            fileobj.write(f"{i + 1:>6} {v[0]:23.17g} {v[1]:23.17g} {v[2]:23.17g}\n")

    fileobj.flush()


def write_gromos(fileobj, atoms: Atoms):
    """Write gromos geometry files (.g96), given an ASE atoms object.
    Arrays in the Atoms object that are used: atomtypes, residuenames and molecule_ids.

    If residuenames are not provided, all atoms get a residuename of 1DUM
    If atomtypes are not provided, chemical symbols from ASE are used
    If residuenumbers (molecule IDs) are not provided, every molecule ID is set as the index+1

    Writes out:
    atom positions,
    and simulation cell (if present)
    Positions are written out in nm (default in GROMACS) not Angstrom (default in ASE)
    """

    from ase import units

    natoms = len(atoms)
    try:
        gromos_residuenames = atoms.get_array("residuenames")
    except KeyError:
        gromos_residuenames = []
        for _ in range(natoms):
            gromos_residuenames.append("1DUM")
    try:
        gromos_atomtypes = atoms.get_array("atomtypes")
    except KeyError:
        gromos_atomtypes = atoms.get_chemical_symbols()
    try:
        gromos_molecule_ids = atoms.get_array("residuenumbers")
    except KeyError:
        gromos_molecule_ids = [i + 1 for i in range(natoms)]

    pos = atoms.get_positions()
    pos = pos / units.nm  # Convert units from Angstrom to nm

    vel = atoms.get_velocities()
    if vel is None:
        vel = pos * 0.0
    else:
        vel *= 1000.0 * units.fs / units.nm

    fileobj.write("TITLE\n")
    fileobj.write("Gromos96 structure file written by ASE \n")
    fileobj.write("END\n")
    fileobj.write("POSITION\n")
    count = 1
    rescount = 0
    old_molid = 0
    for resname, atomtype, mol_id, xyz in zip(
        gromos_residuenames, gromos_atomtypes, gromos_molecule_ids, pos
    ):
        if mol_id != old_molid:
            old_molid = mol_id
            rescount = rescount + 1
        okresname = resname.lstrip("0123456789 ")
        fileobj.write(
            "%5d %-5s %-5s%7d%15.9f%15.9f%15.9f\n"
            % (rescount, okresname, atomtype, count, xyz[0], xyz[1], xyz[2])
        )
        count = count + 1

    fileobj.write("END\n")

    if atoms.get_pbc().any():
        fileobj.write("BOX\n")
        mycell = atoms.get_cell()
        grocell = mycell.flat[[0, 4, 8, 1, 2, 3, 5, 6, 7]] * 0.1
        fileobj.write("".join([f"{x:15.9f}" for x in grocell]))
        fileobj.write("\nEND\n")
