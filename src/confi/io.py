from ase import Atoms, units


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
