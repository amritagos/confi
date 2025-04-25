from collections import Counter
import re
import pytest
import confi
from pathlib import Path
import subprocess
from ase.io import read, lammpsdata

import confi.parameters


def check_xyz(
    xyz_file: Path,
    cation_anion_bond_length: float = 2.34,
    cation_symbol: str = "Fe",
    anion_symbol: str = "Cl",
    oxygen_symbol: str = "O",
):
    atoms = read(xyz_file)
    total_anions = 0
    total_cations = 0
    free_cations = 0
    n_monomers = 0
    n_dimers = 0
    n_waters = 0

    for atom in atoms:
        if atom.symbol == oxygen_symbol:
            n_waters += 1
        elif atom.symbol == cation_symbol:
            total_cations += 1
            # Search for neighbouring anions
            n_anion_neigh = 0
            for j in range(len(atoms)):
                if atoms[j].symbol != anion_symbol:
                    continue
                if (
                    atoms.get_distance(atom.index, j, mic=True)
                    <= cation_anion_bond_length
                ):
                    n_anion_neigh += 1
            if n_anion_neigh == 1:
                n_monomers += 1
            elif n_anion_neigh == 2:
                n_dimers += 1
            else:
                free_cations += 1
        elif atom.symbol == anion_symbol:
            total_anions += 1

    # Calculate the number of free anions, using total_anions and the number of monomers and dimers
    free_anions = total_anions - n_monomers - 2 * n_dimers

    return (
        n_waters,
        total_cations,
        total_anions,
        free_cations,
        free_anions,
        n_monomers,
        n_dimers,
    )


def count_bond_types(bonds_array):
    """
    Count the number of bonds for each bond type based on a bonds array.

    Each element in bonds_array is either an underscore "_" (meaning no bond) or
    a string of comma‑separated bond specifications like "1(2)" where the number
    in parentheses is the bond type.

    Parameters
    ----------
    bonds_array : array-like of str
        The bonds information from the ASE Atoms object.

    Returns
    -------
    dict
        A dictionary mapping bond type (int) to the number of occurrences.
    """
    bond_counter = Counter()
    for entry in bonds_array:
        # Skip blank or underscore entries (underscore denotes no bond)
        if not entry or entry.strip() == "_":
            continue
        # Split multiple bonds separated by commas
        bonds = entry.split(",")
        for bond in bonds:
            # Use regex to extract the bond type from inside parentheses
            match = re.search(r"\((\d+)\)", bond)
            if match:
                bond_type = int(match.group(1))
                bond_counter[bond_type] += 1
    return dict(bond_counter)


@pytest.fixture
def test_packmol_input():
    """Create the input file for a PACKMOL input file and run PACKMOL to create an XYZ file for a system with
    5 water molecules, 2 free cations (Fe), 2 monomers, 2 dimer (and adjust the anions, Cl)
    """
    package_root = Path(__file__).parent.parent.resolve()  # top-level directory
    output_packmol_inp = package_root / Path(
        "tests/output/packmol.inp"
    )  # relative to the top-level directory

    # Create the output directory
    output_packmol_inp.parent.mkdir(parents=True, exist_ok=True)

    # We will create a system with 5 water molecules, 2 free cations, 2 monomers (Fe-Cl) and 1 dimer (Cl-Fe-Cl)
    # We are simulating FeCl2, so the total number of anions is double the number of cations
    n_wat = 5  # H2O
    free_cations = 2  # free Fe
    n_monomers = 2  # Fe-Cl
    n_dimers = 1  # Cl-Fe-Cl
    total_cations = free_cations + n_monomers + n_dimers
    total_anions = 2 * total_cations  # FeCl2 being simulated
    free_anions = total_anions - n_monomers - 2 * n_dimers
    box_length = 10  # in Angstrom, as are all distances

    # Create PackmolParams and PackmolInput for the parameters and inputs required by the Jinja2 renderer
    # Use absolute paths wherever possible
    input = confi.parameters.PackmolInput(
        cation_file=package_root / Path(f"resources/packmol/fe.xyz"),
        anion_file=package_root / Path(f"resources/packmol/cl.xyz"),
        water_file=package_root
        / Path(f"resources/packmol/lammps/tip4p_2005/tip4p_2005_water.xyz"),
    )
    params = confi.parameters.PackmolParams(
        packmol_input=input,
        monomer_file=package_root / Path(f"resources/packmol/fe_cl.xyz"),
        dimer_file=package_root / Path("resources/packmol/fe_cl2.xyz"),
        system_file=package_root
        / Path(
            "tests/output/system.xyz"
        ),  # the actual XYZ that will be created by PACKMOL if you run it with the input file you will create
        n_wat=n_wat,
        n_monomer=n_monomers,
        n_dimer=n_dimers,
        n_free_cations=free_cations,
        n_free_anions=free_anions,
        x_box_length=box_length,
        y_box_length=box_length,
        z_box_length=box_length,
    )

    # Actually render the input file :
    confi.render.render_packmol_input(output_packmol_inp, params, input)

    # Check that the packmol input file exists
    assert output_packmol_inp.exists(), "Packmol input file was not created."

    # 'packmol' should be in your PATH if you're using the micromamba environment
    command = f"packmol < {output_packmol_inp}"

    # Run the command using subprocess
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        pytest.fail(f"Packmol command failed with error: {e}")

    # Check the output file created
    (
        n_wat_read,
        total_cations_read,
        total_anions_read,
        free_cations_read,
        free_anions_read,
        n_monomers_read,
        n_dimers_read,
    ) = check_xyz(xyz_file=params.system_file)

    # The total number of cations etc should be the same as we put in
    assert n_wat_read == n_wat
    assert total_cations_read == total_cations
    assert total_anions_read == total_anions
    assert free_cations_read == free_cations
    assert free_anions_read == free_anions
    assert n_monomers_read == n_monomers
    assert n_dimers_read == n_dimers

    return params


def test_moltemplate_input(test_packmol_input):
    """Using the PACKMOL generated XYZ file as input, create a LAMMPS data file"""
    package_root = Path(__file__).parent.parent.resolve()  # top-level directory
    output_moltemplate_inp = package_root / Path(
        "tests/output/moltemplate/system.lt"
    )  # relative to the top-level directory
    output_moltemplate_inp.parent.mkdir(parents=True, exist_ok=True)

    packmol_params = test_packmol_input

    # Dictionary from the PackmolParams class
    packmol_param_dict = packmol_params.model_dump()

    # Create MoltemplateParams and MoltemplateInput for the parameters and inputs required by the Jinja2 renderer
    # Use absolute paths wherever possible
    input = confi.parameters.MoltemplateInput(
        cation_file=package_root / Path(f"resources/moltemplate/fe_ions.lt"),
        anion_file=package_root / Path(f"resources/moltemplate/cl.lt"),
        water_file=package_root / Path(f"resources/moltemplate/tip4p_2005.lt"),
    )
    params = confi.parameters.MoltemplateParams(
        moltemplate_input=input,
        n_free_cations=packmol_param_dict["n_free_cations"],
        n_free_anions=packmol_param_dict["n_free_anions"],
        n_wat=packmol_param_dict["n_wat"],
        n_monomer=packmol_param_dict["n_monomer"],
        n_dimer=packmol_param_dict["n_dimer"],
        x_box_length=packmol_param_dict["x_box_length"],
        y_box_length=packmol_param_dict["y_box_length"],
        z_box_length=packmol_param_dict["z_box_length"],
    )

    # Render the moltemplate input file using the default moltemplate template:
    confi.render.render_moltemplate_input(
        output_moltemplate_inp, params=params, input=input
    )

    command = f'moltemplate.sh -atomstyle "full" {output_moltemplate_inp} -xyz {packmol_param_dict["system_file"]}'

    # Run the command using subprocess
    try:
        subprocess.run(
            command, shell=True, check=True, cwd=str(output_moltemplate_inp.parent)
        )
    except subprocess.CalledProcessError as e:
        pytest.fail(f"Moltemplate command failed with error: {e}")

    # Read the lammps data file
    Z_of_type = {1: 8, 2: 1, 3: 26, 4: 17}  # Mapping from LAMMPS type to atom number
    data_file_path = output_moltemplate_inp.parent / f"system.data"
    with open(data_file_path) as f:
        atoms = lammpsdata.read_lammps_data(
            f, Z_of_type=Z_of_type, units="metal", atom_style="full"
        )

    # Check the total number of atoms
    total_atom_exp = (
        packmol_param_dict["n_free_cations"]
        + packmol_param_dict["n_free_anions"]
        + 3 * packmol_param_dict["n_wat"]
        + 2 * packmol_param_dict["n_monomer"]
        + 3 * packmol_param_dict["n_dimer"]
    )
    assert (
        len(atoms) == total_atom_exp
    ), f"Expected {total_atom_exp} atoms, but got {len(atoms)}"

    # Test bonds: ASE's read_lammps_data stores bonds as a numpy array of strings
    # each string corresponds to the atom bonded to, with the type in brackets
    bonds = atoms.arrays.get("bonds", None)
    assert bonds is not None  # Bonds should have been written out
    # Count non-empty entries (where '_' means no bond) and bond types
    bond_info = count_bond_types(bonds)
    n_bonds = sum(bond_info.values())
    num_bond_types = len(bond_info)
    # Check the total number of bonds
    n_bonds_exp = (
        2 * packmol_param_dict["n_wat"]
        + packmol_param_dict["n_monomer"]
        + 2 * packmol_param_dict["n_dimer"]
    )  # expected number of bonds
    assert n_bonds == n_bonds_exp, f"Expected {n_bonds_exp} bonds, got {n_bonds}"
    # Check the number of bond types
    assert num_bond_types == 2, f"Expected 2 bond types, got {num_bond_types}"

    # There should be angles for all the water molecules in the system

    angles = atoms.arrays.get("angles", None)
    assert angles is not None  # There should be one angle for each water molecule
    n_angles = sum(1 for a in angles if a != "_")
    n_angles_exp = packmol_param_dict["n_wat"]
    assert n_angles == n_angles_exp, f"Expected {n_angles_exp} angles, got {n_angles}"
