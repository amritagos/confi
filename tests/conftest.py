import pytest
import confi
from pathlib import Path
import subprocess
from ase.io import read

from confi.parameters import FileType


def check_xyz(
    xyz_file: Path,
    cation_anion_bond_length: float = 2.33,
    cation_symbol: str = "Fe",
    anion_symbol: str = "Cl",
    oxygen_symbol: str = "O",
    tolerance: float = 0.001,
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
                    <= cation_anion_bond_length + tolerance
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


@pytest.fixture
def test_packmol_input_lammps():
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


@pytest.fixture
def test_packmol_input_gromacs_xyz():
    """Create the input file for a PACKMOL input file and run PACKMOL to create an XYZ file for a system with
    5 water molecules, 2 free cations (Fe), 2 monomers, 2 dimer (and adjust the anions, Cl)
    """
    package_root = Path(__file__).parent.parent.resolve()  # top-level directory
    output_packmol_inp = package_root / Path(
        "tests/output/packmol_xyz.inp"
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
        / Path(f"resources/packmol/gromacs/tip4p_2005/tip4p_2005_water.xyz"),
    )
    params = confi.parameters.PackmolParams(
        packmol_input=input,
        monomer_file=package_root / Path(f"resources/packmol/fe_cl.xyz"),
        dimer_file=package_root / Path("resources/packmol/fe_cl2.xyz"),
        system_file=package_root
        / Path(
            "tests/output/system_gromacs.xyz"
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


@pytest.fixture
def test_packmol_input_gromacs_pdb():
    """Create the input file for a PACKMOL input file and run PACKMOL to create a PDB file for a system with
    5 water molecules, 2 free cations (Fe), 2 monomers, 2 dimer (and adjust the anions, Cl)
    """
    package_root = Path(__file__).parent.parent.resolve()  # top-level directory
    output_packmol_inp = package_root / Path(
        "tests/output/packmol_pdb.inp"
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
        filetype=FileType.PDB,
        cation_file=package_root / Path(f"resources/packmol/fe.pdb"),
        anion_file=package_root / Path(f"resources/packmol/cl.pdb"),
        water_file=package_root
        / Path(f"resources/packmol/gromacs/tip4p_2005/tip4p_2005_water.pdb"),
    )
    params = confi.parameters.PackmolParams(
        packmol_input=input,
        monomer_file=package_root / Path(f"resources/packmol/fe_cl.pdb"),
        dimer_file=package_root / Path("resources/packmol/fe_cl2.pdb"),
        system_file=package_root
        / Path(
            "tests/output/system_gromacs.pdb"
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

    return params
