import subprocess
from snakemake.script import snakemake
from pathlib import Path
import confi
from confi.parameters import FileType
from ase.io import read


def main(
    snmk_params: dict, snmk_input: dict, packmol_pdb: Path, g96_file: Path
) -> None:

    # Create the folders required for PACKMOL and ASE
    packmol_pdb = packmol_pdb.resolve()
    packmol_pdb.parent.mkdir(parents=True, exist_ok=True)
    g96_file.parent.mkdir(parents=True, exist_ok=True)
    # Intermediate output files for PACKMOL and Gromacs
    output_packmol_inp = packmol_pdb.parent / f"packmol.inp"

    # Use absolute paths wherever possible
    packmol_input = confi.parameters.PackmolInput(
        filetype=FileType.PDB,
        cation_file=snmk_input.get("packmol_cation_file"),
        anion_file=snmk_input.get("packmol_anion_file"),
        water_file=snmk_input.get("packmol_water_file"),
    )
    packmol_params = confi.parameters.PackmolParams(
        packmol_input=packmol_input,
        monomer_file=snmk_params.get("packmol_monomer_file"),
        dimer_file=snmk_params.get("packmol_dimer_file"),
        system_file=packmol_pdb,  # the actual PDB that will be created by PACKMOL if you run it with the input file you will create
        resnumbers=3,
        **snmk_params
    )

    # Actually render the PACKMOL input file :
    confi.render.render_packmol_input(output_packmol_inp, packmol_params, packmol_input)

    # Run PACKMOL
    command = f"packmol < {output_packmol_inp}"
    # Run the command using subprocess
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Packmol command failed with error: {e}")

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

    # Write out the gromacs files
    with open(g96_file, "w") as f:
        confi.io.write_gromos(f, atoms)


if __name__ == "__main__":

    packmol_pdb = Path(
        snakemake.output.get("packmol")
    )  # PDB file that will be written out by PACKMOL
    g96_file = Path(
        snakemake.output.get("g96_file")
    )  # g96 file that will be created by moltemplate
    params = snakemake.params
    input = snakemake.input

    main(
        snmk_params=params,
        snmk_input=input,
        packmol_pdb=packmol_pdb,
        g96_file=g96_file,
    )
