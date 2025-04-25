import subprocess
from snakemake.script import snakemake
from pathlib import Path
import confi

def get_prefix(word:str):
    tokens = word.split("_")
    prefix = word.split("_")[0]
    stripped_word = ""
    for i in range(1, len(tokens)):
        stripped_word += tokens[i]
    return word.split("_")[0]

def create_pydantic_model(model_cls, params : dict, prefix : str):
    prog_specific_prefixes = ["packmol", "moltemplate"]
    params = dict()
    for k,v in params.items():
        key_prefix = get_prefix(k)
        if key_prefix in prog_specific_prefixes:
            # Add the key if the key_prefix is the desired prefix
            if key_prefix == prefix:
                model_cls
            

def main(
    snmk_params: dict, snmk_input: dict, packmol_xyz: Path, moltemplate_data: Path
) -> None:

    # Create the folders required for PACKMOL and moltemplate
    packmol_xyz = packmol_xyz.resolve()
    packmol_xyz.parent.mkdir(parents=True, exist_ok=True)
    moltemplate_data.parent.mkdir(parents=True, exist_ok=True)
    # Intermediate output files for PACKMOL and Moltemplate
    output_packmol_inp = packmol_xyz.parent / f"packmol.inp"
    output_moltemplate_inp = moltemplate_data.parent / f"system.lt"

    # Create the classes needed to create the packmol input
    packmol_input = confi.parameters.PackmolInput(
        cation_file=snmk_input.get("packmol_cation_file"),
        anion_file=snmk_input.get("packmol_anion_file"),
        water_file=snmk_input.get("packmol_water_file"),
    )
    packmol_params = confi.parameters.PackmolParams(
        packmol_input=packmol_input,
        monomer_file=snmk_params.get("packmol_monomer_file"),
        dimer_file=snmk_params.get("packmol_dimer_file"),
        system_file=packmol_xyz,  # the actual XYZ that will be created by PACKMOL if you run it with the input file you will create
        n_wat=snmk_params.get("n_wat"),
        n_monomer=snmk_params.get("n_monomer"),
        n_dimer=snmk_params.get("n_dimer"),
        n_free_cations=snmk_params.get("n_free_cations"),
        n_free_anions=snmk_params.get("n_free_anions"),
        x_box_length=snmk_params.get("x_box_length"),
        y_box_length=snmk_params.get("y_box_length"),
        z_box_length=snmk_params.get("z_box_length"),
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

    moltemplate_input = confi.parameters.MoltemplateInput(
        cation_file=snmk_input.get("moltemplate_cation_file"),
        anion_file=snmk_input.get("moltemplate_anion_file"),
        water_file=snmk_input.get("moltemplate_water_file"),
    )
    moltemplate_params = confi.parameters.MoltemplateParams(
        moltemplate_input=moltemplate_input,
        n_wat=snmk_params.get("n_wat"),
        n_monomer=snmk_params.get("n_monomer"),
        n_dimer=snmk_params.get("n_dimer"),
        n_free_cations=snmk_params.get("n_free_cations"),
        n_free_anions=snmk_params.get("n_free_anions"),
        x_box_length=snmk_params.get("x_box_length"),
        y_box_length=snmk_params.get("y_box_length"),
        z_box_length=snmk_params.get("z_box_length"),
    )

    # Render the moltemplate input file using the default moltemplate template:
    confi.render.render_moltemplate_input(
        output_moltemplate_inp, params=moltemplate_params, input=moltemplate_input
    )

    # Run Moltemplate and profit
    command = f'moltemplate.sh -atomstyle "full" {output_moltemplate_inp.name} -xyz {packmol_xyz} && rm -rf output_ttree system.in.init system.in.settings run*'

    # Run the moltemplate command using subprocess
    try:
        subprocess.run(
            command, shell=True, check=True, cwd=str(output_moltemplate_inp.parent)
        )
    except subprocess.CalledProcessError as e:
        raise Exception(f"Moltemplate command failed with error: {e}")


if __name__ == "__main__":

    packmol_xyz = Path(
        snakemake.output.get("packmol")
    )  # XYZ file that will be written out by PACKMOL
    moltemplate_data = Path(
        snakemake.output.get("moltemplate")
    )  # LAMMPS data file that will be created by moltemplate
    params = snakemake.params
    input = snakemake.input

    main(
        snmk_params=params,
        snmk_input=input,
        packmol_xyz=packmol_xyz,
        moltemplate_data=moltemplate_data,
    )
