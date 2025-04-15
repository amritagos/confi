from snakemake.script import snakemake
from pathlib import Path
import confi


def main(output: Path, snmk_params: dict, snmk_input: dict):

    input = confi.parameters.MoltemplateInput(**snmk_input)
    params = confi.parameters.MoltemplateParams(moltemplate_input=input, **snmk_params)
    # Use the default Moltemplate template
    confi.render.render_moltemplate_input(output, params, input)


if __name__ == "__main__":

    output = Path(snakemake.output[0])
    params = snakemake.params
    input = snakemake.input

    main(
        output=output,
        snmk_params=params,
        snmk_input=input,
    )
