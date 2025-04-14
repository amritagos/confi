from snakemake.script import snakemake
from pathlib import Path
import confi


def main(template: Path, output: Path, params: dict, input: dict):

    variables = dict(params=params, input=input)
    confi.render_jinja2(template, variables, output)


if __name__ == "__main__":

    template = Path(snakemake.input.template)
    output = Path(snakemake.output[0])
    params = snakemake.params
    input = snakemake.input

    main(
        template=template,
        output=output,
        params=params,
        input=input,
    )
