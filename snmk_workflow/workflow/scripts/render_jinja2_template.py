from typing import Union
from snakemake.script import snakemake
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
import confi


def main(template: Path, output: Path, params: dict, input: dict):

    confi.render_jinja2(template, output, params, input)


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
