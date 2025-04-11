from typing import Union
from snakemake.script import snakemake
from jinja2 import Environment, FileSystemLoader
from pathlib import Path


def main(template: Path, output: Path, params: dict, input: dict):
    # Absolute path to where the template is
    template_path = Path(template).resolve()

    env = Environment(loader=FileSystemLoader(str(template_path.parent)))
    # Load the template using just its name (filename)
    template = env.get_template(template_path.name)

    with open(output, "w") as f:
        variables = dict(params=params, input=input)
        f.write(template.render(**variables))


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
