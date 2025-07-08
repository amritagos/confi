from snakemake.script import snakemake
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
import shutil

import sys


def render_jinja2(template_folder: Path, variables: dict, output: Path):

    template_file = variables["input"]["template"]
    potential_file = variables["input"]["potential_file"]
    shutil.copy(template_file, template_folder)
    shutil.copy(potential_file, template_folder)

    variables["input"]["potential_file"] = Path(
        variables["input"]["potential_file"]
    ).name

    env = Environment(loader=FileSystemLoader(str(template_folder)))
    # Load the template using just its name (filename)
    template = env.get_template(Path(template_file).name)

    with open(output, "w") as f:
        f.write(template.render(**variables))


variables = dict(params=dict(snakemake.params), input=dict(snakemake.input))

output = snakemake.output[0]  # the actual rendered file
template_folder = snakemake.params[
    "template_folder"
]  # path to where the template file is. I guess this could be relative?

# Make the folder for the templates
Path(template_folder).mkdir(parents=True, exist_ok=True)

render_jinja2(template_folder, variables, output)
