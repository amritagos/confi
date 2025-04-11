from typing import Union
from jinja2 import Environment, FileSystemLoader
from pathlib import Path


def render_jinja2(template: Path, output: Path, params: dict, input: dict):
    # Absolute path to where the template is
    template_path = Path(template).resolve()

    env = Environment(loader=FileSystemLoader(str(template_path.parent)))
    # Load the template using just its name (filename)
    template = env.get_template(template_path.name)

    with open(output, "w") as f:
        variables = dict(params=params, input=input)
        f.write(template.render(**variables))
