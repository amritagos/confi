from typing import Optional
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
from confi.parameters import PackmolInput, PackmolParams


def render_jinja2(templ: Path, variables: dict, output: Path):
    # Absolute path to where the template is
    template_path = Path(templ).resolve()

    env = Environment(loader=FileSystemLoader(str(template_path.parent)))
    # Load the template using just its name (filename)
    template = env.get_template(template_path.name)

    with open(output, "w") as f:
        f.write(template.render(**variables))


def render_packmol_input(
    output: Path,
    params: PackmolParams,
    input: PackmolInput,
    template: Optional[Path] = None,
    **kwargs,
):
    """Renders a PACKMOL input file that can then be run with PACKMOL.

    Args:
        output (Path): The output path where the PACKMOL input file will be written out
        params (PackmolParams): Parameters that are used in the PACKMOL input jinja2 template.
        input (PackmolInput): Input parameters used in the PACKMOL input jinja2 template
        template (Optional[Path], optional): The Jinja2 template used for PACKMOL. Defaults to None (essentially the template file for PACKMOL inside confi/templates).
        kwargs: These are extra keywords that can be used in a user-defined jinja2 template. Inside the template, use like so: extra.kwarg
    """
    if template is None:
        src_dir = Path(__file__).parent.resolve()
        template = src_dir / Path(f"templates/packmol_template.inp")

    variables = dict(params=params.model_dump(), input=input.model_dump(), extra=kwargs)

    render_jinja2(template, variables, output)
