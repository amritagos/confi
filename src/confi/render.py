from typing import Optional, Union
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
from pydantic import BaseModel


class PackmolInput(BaseModel):
    cation_file: Path
    anion_file: Path
    water_file: Path


class PackmolParams(BaseModel):
    system_file: Path  # output path to which PACKMOL will create the file if you actually run it with the created input file
    x_box_length: float  # Box length in the x-dimension
    y_box_length: float  # Box length in the y-dimension
    z_box_length: float  # Box length in the z-dimension
    n_wat: Optional[int] = None  # Number of water molecules
    water_file: Optional[Path] = (
        None  # if n_wat is greater than 0 then this cannot be None
    )
    n_free_cations: Optional[int] = None
    cation_file: Optional[Path] = (
        None  # if n_free_cations is greater than 0 then this cannot be None
    )
    cation_radius: Optional[float] = 1.75  # overrides the tolerance
    n_free_anions: Optional[int] = None  # free anions
    anion_file: Optional[Path] = (
        None  # if n_free_cations is greater than 0 then this cannot be None
    )
    anion_radius: Optional[float] = 1.75  # overrides the tolerance
    n_monomer: Optional[int] = None  # Number of cation-anion units
    monomer_file: Optional[Path] = (
        None  # if n_monomer is greater than 0 you must input this or PACKMOL will fail
    )
    n_dimer: Optional[int] = None  # Number of anion-cation-anion units
    dimer_file: Optional[Path] = (
        None  # if n_dimer is greater than 0 you must input this or PACKMOL will fail
    )
    tolerance: float = (
        2.0  # minimum distance between subunits in PACKMOL. can be overriden by radius
    )
    seed: int = (
        -1
    )  # PACKMOL uses -1 to generate a seed automatically from the computer time.


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
