from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Optional

import typer
from ase.io import read

from confi import io, render
from confi.parameters import (
    FileType,
    MoltemplateInput,
    MoltemplateParams,
    PackmolInput,
    PackmolParams,
)

app = typer.Typer(
    help="CLI helpers for generating Confi inputs and outputs without Snakemake."
)


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _require_file(path: Optional[Path], label: str) -> None:
    if path is None:
        return
    if not path.exists():
        raise typer.BadParameter(f"{label} does not exist: {path}")
    if not path.is_file():
        raise typer.BadParameter(f"{label} is not a file: {path}")


def _require_command(name: str) -> None:
    if shutil.which(name) is None:
        raise typer.BadParameter(
            f"Required external command '{name}' was not found in PATH."
        )


def _run_shell(command: str, cwd: Optional[Path] = None) -> None:
    try:
        subprocess.run(
            command,
            shell=True,
            check=True,
            cwd=str(cwd) if cwd else None,
        )
    except subprocess.CalledProcessError as exc:
        raise typer.Exit(code=exc.returncode) from exc


def _validate_common_packmol_inputs(
    filetype: FileType,
    cation_file: Optional[Path],
    anion_file: Optional[Path],
    water_file: Optional[Path],
    monomer_file: Optional[Path],
    dimer_file: Optional[Path],
    trimer_file: Optional[Path],
    tetrahedral_file: Optional[Path],
) -> PackmolInput:
    for label, path in {
        "cation file": cation_file,
        "anion file": anion_file,
        "water file": water_file,
        "monomer file": monomer_file,
        "dimer file": dimer_file,
        "trimer file": trimer_file,
        "tetrahedral file": tetrahedral_file,
    }.items():
        _require_file(path, label)

    return PackmolInput(
        filetype=filetype,
        cation_file=cation_file,
        anion_file=anion_file,
        water_file=water_file,
    )


@app.command("build-lammps")
def build_lammps(
    output_dir: Path = typer.Option(
        Path("results/lammps"),
        help="Directory for all generated files.",
    ),
    cation_file: Optional[Path] = typer.Option(
        None,
        help="PACKMOL cation structure file.",
    ),
    anion_file: Optional[Path] = typer.Option(
        None,
        help="PACKMOL anion structure file.",
    ),
    water_file: Optional[Path] = typer.Option(
        None,
        help="PACKMOL water structure file.",
    ),
    monomer_file: Optional[Path] = typer.Option(
        None,
        help="PACKMOL monomer structure file.",
    ),
    dimer_file: Optional[Path] = typer.Option(
        None,
        help="PACKMOL dimer structure file.",
    ),
    moltemplate_cation_file: Optional[Path] = typer.Option(
        None,
        help="Moltemplate cation LT file.",
    ),
    moltemplate_anion_file: Optional[Path] = typer.Option(
        None,
        help="Moltemplate anion LT file.",
    ),
    moltemplate_water_file: Optional[Path] = typer.Option(
        None,
        help="Moltemplate water LT file.",
    ),
    x_box_length: float = typer.Option(..., min=0.0),
    y_box_length: float = typer.Option(..., min=0.0),
    z_box_length: float = typer.Option(..., min=0.0),
    n_wat: int = typer.Option(0, min=0),
    n_free_cations: int = typer.Option(0, min=0),
    n_free_anions: int = typer.Option(0, min=0),
    n_monomer: int = typer.Option(0, min=0),
    n_dimer: int = typer.Option(0, min=0),
    tolerance: float = typer.Option(2.0, min=0.0),
    cation_radius: float = typer.Option(1.5, min=0.0),
    anion_radius: float = typer.Option(1.5, min=0.0),
    monomer_radius: float = typer.Option(1.5, min=0.0),
    dimer_radius: float = typer.Option(1.5, min=0.0),
    seed: int = typer.Option(-1),
    keep_intermediate: bool = typer.Option(
        False,
        help="Keep packmol.inp and system.lt.",
    ),
) -> None:
    """Build a LAMMPS data file by running PACKMOL and moltemplate."""
    _require_command("packmol")
    _require_command("moltemplate.sh")

    packmol_input = _validate_common_packmol_inputs(
        FileType.XYZ,
        cation_file,
        anion_file,
        water_file,
        monomer_file,
        dimer_file,
        None,
        None,
    )

    for label, path in {
        "moltemplate cation file": moltemplate_cation_file,
        "moltemplate anion file": moltemplate_anion_file,
        "moltemplate water file": moltemplate_water_file,
    }.items():
        _require_file(path, label)

    output_dir = output_dir.resolve()
    packmol_xyz = output_dir / "system.xyz"
    packmol_inp = output_dir / "packmol.inp"
    moltemplate_lt = output_dir / "system.lt"
    lammps_data = output_dir / "system.data"
    _ensure_parent(packmol_xyz)

    packmol_params = PackmolParams(
        packmol_input=packmol_input,
        system_file=packmol_xyz,
        x_box_length=x_box_length,
        y_box_length=y_box_length,
        z_box_length=z_box_length,
        n_wat=n_wat,
        n_free_cations=n_free_cations,
        n_free_anions=n_free_anions,
        n_monomer=n_monomer,
        n_dimer=n_dimer,
        water_file=water_file,
        cation_file=cation_file,
        anion_file=anion_file,
        monomer_file=monomer_file,
        dimer_file=dimer_file,
        tolerance=tolerance,
        cation_radius=cation_radius,
        anion_radius=anion_radius,
        monomer_radius=monomer_radius,
        dimer_radius=dimer_radius,
        seed=seed,
    )
    render.render_packmol_input(packmol_inp, packmol_params, packmol_input)
    _run_shell(f"packmol < {packmol_inp}")

    moltemplate_input = MoltemplateInput(
        cation_file=moltemplate_cation_file,
        anion_file=moltemplate_anion_file,
        water_file=moltemplate_water_file,
    )
    moltemplate_params = MoltemplateParams(
        moltemplate_input=moltemplate_input,
        n_free_cations=n_free_cations,
        n_free_anions=n_free_anions,
        n_wat=n_wat,
        n_monomer=n_monomer,
        n_dimer=n_dimer,
        x_box_length=x_box_length,
        y_box_length=y_box_length,
        z_box_length=z_box_length,
    )
    render.render_moltemplate_input(
        moltemplate_lt,
        params=moltemplate_params,
        input=moltemplate_input,
    )

    _run_shell(
        f'moltemplate.sh -atomstyle "full" {moltemplate_lt.name} -xyz {packmol_xyz}',
        cwd=output_dir,
    )

    if not keep_intermediate:
        for path in [
            output_dir / "output_ttree",
            output_dir / "system.in.init",
            output_dir / "system.in.settings",
        ]:
            if path.is_dir():
                shutil.rmtree(path, ignore_errors=True)
            elif path.exists():
                path.unlink()

        for extra in output_dir.glob("run*"):
            if extra.is_dir():
                shutil.rmtree(extra, ignore_errors=True)
            else:
                extra.unlink()

    typer.echo(f"Created {lammps_data}")


@app.command("build-gromacs")
def build_gromacs(
    output_dir: Path = typer.Option(
        Path("results/gromacs"),
        help="Directory for all generated files.",
    ),
    cation_file: Optional[Path] = typer.Option(
        None,
        help="PACKMOL cation structure file.",
    ),
    anion_file: Optional[Path] = typer.Option(
        None,
        help="PACKMOL anion structure file.",
    ),
    water_file: Optional[Path] = typer.Option(
        None,
        help="PACKMOL water structure file.",
    ),
    monomer_file: Optional[Path] = typer.Option(
        None,
        help="PACKMOL monomer structure file.",
    ),
    dimer_file: Optional[Path] = typer.Option(
        None,
        help="PACKMOL dimer structure file.",
    ),
    trimer_file: Optional[Path] = typer.Option(
        None,
        help="PACKMOL trimer structure file.",
    ),
    tetrahedral_file: Optional[Path] = typer.Option(
        None,
        help="PACKMOL tetrahedral-shell structure file.",
    ),
    x_box_length: float = typer.Option(..., min=0.0),
    y_box_length: float = typer.Option(..., min=0.0),
    z_box_length: float = typer.Option(..., min=0.0),
    n_wat: int = typer.Option(0, min=0),
    n_free_cations: int = typer.Option(0, min=0),
    n_free_anions: int = typer.Option(0, min=0),
    n_monomer: int = typer.Option(0, min=0),
    n_dimer: int = typer.Option(0, min=0),
    n_trimer: int = typer.Option(0, min=0),
    n_tetrahedral: int = typer.Option(0, min=0),
    tolerance: float = typer.Option(2.0, min=0.0),
    cation_radius: float = typer.Option(1.5, min=0.0),
    anion_radius: float = typer.Option(1.5, min=0.0),
    monomer_radius: float = typer.Option(1.5, min=0.0),
    dimer_radius: float = typer.Option(1.5, min=0.0),
    trimer_radius: float = typer.Option(1.5, min=0.0),
    tetrahedral_radius: float = typer.Option(1.5, min=0.0),
    seed: int = typer.Option(-1),
    resnumbers: int = typer.Option(
        3,
        min=0,
        help="PACKMOL resnumbers setting for PDB output.",
    ),
    keep_intermediate: bool = typer.Option(
        False,
        help="Keep packmol.inp and system.pdb.",
    ),
) -> None:
    """Build a GROMACS .g96 file by running PACKMOL and ASE."""
    _require_command("packmol")

    packmol_input = _validate_common_packmol_inputs(
        FileType.PDB,
        cation_file,
        anion_file,
        water_file,
        monomer_file,
        dimer_file,
        trimer_file,
        tetrahedral_file,
    )

    output_dir = output_dir.resolve()
    packmol_pdb = output_dir / "system.pdb"
    packmol_inp = output_dir / "packmol.inp"
    g96_file = output_dir / "system.g96"
    _ensure_parent(packmol_pdb)

    packmol_params = PackmolParams(
        packmol_input=packmol_input,
        system_file=packmol_pdb,
        x_box_length=x_box_length,
        y_box_length=y_box_length,
        z_box_length=z_box_length,
        n_wat=n_wat,
        n_free_cations=n_free_cations,
        n_free_anions=n_free_anions,
        n_monomer=n_monomer,
        n_dimer=n_dimer,
        n_trimer=n_trimer,
        n_tetrahedral=n_tetrahedral,
        water_file=water_file,
        cation_file=cation_file,
        anion_file=anion_file,
        monomer_file=monomer_file,
        dimer_file=dimer_file,
        trimer_file=trimer_file,
        tetrahedral_file=tetrahedral_file,
        tolerance=tolerance,
        cation_radius=cation_radius,
        anion_radius=anion_radius,
        monomer_radius=monomer_radius,
        dimer_radius=dimer_radius,
        trimer_radius=trimer_radius,
        tetrahedral_radius=tetrahedral_radius,
        seed=seed,
        resnumbers=resnumbers,
    )
    render.render_packmol_input(packmol_inp, packmol_params, packmol_input)
    _run_shell(f"packmol < {packmol_inp}")

    atoms = read(packmol_pdb, format="proteindatabank")
    atoms.set_cell([x_box_length, y_box_length, z_box_length])
    atoms.set_pbc([True, True, True])
    with open(g96_file, "w") as handle:
        io.write_gromos(handle, atoms)

    if not keep_intermediate:
        for path in [packmol_inp, packmol_pdb]:
            if path.exists():
                path.unlink()

    typer.echo(f"Created {g96_file}")


if __name__ == "__main__":
    app()