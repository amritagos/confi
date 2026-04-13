import shutil
from pathlib import Path

import pytest
from typer.testing import CliRunner

from confi.cli import app

runner = CliRunner()


def _has_command(name: str) -> bool:
    return shutil.which(name) is not None


@pytest.mark.skipif(
    not _has_command("packmol"),
    reason="PACKMOL is required for CLI test",
)
def test_build_gromacs_cli_real(tmp_path):
    """
    This test actually runs PACKMOL.
    It uses existing resource files from the repo.
    """
    repo_root = Path(__file__).parent.parent.resolve()

    output_dir = tmp_path / "gromacs"

    result = runner.invoke(
        app,
        [
            "build-gromacs",
            "--output-dir",
            str(output_dir),
            "--cation-file",
            str(repo_root / "resources/packmol/fe.pdb"),
            "--anion-file",
            str(repo_root / "resources/packmol/cl.pdb"),
            "--water-file",
            str(
                repo_root
                / "resources/packmol/gromacs/tip4p_2005/tip4p_2005_water.pdb"
            ),
            "--monomer-file",
            str(repo_root / "resources/packmol/fe_cl.pdb"),
            "--dimer-file",
            str(repo_root / "resources/packmol/fe_cl2.pdb"),
            "--x-box-length",
            "15",
            "--y-box-length",
            "15",
            "--z-box-length",
            "15",
            "--n-wat",
            "10",
            "--n-free-cations",
            "2",
            "--n-free-anions",
            "4",
            "--n-monomer",
            "1",
            "--n-dimer",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output

    g96 = output_dir / "system.g96"
    assert g96.exists()

    # Basic sanity check
    content = g96.read_text()
    assert "POSITION" in content


@pytest.mark.skipif(
    not (_has_command("packmol") and _has_command("moltemplate.sh")),
    reason="PACKMOL + moltemplate required",
)
def test_build_lammps_cli_real(tmp_path):
    """
    Full LAMMPS pipeline (PACKMOL + moltemplate).
    """
    repo_root = Path(__file__).parent.parent.resolve()

    output_dir = tmp_path / "lammps"

    result = runner.invoke(
        app,
        [
            "build-lammps",
            "--output-dir",
            str(output_dir),
            "--cation-file",
            str(repo_root / "resources/packmol/fe.xyz"),
            "--anion-file",
            str(repo_root / "resources/packmol/cl.xyz"),
            "--water-file",
            str(
                repo_root
                / "resources/packmol/lammps/tip4p_2005/tip4p_2005_water.xyz"
            ),
            "--monomer-file",
            str(repo_root / "resources/packmol/fe_cl.xyz"),
            "--dimer-file",
            str(repo_root / "resources/packmol/fe_cl2.xyz"),
            "--moltemplate-cation-file",
            str(repo_root / "resources/moltemplate/fe_ions.lt"),
            "--moltemplate-anion-file",
            str(repo_root / "resources/moltemplate/cl.lt"),
            "--moltemplate-water-file",
            str(repo_root / "resources/moltemplate/tip4p_2005.lt"),
            "--x-box-length",
            "15",
            "--y-box-length",
            "15",
            "--z-box-length",
            "15",
            "--n-wat",
            "10",
            "--n-free-cations",
            "2",
            "--n-free-anions",
            "4",
            "--n-monomer",
            "1",
            "--n-dimer",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output

    data_file = output_dir / "system.data"
    assert data_file.exists()

    # Sanity check
    content = data_file.read_text()
    assert "Atoms" in content