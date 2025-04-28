from pydantic import BaseModel, model_validator
from pydantic.fields import Field
from typing import Optional, Self
from pathlib import Path
from enum import Enum


class FileType(Enum):
    """Input file types allowed for PACKMOL

    Args:
        Enum (str): File extensions for allowed input file types in PACKMOL
    """

    XYZ = "xyz"
    PDB = "pdb"


class PackmolInput(BaseModel):
    filetype: FileType = FileType.XYZ  # Can be either "xyz" or "pdb"
    cation_file: Optional[Path] = None
    anion_file: Optional[Path] = None
    water_file: Optional[Path] = None


class PackmolParams(BaseModel):

    packmol_input: PackmolInput
    system_file: Path  # output path to which PACKMOL will create the file if you actually run it with the created input file
    x_box_length: float = Field(gt=0)  # Box length in the x-dimension
    y_box_length: float = Field(gt=0)  # Box length in the y-dimension
    z_box_length: float = Field(gt=0)  # Box length in the z-dimension
    n_wat: int = Field(ge=0, default=0)  # Number of water molecules
    water_file: Optional[Path] = (
        None  # if n_wat is greater than 0 then this cannot be None
    )
    n_free_cations: int = Field(ge=0, default=0)
    cation_file: Optional[Path] = (
        None  # if n_free_cations is greater than 0 then this cannot be None
    )
    cation_radius: float = Field(gt=0, default=1.75)  # overrides the tolerance
    n_free_anions: int = Field(ge=0, default=0)  # free anions
    anion_file: Optional[Path] = (
        None  # if n_free_cations is greater than 0 then this cannot be None
    )
    anion_radius: float = Field(gt=0, default=1.75)  # overrides the tolerance
    n_monomer: int = Field(ge=0, default=0)  # Number of cation-anion units
    monomer_file: Optional[Path] = (
        None  # if n_monomer is greater than 0 you must input this or PACKMOL will fail
    )
    n_dimer: int = Field(ge=0, default=0)  # Number of anion-cation-anion units
    dimer_file: Optional[Path] = (
        None  # if n_dimer is greater than 0 you must input this or PACKMOL will fail
    )
    tolerance: float = Field(
        gt=0, default=2.0
    )  # minimum distance between subunits in PACKMOL. can be overriden by radius
    seed: int = (
        -1
    )  # PACKMOL uses -1 to generate a seed automatically from the computer time.

    @model_validator(mode="after")
    def check_required_files(self) -> Self:
        errors = []
        # Check free cations in packmol_input
        if self.n_free_cations > 0 and self.packmol_input.cation_file is None:
            errors.append(
                "n_free_cations > 0 but packmol_input.cation_file has not been provided."
            )
        # Check free anions in packmol_input
        if self.n_free_anions > 0 and self.packmol_input.anion_file is None:
            errors.append(
                "n_free_anions > 0 but packmol_input.anion_file has not been provided."
            )
        # Check monomer_file if n_monomer > 0
        if self.n_monomer > 0 and self.monomer_file is None:
            errors.append("n_monomer > 0 but monomer_file has not been provided.")
        # Check dimer_file if n_dimer > 0
        if self.n_dimer > 0 and self.dimer_file is None:
            errors.append("n_dimer > 0 but dimer_file has not been provided.")

        # Check water_file if n_wat>0
        if self.n_wat > 0 and self.packmol_input.water_file is None:
            errors.append("n_wat > 0 but water_file has not been provided.")

        if errors:
            raise ValueError("; ".join(errors))
        return self


class MoltemplateInput(BaseModel):
    cation_file: Optional[Path] = None
    anion_file: Optional[Path] = None
    water_file: Optional[Path] = None


class MoltemplateParams(BaseModel):

    moltemplate_input: MoltemplateInput
    n_free_cations: int = Field(ge=0, default=0)
    n_free_anions: int = Field(ge=0, default=0)  # free anions
    n_wat: int = Field(ge=0, default=0)  # Number of water molecules
    n_monomer: int = Field(
        ge=0, default=0
    )  # Number of cation-anion units (should be inside the cation file)
    n_dimer: int = Field(
        ge=0, default=0
    )  # Number of anion-cation-anion units (definition should be inside cation file)
    x_box_length: float = Field(gt=0)  # Box length in the x-dimension
    y_box_length: float = Field(gt=0)  # Box length in the y-dimension
    z_box_length: float = Field(gt=0)  # Box length in the z-dimension

    @model_validator(mode="after")
    def check_required_files(self) -> Self:
        errors = []
        # Check free cations
        if self.n_free_cations > 0 and self.moltemplate_input.cation_file is None:
            errors.append(
                "n_free_cations > 0 but moltemplate_input.cation_file has not been provided."
            )
        # Check free anions
        if self.n_free_anions > 0 and self.moltemplate_input.anion_file is None:
            errors.append(
                "n_free_anions > 0 but moltemplate_input.anion_file has not been provided."
            )

        # Check water_file if n_wat>0
        if self.n_wat > 0 and self.moltemplate_input.water_file is None:
            errors.append("n_wat > 0 but water_file has not been provided.")

        if errors:
            raise ValueError("; ".join(errors))
        return self
