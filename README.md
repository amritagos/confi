# confi
A package that can create configurations and data files for LAMMPS and GROMACS.

## Installation

In the top-level directory: 

```bash
micromamba create -f environment.yml # for the first time to create the environment 
micromamba activate confienv # always 
```

Install the python package `confi` using `pip`: 

```bash
pip install .
```

The folder `snmk_worflow` contains an example of a workflow that creates LAMMPS data files. 
In order to use it: 

```bash
cd snmk_workflow
snakemake --cores 1 -n # dry run 
snakemake --cores 1 # actual run 
```

### External Dependencies

The following are external dependencies that need to be in your path: 

- moltemplate
- packmol

I recommend installing them via `micromamba` or `conda`.

## CLI usage

As an alternative to the Snakemake workflows, Confi also provides a CLI.

You can see the available commands with:

```bash
confi --help
python -m confi --help
```

### Quick examples using the CLI

These examples mirror the Snakemake workflows but can be run directly.

#### GROMACS example

This generates a `system.g96` file using PACKMOL and ASE.

```bash
confi build-gromacs \
  --output-dir results/gromacs/example \
  --cation-file resources/packmol/fe.pdb \
  --anion-file resources/packmol/cl.pdb \
  --water-file resources/packmol/gromacs/tip4p_2005/tip4p_2005_water.pdb \
  --monomer-file resources/packmol/fe_cl.pdb \
  --dimer-file resources/packmol/fe_cl2.pdb \
  --x-box-length 20 \
  --y-box-length 20 \
  --z-box-length 20 \
  --n-wat 50 \
  --n-free-cations 5 \
  --n-free-anions 10 \
  --n-monomer 2 \
  --n-dimer 1
```

#### LAMMPS example

```bash
confi build-lammps \
  --output-dir results/lammps/example \
  --cation-file resources/packmol/fe.xyz \
  --anion-file resources/packmol/cl.xyz \
  --water-file resources/packmol/lammps/tip4p_2005/tip4p_2005_water.xyz \
  --monomer-file resources/packmol/fe_cl.xyz \
  --dimer-file resources/packmol/fe_cl2.xyz \
  --moltemplate-cation-file resources/moltemplate/fe_ions.lt \
  --moltemplate-anion-file resources/moltemplate/cl.lt \
  --moltemplate-water-file resources/moltemplate/tip4p_2005.lt \
  --x-box-length 20 \
  --y-box-length 20 \
  --z-box-length 20 \
  --n-wat 50 \
  --n-free-cations 5 \
  --n-free-anions 10 \
  --n-monomer 2 \
  --n-dimer 1
```