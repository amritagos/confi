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