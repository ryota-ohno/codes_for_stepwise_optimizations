# codes_for_stepwise_optimizations

This repository contains the workflow developed for crystal structure prediction based on stepwise optimization of intermolecular interactions.

The workflow accompanies the manuscript:

> Stepwise quantum chemical optimization reveals the origin of layered herringbone packing and polymorphism in polyacenes
>
> Ryota Ono et al.

The workflow is designed to optimize molecular packing by separating the optimization of intralayer and interlayer structures, enabling efficient crystal structure prediction for organic molecular crystals.
Sub-workflows are also prepared: transfer integral calculations for multiple molecular arrangements, XRD-powder-pattern calculation and total-intermolecular-interaction calculations 

---

## Features

- Stepwise optimization of molecular packing
- Transfer integral calculations 
- XRD-powder-pattern calculations
- Calculation of total intermolecular interaction energy
- atomic coordinates of polyacene monomers

---

## Repository structure

```text
stepwise_optimizations/

tcal_csv/

XRD_pattern/

tot_energy/

monomer/

```

Detailed instructions for each module are provided in the corresponding subdirectory README files.

---

## Installation

Install the required Python packages

```bash
pip install -r requirements.txt
```

or

```bash
pip install .
```

---

## Usage

Each module contains its own documentation.

For example,

```
stepwise_optimizations/readme.txt
tcal_csv/readme.txt
XRD_pattern/readme.txt
tot_energy/readme.txt
```

describe the required input files, command-line options, and output formats.

---

## Citation

If you use this workflow in your research, please cite both the accompanying paper and the software archived on Zenodo.

The citation information is also available in `CITATION.cff`.

---

## License

This project is distributed under the MIT License.
