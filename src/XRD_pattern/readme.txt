This workflow calculates X-ray diffraction (XRD) powder patterns using optimized lattice parameters. We also provide a workflow for experimental crystal structures, allowing direct comparison between calculated and experimental XRD patterns.

The workflow is provided as Jupyter notebooks (.ipynb), as the code is primarily intended for visualization and data analysis. Simply execute the notebook cells in order, then provide your input parameters and molecular coordinate files.

The optimized lattice parameters used in our study are reported in our paper. We also include the corresponding molecular coordinate files (.xyz).

The workflow performs the following steps:
1. Calculates the reciprocal lattice vectors from the optimized or experimental lattice vectors.
2. Computes atomic scattering factors and the crystal structure factor.
3. Combines these quantities with a finite-size Laue function (N = 20,000) to calculate the XRD powder pattern.
4. Applies the Lorentz–polarization correction factor to obtain the final diffraction intensity.

Directry structures
XRD_pattern/
├── src/
│   └── powder_pattern_model.ipynb
├── opt_model/
│   └── molecular coordinate files(.xyz) of dimer in the optimized molecular arrangements
├── cryst/
│   └── molecular arrangements(.xyz) extracted from crystal information file.
└── atom_scattering_factor_list
    └── X.txt contains coeficients of atom scattering factors.

Acknowledgement
The workflow included in this repository is based on code originally developed by Dr. Shunto Arai and Mr. Kanata Koyama. 
We gratefully acknowledge their contribution to the development of the XRD calculation workflow.