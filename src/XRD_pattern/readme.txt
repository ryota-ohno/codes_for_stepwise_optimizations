This is workflow where we calculate XRD pattern with optimized lattice parameters.
We also prepared the code for experimental structures so that we can compare with each other.

The file format is .ipynb for Jupyter notebook since we use this code just for visuallizing.
You excuse each cell in the order then you input your parameters and molecular coordinates files.
Parameters are shown in our paper so we add molecular coordinates file (.xyz).

We first calculate reciplocal lattice vector from optimized or experimental lattice vector.
We also calculate atom scattering factors and the sturucture factor.
By combining these values and sharp Laue's function where N=20000 we can calculate xrd powder pattern.
Lorentz's polarization factor is added in this calculation.

directry tree
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
This workflow was developved by Mr. Arai and Mr. Koyama.