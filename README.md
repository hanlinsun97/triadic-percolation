# Percolation on networks with triadic interactions

Codes for triadic percolation on Poisson networks, scale-free networks and real networks.

Link to the paper: https://arxiv.org/pdf/2204.13067.pdf

This repository contains implementation of triadic percolation in both Julia and C:

- `functions.jl`: The codes include all the functions used for producing orbit diagrams from both theory and Monte Carlo simulation.
- `poi_simulation_theory.ipynb`: Examples of triadic percolation on Poisson structural networks with both theoretical and simulation results.
- `sf_simulation_theory.ipynb`: Examples of triadic percolation on scale-free structural networks with both theoretical and simulation results.


# How to use
The Julia code is implemented in Julia 1.8. Packages used are `LinearAlgebra`, `SpecialFunctions` and `Random`.

# Citing
If you find the code useful in your research, please cite the following paper:

```latex

@article{sun2022triadic,
  title={Triadic interactions induce blinking and chaos in the connectivity of higher-order networks},
  author={Sun, Hanlin and Radicchi, Filippo and Kurths, Juergen and Bianconi, Ginestra},
  journal={arXiv preprint arXiv:2204.13067},
  year={2022}
}
```
# License
This code can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  
This program is distributed by the authors in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
