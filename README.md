# Semi-Classical 2Spin System

This repository contains C++ and Python code for analysis of a semi-classical two-spin system. We consider the Hamiltonian

$$H_0 = \sum_{\alpha \in\{x,y,z\}} \left[-J_{\alpha} S_{1\alpha} S_{2\alpha} + \frac{1}{2} A_\alpha \left(S_{1\alpha}^2 + S_{2\alpha}^2\right)\right],$$

where 

$$\mathbf{S}_{1} = (S_{1x},S_{1y},S_{1z})$$

and 

$$\mathbf{S}_{2} = (S_{2x},S_{2y},S_{2z})$$ 

are spin- $S$ degrees of freedom in the quantum version, and are unit vectors in the classical case. Depending on the choice of parameters 

$$\mathbf{J} = (J_x,J_y,J_z)$$ 

and 

$$\mathbf{A}=(A_x,A_y,A_z),$$ 

this model can be integrable or chaotic.

## Speed-Fisher Information

The chaotic/non-chaotic nature of the system is captured by the so called "speed-Fisher information". This quantity can be estimated by measuring thermodynamic drag coefficient under a cyclic perturbation of the system. Concretely, in a Gibbs state, the speed-Fisher information is given by

$$\mathscr{I}_{\bar{v}}(\mu) = \frac{2}{k_B T \bar{v}^2}\langle\Delta E\rangle,$$

where $\bar{v}$ is the average speed of the protocol and $\langle\Delta E\rangle$ is the average heating induced by the perturbation. Specifically, here we consider the perturbation $\lambda(t) V$, where

$$V = S_{1z}S_{2z}$$

and

$$\lambda(t) = \begin{cases}
\frac{\pi\bar{v}}{2\mu}\sin^2 \mu t, & 0 \leq t \leq \frac{2\pi}{\mu} \\ 
0, & \text{otherwise.}
\end{cases}$$

## Repository structure

- `src/` - C++ source code
- `include/` - C++ header files
- `python/` - Python analysis and plotting scripts
- `data/` - dataset and parameter files
- `build/` - compiled object files and binaries (ignored by git)
- `figures/` - plots and generated visuals
- `CITATION.cff` - citation metadata

## Build

Run:

```sh
make
```

Binaries will be created under `build/bin/`.

## Clean

```sh
make clean
```

## Notes

- External dependencies such as Eigen and nlohmann/json are expected under `../Libraries/`.
- Add Python dependencies to `python/requirements.txt` if required.
