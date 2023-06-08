# dune-surfaces_int

Code to produce the results of the paper.

```bibtex
@article{...,
  title={HIGH-ORDER-INTEGRATION FOR REGULAR EMBEDDED MANIFOLDS},
  author={Zavalani, G., Sander, O. and Hecht, M.},
  journal={arXiv preprint arXiv:...},
  year={2023}
}
```


## Preparing the Sources
This project is developed in the [DUNE](https://dune-project.org) framework.

You need to install the following dependencies, all with the branch `releases/2.7`

- [dune-common](https://gitlab.dune-project.org/core/dune-common)
- [dune-geometry](https://gitlab.dune-project.org/core/dune-geometry)
- [dune-localfunctions](https://gitlab.dune-project.org/core/dune-localfunctions)
- [dune-istl](https://gitlab.dune-project.org/core/dune-istl)
- [dune-grid](https://gitlab.dune-project.org/core/dune-grid)
- [dune-typetree](https://gitlab.dune-project.org/staging/dune-typetree)
- [dune-functions](https://gitlab.dune-project.org/staging/dune-typetree)
- [dune-foamgrid](https://gitlab.dune-project.org/extensions/dune-foamgrid)
- [dune-vtk](https://gitlab.dune-project.org/extensions/dune-vtk)
- [dune-curvedgeometry](https://gitlab.mn.tu-dresden.de/iwr/dune-curvedgeometry)
- [dune-curvedgrid](https://gitlab.mn.tu-dresden.de/iwr/dune-curvedgrid)
- [dune-gmsh4](https://gitlab.mn.tu-dresden.de/iwr/dune-gmsh4)

Put all these modules together with this project into a common base directory, called `my_dune`
in the following. This module is assumed to be in a subdirectory called `dune-surface_int`.

As additional external dependencies we have

- CMake >= 3.13
- PkgConfig
- SuiteSparse

that can be installed as packages in your Linux distribution.


## Compiling the Code
If these preliminaries are met, you should run

```bash
cd my_dune
dunecontrol cmake -DCMAKE_BUILD_TYPE=Release
dunecontrol make
```

to configure and compile all modules, including this one. All CMake and make output is
put into the build directory `build-cmake/`.

For further information about the installation process, see also the official documentation
https://dune-project.org/doc/installation/ of DUNE.

## Run the Examples from dune-surface_int
In order to execute the code to run the numerical examples for the biconcave shape, run


```bash
cd my_dune/dune-surface_int./build-cmake/src/biconcave_disc
```
 for the double torus, run

```bash
cd my_dune/dune-surface_int./build-cmake/src/double_torus
```
 for the genus two surface, run

```bash
cd my_dune/dune-surface_int./build-cmake/src/genus_surface
```
