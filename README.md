# ManifoldInterpolation

[![Documentation](https://github.com/Berger-Stefan/ManifoldInterpolation/actions/workflows/documentation.yml/badge.svg)](https://berger-stefan.github.io/ManifoldInterpolation/dev/)

This Julia package was created for my bachelor thesis why can be found under :
```
https://git.rwth-aachen.de/berger.st.11.11/bachelor_thesis
```
## Features
- loading of meshes in the vtk format
- converting of meshes to graphs
- computation of the eigenvalues and eigenvectors of the Laplace operator applied on graphs
- interpolation of the eigenvectors in a one dimensional parameter domain
- export of graphs and the associated eigenvalues and eigenvectors in the vtk format

## Installation
run the following command inside the Julia REPL:
```
    import Pkg
    Pkg.add("https://github.com/Berger-Stefan/ManifoldInterpolation")
```
