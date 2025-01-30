# NativeSolid / RBM
Native solid solver for [CUPyDO](https://github.com/ulgltas/CUPyDO)

## Features
Compute static or dynamic displacements for several 2D configurations.

* Meshes
  - .su2 format (can be built and exported from gmsh)
* Configuration
  - Airfoil (pitch/plunge)
  - Horizontal spring
  - Vertical spring
* Dynamic computation
  - AlphaGen time integration
  - Runge-Kutta (order 4) time integration
  - Mass/damping/stiffness/nonlinear terms
* Static Computation
  - Stiffness term only

## Build

Required packages
```
sudo apt install build-essential cmake libopenblas-dev liblapacke-dev libpython3-dev swig
```

```
mkdir build && cd build
cmake ..
make -j 4
```

