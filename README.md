# What is **ænet**?

The Atomic Energy NETwork (**ænet**) package (http://ann.atomistic.net) is a collection of tools for the construction and application of atomic interaction potentials based on artificial neural networks (ANN). The **ænet** code allows the accurate interpolation of structural energies, e.g., from electronic structure calculations, using ANNs. ANN potentials generated with **ænet** can then be used in larger scale atomistic simulations and in situations where extensive sampling is required, e.g., in molecular dynamics or Monte-Carlo simulations.

# License

Copyright (C) 2012-2025 Nongnuch Artrith (nartrith@atomistic.net)

The **ænet** source code is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, you can obtain one at <http://mozilla.org/MPL/2.0/>.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Mozilla
Public License, v. 2.0, for more details.

# Installation

The **ænet** code is built using the [CMake Build System](https://cmake.org/).

Except for a number of Python scripts, **ænet** is developed in Fortran 95/2003 and has been tested with the GNU Fortran and Intel Fortran compilers on Linux and MacOS systems.

**ænet** requires the following external libraries:

1.  **BLAS** (Basic Linear Algebra Subprograms)
2.  **LAPACK** (Linear Algebra PACKage)
3.  **L-BFGS-B** optimization routines by Nocedal et al.

A BLAS/LAPACK implementation is usually provided by the operating system or compiler. Alternatively, high-performance libraries like [OpenBLAS](https://www.openblas.net/) or [Intel MKL](https://software.intel.com/en-us/mkl) can be used.

The L-BFGS-B sources are included with **ænet** and are compiled automatically during the build process. When using the `bfgs` training method, please cite:

> R. H. Byrd, P. Lu and J. Nocedal, *SIAM J. Sci. Stat. Comp.* **16** (1995) 1190-1208.

The Python interface also requires [NumPy](http://www.numpy.org) and the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase).

## Prerequisites

- **CMake** (version 3.15 or higher)
- **Fortran compiler** (e.g., `gfortran`, `ifort`)
- **BLAS and LAPACK** libraries
- **MPI** (optional, for parallel builds)

## Building with CMake

The build process is designed to be straightforward. All commands should be run from the root directory of the **ænet** project.

### 1. Create a Build Directory

It is best practice to create a separate directory for the build.

```sh
mkdir build
cd build
```

### 2. Configure the Build

To see a list of build options (such as MPI), run `cmake` from the `build` directory on the parent directory without any flags

```sh
# from the build directory
cmake ..
```

To configure ænet with standard options, run `cmake` with the flag `-DBUILD_AENET=ON` 

```sh
# from the build directory
cmake .. -DBUILD_AENET=ON
```

This will configure a standard **serial** build using the default Fortran compiler and system BLAS/LAPACK libraries.

### 3. Build the Code

After configuration, build ænet with

```sh
cmake --build . --target <target>
```

or `make <target>` where `<target>` is one of the following:

- `build_all`: Builds all executables, libraries, and tools.
- `main`: Builds the main executables (`generate.x`, `train.x`, `predict.x`).
- `lib`: Builds the `aenet` static and shared libraries.
- `tools`: Builds the auxiliary tool executables.
- `build_tests`: Builds the unit test executables.
- `test`: Runs the unit tests (after they have been built).

The resulting files will be placed in the `bin`, `lib`, and `tools` directories in the `build` directory.  To install them in the corresponding directories of the root directory, run

```sh
cmake --build . --target install
```

or `make install`.

### Common Build Scenarios

#### Building a Parallel Version with MPI

To build the MPI-parallelized version of **ænet**, enable the `USE_MPI` option during the configuration step:

```sh
cmake .. -DBUILD_AENET=ON -DUSE_MPI=ON
```

Then, run `make <target>` as usual.

#### Creating a Debug Build

To compile with debug flags (e.g., for checking array bounds and backtraces), set `CMAKE_BUILD_TYPE` to `Debug`:

```sh
cmake .. -DBUILD_AENET=ON -DCMAKE_BUILD_TYPE=Debug
```

The resulting executables will have a `_debug` suffix.

#### Specifying a Compiler

CMake will automatically find a Fortran compiler. To specify a different one (e.g., `ifort`), use the `CMAKE_Fortran_COMPILER` variable:

```sh
cmake .. -DBUILD_AENET=ON -DCMAKE_Fortran_COMPILER=ifort
```

#### Using a Specific BLAS/LAPACK Library

You can bias the build system to use a vendor-specific high-performance library for BLAS/LAPACK:

- **Intel MKL:**
  ```sh
  cmake .. -DBUILD_AENET=ON -DUSE_MKL=ON
  ```

- **OpenBLAS:**
  ```sh
  cmake .. -DBUILD_AENET=ON -DUSE_OPENBLAS=ON
  ```

### Installing the Python Interface

After building the Fortran libraries, you can install the Python interface:

```sh
cd ../python3
python setup.py install --user
```

This will install the `aenet` Python module and the scripts `aenet-predict.py` and `aenet-md.py`.
