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
- **Fortran compiler** (e.g., `gfortran`, `ifort`, `ifx`)
- **BLAS and LAPACK** libraries
- **MPI** (optional, for parallel builds)

## Building with CMake

All commands should be run from the root directory of the **ænet** project.

### Build Options

Print available build options with

```sh
cmake .
```

This will print

```sh
-- ==============================================================================
-- Build aenet with CMake
-- ==============================================================================
-- To build:  cmake -S . -B build -DBUILD_AENET=ON && cmake --build build --target build_all
-- Or pick a different target: main | lib | tools | build_tests | test
--
-- The following options are available to configure the build:
--
--   -DCMAKE_BUILD_TYPE=Release/Debug  (default: Release)
--      Build type.
--   -DUSE_MPI=ON/OFF                  (default: OFF)
--      Enable MPI parallelization.
--   -DUSE_MKL=ON/OFF                  (default: OFF)
--      Use Intel MKL for BLAS/LAPACK.
--   -DUSE_OPENBLAS=ON/OFF             (default: OFF)
--      Use OpenBLAS for BLAS/LAPACK.
--   -DCMAKE_Fortran_COMPILER=<path>   (e.g., ifort, gfortran)
--      Specify the Fortran compiler.
--
-- Example: cmake .. -DBUILD_AENET=ON -DUSE_MPI=ON -DCMAKE_BUILD_TYPE=Debug
--
-- The following custom targets are available:
--
--   make main         (Builds generate.x, train.x, predict.x)
--   make lib          (Builds aenet static and shared libraries)
--   make tools        (Builds all tool executables)
--   make build_all    (Builds all executables, libraries, and tests)
--   make build_tests  (Builds all unit test executables)
--   make test         (Runs the unit tests)
--
-- Configuring done (0.0s)
-- Generating done (0.0s)
```

### Example Build Configurations

In all examples, a new directory named `build` will be created, containing the CMake configuration.  This directory can be safely removed to reset to the initial state.  All commands should be run from the root directory of the **ænet** project.

1. Let CMake decide all parameters:

```sh
cmake  -S . -B build -DBUILD_AENET=On
```

2. Turning on MPI parallelization

```sh
cmake  -S . -B build -DBUILD_AENET=On -DUSE_MPI=On
```

3. Selecting Intel's `ifx` compiler with Intel MPI and Math Kernel Library (MKL)

```sh
cmake  -S . -B build -DBUILD_AENET=On -DUSE_MPI=On -DUSE_MKL=On -DCMAKE_Fortran_COMPILER=ifx
```

4. GNU Fortran compiler with OpenBLAS library

```sh
cmake  -S . -B build -DBUILD_AENET=On -DUSE_OPENBLAS=On -DCMAKE_Fortran_COMPILER=gfortran
```

### Build the Code

Build everything, including the main executables, the library, and the tools with

```sh
cmake --build build --target build_all
```

Build only the executables `generate.x`, `train.x`, and `predict.x` with

```sh
cmake --build build --target main
```

or `cmake --build build --target  <target>` where `<target>` is one of the following:

- `build_all`: Builds all executables, libraries, and tools.
- `main`: Builds the main executables (`generate.x`, `train.x`, `predict.x`).
- `lib`: Builds the `aenet` static and shared libraries.
- `tools`: Builds the auxiliary tool executables.
- `build_tests`: Builds the unit test executables.
- `test`: Runs the unit tests (after they have been built).

The resulting files will be placed in the `bin`, `lib`, and `tools` directories in the `build` directory:

- `build/bin`: `generate.x`, `train.x`, and `predict.x`
- `build/lib`: system-dependent, e.g., `libaenet.so` and `libaenet.a`
- `build/tools`: `fingerprint.x`, `neighbors.x`, `trnset_info.x`, and `trnset2ASCII.x`

### (Optional) Run the Unit Tests

To confirm that the **ænet** command-line tools are working as expected, you can run the test suite with

```sh
ctest --test-dir build --output-on-failure
```

### (Optional) Installing the ænet Files

To install the **ænet** files in the corresponding directories of the root directory, run

```sh
cmake --install build
```

this will install the files in the `bin`, `lib`, and `tools` directories in
the **ænet** project root.  The installation location can be changed with
`--prefix` flag

```sh
cmake --install build --prefix /custom/installation/path
```

### (Optional) Creating a Debug Build

To compile with debug flags (e.g., for checking array bounds and backtraces), set `CMAKE_BUILD_TYPE` to `Debug`:

```sh
cmake .. -DBUILD_AENET=ON -DCMAKE_BUILD_TYPE=Debug ...
```

The resulting executables will have a `_debug` suffix.

## Makefile-based Installation

### 1. Compile the L-BFGS-B library

Enter the directory “./lib”

    $ cd ./lib

Adjust the compiler settings in the “Makefile”

Compile the library with

    $ make

The library file `liblbfgsb.a`, required for compiling ænet, will be created.

### 2. Compile the ænet package

Enter the directory “./src”

    $ cd ./src

Compile the **ænet** source code with

    $ make -f makefiles/Makefile.XXX

where Makefile.XXX is an approproiate Makefile.

To see a list of available Makefiles just type:

    $ make

The following executables will be generated in “./bin”:

- `generate.x`: generate training sets from atomic structure files
- `train.x`: train new neural network potentials
- `predict.x`: use existing ANN potentials for energy/force prediction

### 3. Compile Optional Components

Compile the **ænet** tools with

    $ make tools -f makefiles/Makefile.XXX

Compile the **ænet** libraries (shared and static) with

    $ make lib -f makefiles/Makefile.XXX

## Installing the Python Interface

**Note:** We now recommend using [aenet-python](https://github.com/atomisticnet/aenet-python).

After building the Fortran libraries, you can install a basic Python interface with

```sh
cd ../python3
python setup.py install --user
```

This will install the `aenet` Python module and the scripts `aenet-predict.py` and `aenet-md.py`.
