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

The **ænet** code can be compiled automatically using the [CMake Build System](https://cmake.org/) or by manually adjusting platform-dependent Makefiles (see *Legacy Build* below).

Except for a number of Python scripts, **ænet** is developed in Fortran
95/2003. Generally, the source code is tested with the free GNU Fortran
compiler and the Intel Fortran compiler. While the **ænet** source code should be platform independent, it has only been tested on Linux and MacOS systems.

**ænet** requires three external libraries:

1.  BLAS (Basic Linear Algebra Subprograms),
2.  LAPACK (Linear Algebra PACKage),
3.  And the L-BFGS-B optimization routines by Nocedal et al.

Usually, an implementation of BLAS and LAPACK comes with the operating
system or the compiler. If that is not the case, the libraries can be
obtained from [Netlib.org](http://www.netlib.org/). `libblas.a` and
`liblapack.a` have to be in the system library path.

The L-BFGS-B routines, an implementation of the bounded limited-memory
Broyden-Fletcher-Goldfarb-Shanno algorithm, is distributed on the
[homepage of the
authors](http://www.ece.northwestern.edu/~nocedal/lbfgsb.html) (Nocedal
et al.). For the user’s convenience we have decided to distribute the
original L-BFGS-B files along with **ænet** package, so you do not have
to download the library yourself. However, when you use the `bfgs` training method with **ænet**, you should acknowledge the use of the L-BFGS-B library by citing:

R. H. Byrd, P. Lu and J. Nocedal, *SIAM J. Sci. Stat. Comp.* **16**
(1995) 1190-1208.

**ænet**’s Python interface further relies on
[NumPy](http://www.numpy.org) and on the [Atomic simulation
Environment](https://wiki.fysik.dtu.dk/ase), so these dependencies have
to be available when the **ænet** Python module is set up.


## Prerequisites

- **Fortran compiler** (e.g., gfortran, ifort)
- **BLAS and LAPACK** libraries (system, OpenBLAS, or Intel MKL)
- **MPI** (optional, for parallel builds)
- **OpenBLAS** or **Intel MKL** (optional, for optimized linear algebra)


## Building with CMake

> **Note:** When CMake is used, the L-BFGS-B library required by **ænet** is built automatically as part of the CMake process. No manual steps are needed.
 
### Short installation summary

In **ænet**'s root directory, run

```sh
$ mkdir build && cd build && cmake ../src/ && make
```

This will build a serial version of **æenet** with default options.  The main **ænet** executables, `generate.x`, `train.x`, and `predict.x` will be placed in the directory `./bin/`. Additional tools will be located in `./tools/` and the **ænet** library for interfaces with other simulations software will be in `./lib/`.

To build with support for parallelization via the message-passing interface (MPI), use

```sh
mkdir build && cd build && cmake ../src/ -DUSE_MPI=ON && make
```

Further options are described in the following.


### Additional Prerequisites

- **CMake** (version 3.15 or higher)

### Recommended Environment

- **MacOS (Apple Silicon or Intel):**
  - Compiler: `gfortran` (from Homebrew or MacPorts)
  - BLAS/LAPACK: Use the default (Apple Accelerate framework, detected automatically)
  - CMake options: No special options needed unless you want to enable MPI

- **Linux (Intel CPUs):**
  - Compiler: The Intel Fortran compiler is recommended, but the GNU Fortran compiler (`gfortran`) also performs well.
  - BLAS/LAPACK: The Intel Math Kernel Library (MKL) is strongly recommended.
  - CMake options:
    - For Intel MKL: `-DUSE_MKL=ON`
    - For MPI: `-DUSE_MPI=ON`

- **Linux (AMD or ARM CPUs):**
  - Compiler: The GNU Fortran compiler (`gfortran`) is widely supported and recommended.
  - BLAS/LAPACK:
    - The OpenBLAS library is recommended (`-DUSE_OPENBLAS=ON`)
    - Otherwise, system BLAS/LAPACK can be used
  - CMake options:
    - For OpenBLAS: `-DUSE_OPENBLAS=ON`
    - For MPI: `-DUSE_MPI=ON` (if you want parallel execution)

### Configuration Options

The CMake build supports the following **ænet**-specific options:

- `USE_MPI` (default: OFF): Enable MPI parallelization
- `USE_OPENBLAS` (default: OFF): Use OpenBLAS instead of system BLAS/LAPACK
- `USE_MKL` (default: OFF): Use Intel MKL instead of system BLAS/LAPACK
- `CMAKE_BUILD_TYPE` (default: Release): Set to `Debug` for debug builds

Set these options with the `-D` flag when running CMake, e.g.:

```sh
cmake ../src/ -DUSE_MPI=ON -DUSE_OPENBLAS=ON -DCMAKE_BUILD_TYPE=Debug
```

### Detailed Build Instructions

#### 1. **Create a build directory in the project root and configure the project**

From the project's root directory (where this `README.md` is located), run:

```sh
mkdir build
cd build
cmake ../src
```

This creates a `build` directory separate from the source code (in `./src/`) and configures the build inside it. Add any configuration options to the `cmake` command as needed (see above), for example: `cmake ../src -DUSE_MPI=ON`.

> Note: The name of the build directory is arbitrary. Multiple build directories can be created to build multiple versions of **ænet** (e.g., serial, parallel, and debug versions).

#### 2. **Build the main executables and libraries (from within the `build` directory)**

```sh
make
```

Or build specific targets, e.g.:

```sh
make main    # Builds generate.x train.x predict.x
make tools   # Builds all tool executables
make lib     # Builds static and shared libraries
```

This will build the following:

- **`./bin/`:** `generate.x`, `train.x`, `predict.x`
- **`./tools/`:** `fingerprint.x`, `neighbors.x`, `trnset_info.x`, `trnset2ASCII.x`
- **`./lib/`:** `libaenet.a` (static), `libaenet.so` or `libaenet.dylib` (shared)


#### 3. **(Optional) Run the unit tests**

Build and run all unit tests with:

```sh
make test
```

## Legacy Makefile-Based Build


### Short installation summary

1.  Compile the L-BFGS-B library

      - Enter the directory “./lib”

        `$ cd ./lib`

      - Adjust the compiler settings in the “Makefile”

      - Compile the library with

        `$ make`

    The library file `liblbfgsb.a`, required for compiling **ænet**,
    will be created.

2.  Compile the **ænet** package

      - Enter the directory “./src”

        `$ cd ./src`

      - Compile the ænet source code with

        `$ make -f makefiles/Makefile.XXX`

        where `Makefile.XXX` is an approproiate Makefile.

        To see a list of available Makefiles just type:

        `$ make`

    The following executables will be generated in “./bin”:

      - `generate.x`: generate training sets from atomic structure files
      - `train.x`: train new neural network potentials
      - `predict.x`: use existing ANN potentials for energy/force
        prediction

    *Alternatively, you may build ænet using CMake. See below for details.*

3.  (Optional) Install the Python interface

      - Enter the directory “./python”

        `$ cd ./python`

      - Install the Python module with

        `$ python setup.py install --user`

    This will set up the Python **ænet** module for the current user,
    and it will also install the user scripts `aenet-predict.py` and
    `aenet-md.py`.

### Detailed installation instructions

#### Compilation of external libraries that are distributed with **ænet**

All external libraries needed by the ænet code are in the directory
“./lib”. Currently, only one external library is distributed with
**ænet**, the L-BFGS-B library (see above).

To compile the external libraries

1.  Enter the directory “./lib”

    `$ cd ./lib`

2.  Adjust the compiler settings in the “Makefile”

    The Makefile contains settings for the GNU Fortran compiler
    (`gfortran`) and the Intel Fortran compiler (`ifort`). Uncomment the
    section that is appropriate for your system.

3.  Compile the library with

    `$ make`

The static library “liblbfgsb.a”, required to build **ænet**, will be
created.

#### Build **ænet**

The **ænet** source code is located in “./src”.

1.  Enter “./src”

    `$ cd ./src`

2.  To see a short explanation of the Makefiles that come with **ænet**,
    just run `make` without any options.

    `$ make`

    Select the Makefile that is appropriate for your computer.

3.  Compile with

    `$ make -f makefiles/Makefile.XXX`

    where `Makefile.XXX` is the selected Makefile.

Three executables will be generated and stored in “./bin”:

  - `generate.x`: generate training sets from atomic structure files
  - `train.x`: train new neural network potentials
  - `predict.x`: use existing ANN potentials for energy/force prediction

#### Set up the Python interface

1.  Enter the directory “./python”

    `$ cd ./python`

2.  Install the Python module with

    `$ python setup.py install --user`

This will set up the Python **ænet** module for the current user, and it
will also install the user scripts `aenet-predict.py` and `aenet-md.py`.
