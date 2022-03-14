# What is **ænet**?

<span id="sec:about"></span>

The Atomic Energy NETwork (**ænet**) package (http://ann.atomistic.net) is a collection of tools
for the construction and application of atomic interaction potentials
based on artificial neural networks (ANN). The **ænet** code allows the
accurate interpolation of structural energies, e.g., from electronic
structure calculations, using ANNs. ANN potentials generated with
**ænet** can then be used in larger scale atomistic simulations and in
situations where extensive sampling is required, e.g., in molecular
dynamics or Monte-Carlo simulations.

# License

Copyright (C) 2012-2022 Nongnuch Artrith (nartrith@atomistic.net)

The **aenet** source code is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at <http://mozilla.org/MPL/2.0/>.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Mozilla
Public License, v. 2.0, for more details.

# Installation

<span id="sec:installation"></span>

## Short installation summary

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

3.  (Optional) Install the Python interface
    
      - Enter the directory “./python”
        
        `$ cd ./python`
    
      - Install the Python module with
        
        `$ python setup.py install --user`
    
    This will set up the Python **ænet** module for the current user,
    and it will also install the user scripts `aenet-predict.py` and
    `aenet-md.py`.

## Detailed installation instructions

Except for a number of Python scripts, **ænet** is developed in Fortran
95/2003. Generally, the source code is tested with the free GNU Fortran
compiler and the commercial Intel Fortran compiler, and the Makefile
settings for these two compilers are provided. While the **ænet** source
code should be platform independent, we mainly target Linux and Unix
clusters and **ænet** has not been tested on other operating systems.

**ænet** requires three external libraries:

1.  BLAS (Basic Linear Algebra Subprograms),
2.  LAPACK (Linear Algebra PACKage),
3.  And the L-BFGS-B optimization routines by Nocedal et al.

Usually, some implementation of BLAS and LAPACK comes with the operating
system or the compiler. If that is not the case, the libraries can be
obtained from [Netlib.org](http://www.netlib.org/). `libblas.a` and
`liblapack.a` have to be in the system library path in order to compile
**ænet**.

The L-BFGS-B routines, an implementation of the bounded limited-memory
Broyden-Fletcher-Goldfarb-Shanno algorithm, is distributed on the
[homepage of the
authors](http://www.ece.northwestern.edu/~nocedal/lbfgsb.html) (Nocedal
et al.). For the user’s convenience we have decided to distribute the
original L-BFGS-B files along with **ænet** package, so you do not have
to actually download the library yourself. However, each application of
**ænet** should also acknowledge the use of the L-BFGS-B library by
citing:

R. H. Byrd, P. Lu and J. Nocedal, *SIAM J. Sci. Stat. Comp.* **16**
(1995) 1190-1208.

**ænet**’s Python interface further relies on
[NumPy](http://www.numpy.org) and on the [Atomic simulation
Environment](https://wiki.fysik.dtu.dk/ase), so these dependencies have
to available when the **ænet** Python module is set up.

### Compilation of external libraries that are distributed with **ænet**

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

### Build **ænet**

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

### Set up the Python interface

1.  Enter the directory “./python”
    
    `$ cd ./python`

2.  Install the Python module with
    
    `$ python setup.py install --user`

This will set up the Python **ænet** module for the current user, and it
will also install the user scripts `aenet-predict.py` and `aenet-md.py`.
