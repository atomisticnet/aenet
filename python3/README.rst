------------------------------------------------------------------------
            Python interface to the symmetry function module
------------------------------------------------------------------------

Prerequisites:

  Python 3.6
  Cython


Installation:

(1) Compile the symmetry function C library in ../src

$ cd ../src
$ make -f ./makefiles/Makefile.XXX lib
$ cd -

Select the Makefile.XXX that is appropriate for your system.

(2) Build Python extension module

$ python setup.py build_ext --inplace
