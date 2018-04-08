Chebyshev Descriptor
====================

This directory contains an example that showcases the Chebyshev
descriptor for local atomic environments proposed in
[*Phys. Rev. B* **96** (2017)
014112](http://dx.doi.org/10.1103/PhysRevB.96.014112).

Owing for the large amount of data the reference data set and the
training set files are not included in the example.  Otherwise the
example is self-contained.  The reference data set used is the
TiO<sub>2</sub> data set from [*Comput. Mater. Sci.* **114** (2016)
135-150](http://dx.doi.org/10.1016/j.commatsci.2015.11.047), which can
be downloaded from
[ann.atomistic.net](http://ann.atomistic.net/download/).  Using that
data set, the training set file can be generated as described in the
first subdirectory `01-generate`.

Output files generated at each step are contained in `output`
subdirectories.

01-generate
-----------

Evaluation of the Chebyshev descriptor for all structures in the
reference data set.  The resulting feature vectors are written to a
*training set* file.

- `generate.in`: Input file for `generate.x`
- `O.fingerprint.stp` and `Ti.fingerprint.stp`: Descriptor definitions
  for the atomic species O and Ti.  For both species relatively small
  descriptor sizes with a radial expansion order of 16 and an angular
  order of 4 are used.

02-train
--------

Training examples using the *training set* files generated in the
previous step (not included because of the file size).  The results of
two training runs are shown in the subdirectories `set001` and `set002`.

- `train.in`: Input file for `train.x`
- `get-energies`: Subdirectory demonstrating how to write out the
  energies and errors of all training and testing samples (after
  training has completed).

Training will generate the ANN potential files `O.15t-15t.nn` and
`Ti.15t-15t.nn` which are also provided in the `output` subdirectory.

03-predict
----------

Usage of the ANN potentials trained in the previous step for the
prediction of the entire reference data set.

- `predict.in`: Input file for `predict.x`
- `O.15t-15t.nn` and `Ti.15t-15t.nn`: ANN potential files

The output (energies and atomic forces) can be found in the `output`
subdirectory.
