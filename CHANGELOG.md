# Changelog

## Version 2.0.0 (2018-04-08)

### Added since v1.0.0

* Chebyshev descriptor for local atomic environments ([N. Artrith,
  A. Urban, and G. Ceder, Phys. Rev. B 96, 2017,
  014112.](http://dx.doi.org/10.1103/PhysRevB.96.014112))

### Changed since v1.0.0

* Implementation is now thread-safe so that calls to the aenet library
  can be parallelized using shared-memory parallelization.

### Fixed since v1.0.0

* Fixed a bug in the neighbor list that could result in crashes for
  periodic structures with atoms exactly on the cell boundaries.
