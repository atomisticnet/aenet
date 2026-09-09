---
name: build-test
description: Configure, build, and validate the aenet Fortran backend with CMake and CTest, including compiler, MPI, and BLAS variants relevant to the change.
---

# Build and Test

Read root and src/CMakeLists.txt and any config/local.toml before choosing
commands. CMake is authoritative for options and targets. Local TOML records
tool locations for the agent; CMake does not read it. Do not eval local
configuration as shell code; construct commands with properly quoted arguments.

## Select a configuration

Use an out-of-source directory per compiler, build type, MPI, and BLAS
combination. Do not reuse an unrelated existing build cache or change its
compiler. Preserve source-tree artifacts and other work already present.

For a GNU serial Debug baseline, with tools available on PATH:

```sh
cmake -S . -B build-gnu-debug -DBUILD_AENET=ON \
  -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Debug \
  -DUSE_MPI=OFF -DUSE_MKL=OFF -DUSE_OPENBLAS=OFF
cmake --build build-gnu-debug --target build_tests
ctest --test-dir build-gnu-debug -N
ctest --test-dir build-gnu-debug --output-on-failure
```

Choose a new directory if this example name already belongs to other work.
Use `ctest -R '<pattern>'` for focused tests after checking registration with
`-N`; ensure the pattern actually selects the intended cases. The registered
build_test_binaries fixture can also build tests automatically. CTest's
`--test-dir` needs CTest 3.20 or newer; with older supported CMake versions,
run CTest from inside the chosen build directory instead.

The default build excludes the main artifacts. Use explicit targets:
`main`, `lib`, `tools`, `build_tests`, or `build_all`. New Fortran tests belong
under src/tests/ and must be registered in src/CMakeLists.txt, including the
existing test-build target and fixture dependencies as appropriate.

For a test-first change, add or extend a focused case using src/ext/unittest.f90
and the existing test-driver pattern. Build it and confirm CTest reports the
intended failure before implementing the fix. Ensure assertion failures reach
a nonzero process exit status (see tst_exit_nonzero_if_failed); printed failure
messages alone are not sufficient. Then implement, rebuild, rerun the focused
case, and broaden validation according to risk. CTest runs the test programs;
it does not require a different development cycle or a new test framework.

## Validate according to the change

- Numerical changes: focused module tests, then the broader suite for shared
  routines. Debug checks help detect bounds and initialization errors; test
  Release too when optimization or floating-point behavior is relevant.
- Compiler/build changes: configure and build the affected GNU/Intel variants
  when available, using CMAKE_Fortran_COMPILER and separate directories.
- MPI changes: build with USE_MPI=ON and run relevant serial/multi-rank
  comparisons with bounded execution. Current CTest registrations do not
  launch mpiexec; an MPI build alone is not multi-rank validation.
- BLAS changes: exercise the affected system/OpenBLAS/MKL configuration;
  choose only one vendor option and verify actual resolved linkage.
- CLI/library/install changes: build the affected artifacts, exercise their
  public boundary, and test installation into a temporary prefix when needed.
  Inspect actual output names; build configuration currently adds suffixes.
- Policy-only changes: validate links, examples against configuration, and
  skill structure. A full Fortran rebuild is not required for prose alone.

Start with serial test execution. Enable CTest parallelism only when selected
cases have independent output paths. Do not change tests or numerical
thresholds merely to obtain a pass. Report baseline failures separately from
failures caused by the change.

Record exact commands, compiler/version, configuration, results, and unrun
required checks. An unavailable Intel compiler, MPI runtime, or vendor
library limits coverage; it does not demonstrate portability. Resolve or
record the gap according to the task's acceptance criteria before closure.
