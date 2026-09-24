# Binary release feasibility and contract

Status: feasibility demonstrated on both targets; contract accepted for
implementation. This is not a released-platform support guarantee. Source
examined: `9be145096da995cbf3c68346cc8d3ebf7a0a07b3`.
No production Fortran/CMake source was changed during the experiments.

## Scope

The first release targets native macOS arm64 and Linux x86_64 GNU serial
executables, tools, native libraries, and the C header. MPI/Intel binaries,
additional architectures, and Python download automation remain deferred.

## Evidence so far

| Environment | Build/test result | Runtime result |
| --- | --- | --- |
| Local macOS 26.6.2 arm64, GNU 14.2, CMake 4.1.0-rc1, OpenBLAS 0.3.29 | All targets built with one job; 19/19 CTest entries passed | Relocated CLI and C API checks passed with toolchain-directory reads denied |
| GitHub Ubuntu 22.04 x86_64, GNU 11.4, CMake 3.31.6 | All targets built with one job; 19/19 CTest entries passed | Bundled runtime checks passed in a fresh Ubuntu 22.04 container with Python but no Fortran compiler |
| GitHub macOS 14.6 arm64, GNU 14.4, CMake 4.4.3, OpenBLAS 0.3.34, libomp 23.1 | All targets built with one job; 19/19 CTest entries passed | Relocated CLI and C API checks passed with toolchain-directory reads denied |
| GitHub macOS 14.6 arm64, GNU 14.4, Accelerate | All targets built with one job; 19/19 CTest entries passed | Relocated CLI and C API checks passed with toolchain-directory reads denied |
| Local macOS 26.6.2 arm64, GNU 14.2, Accelerate | All targets built with one job; 19/19 CTest entries passed | Relocated CLI smoke test passed; controlled OpenBLAS comparison below |

The final paired run passed on macOS but failed on Linux while downloading
an Ubuntu package index (hash/index consistency failure), before compilation.
The Linux result above comes from earlier successful runs of the same source
and Linux recipe. This infrastructure failure is not counted as a test pass.

The CTest count includes the build-test-binaries fixture. The experimental
CLI check generates eight synthetic Cu dimer structures with analytical
energies/forces, generates a training set, trains for two iterations, and
predicts twice with the resulting model. It checks output existence, explicit
error messages, finite energy, and identical repeated prediction to 1e-8 eV.
Training initializes from the system clock; energies from separate runs are
not compared. This smoke check establishes runtime functionality, not model
accuracy or scientific validation. CTest supplies the existing numerical tests.
The separately compiled C caller checks initialization, unloaded-potential
state, atom-type conversion with known expected mapping, and finalization.

## Build and dependency findings

- Static OpenBLAS does not imply a self-contained binary. Local OpenBLAS needs
  GNU libgomp; the hosted Homebrew OpenBLAS needs LLVM libomp. Use the actual
  library's link requirements rather than adding global Fortran OpenMP flags.
- Global -fopenmp activated legacy conditional code that did not compile.
  That was an unsuitable linking experiment, not an approved source change.
- Both CMake library targets write the same module directory. A parallel build
  hit an io.mod rename race; single-job builds worked. Track a focused fix
  before enabling parallel CI builds.
- Local macOS direct GCC dependencies were libgfortran, libquadmath, and
  libgomp. libgcc_s was a transitive @rpath dependency missed by the first
  copy experiment; denying Homebrew reads exposed it. Including it made the
  restricted API and executable tests pass.
- Hosted macOS GNU 16.2 built and passed CTest, but the C harness failed to
  link aenet_init and aenet_convert_atom_types. GNU 14.4 passed the same
  harness. Pin GNU 14 for the first macOS candidate; the missing exports
  need separate investigation before claiming GNU 16 library compatibility.
- The hosted OpenBLAS variant bundles LLVM libomp (not GNU libgomp), plus libgfortran,
  libquadmath, and libgcc_s. Its reported non-system dependency closure
  resolves within the relocated tree; libSystem remains external.
- Linux experiments bundled libgfortran.so.5, libquadmath.so.0, and libgcc_s.so.1
  while retaining system libc/libm/loader. OpenBLAS was statically linked.
- Developer build trees retain variant suffixes on the main executables.
  Installed trees expose stable executable names and include `aenet.h`.

## Release choices

- Require Ubuntu 22.04/glibc 2.35 or newer for the first Linux release. This
  is the demonstrated runtime baseline; do not claim generic Linux or
  older-glibc compatibility yet.
- Require macOS 14 or newer on arm64 for the first macOS release and use GNU
  14 to build it. The hosted
  experiment ran on macOS 14.6 (Darwin 23.6), and predict reports minos 14.0
  and SDK 14.0. This does not establish every macOS 14 patch level or every
  bundled library deployment target; inspect the full closure and validate
  the advertised minimum before release. Local GNU 14.2 artifacts with minos
  16.0 are not the proposed release artifacts.
- Prefer Accelerate for the macOS candidate on the demonstrated correctness
  and packaging evidence: it removes the OpenBLAS/OpenMP runtime dependency.
  This is not a measured speedup recommendation. Keep OpenBLAS for Linux.
  Bundle non-system compiler runtimes on both platforms; keep system
  libraries/frameworks external.
- Use relative runtime load paths; remove build-directory rpaths and inspect
  transitive dependency closure. Apply ad-hoc signatures after modifying
  Mach-O files where required. This is not Developer ID signing/notarization.
- Release archive names are
  `aenet-<version>-macos-arm64-gnu-serial.tar.gz` and
  `aenet-<version>-linux-x86_64-gnu-serial.tar.gz`, with `bin/`, `tools/`,
  `lib/`, `include/aenet.h`, and notices. Issue 7 adds the runtime libraries
  and notices without changing this top-level layout.

## Accelerate extension: correctness, packaging, and performance

Accelerate selection requires no production source change: configure a fresh
build with USE_OPENBLAS=OFF, USE_MKL=OFF, and BLA_VENDOR=Apple. CMake finds
both BLAS and LAPACK in Accelerate.framework. Remove the experimental OpenMP
linker flags used for OpenBLAS. The existing BLA_STATIC setting does not make
the system framework static; otool confirms the system Accelerate dependency.

The hosted Accelerate build reports minos 14.0/SDK 14.0 and passes the same
CTest and relocated API/CLI checks as OpenBLAS. Its bundled dependency closure
contains libgfortran, libquadmath, and libgcc_s; no libomp or libgomp is needed.
The Apple framework remains an OS dependency. Signing/notarization and final
minimum-patch-level validation are still separate release gates.

The local binaries reference conventional Fortran symbols such as dgemm_,
dgemv_, and dpotrf_. This experiment did not select Apple's newer LAPACK
symbol interface or introduce wrappers/aliases. Performance conclusions apply
only to the interface tested.

### Controlled local timing comparison

Both variants used the same source, local GNU 14.2, Release optimization,
and -fexternal-blas. The OpenBLAS comparator was 0.3.29, not the hosted 0.3.34
build. The synthetic workload contains 64 isolated 27-atom Cu grids with
varying spacing, a 22-32-32-1 network, 10 BFGS training iterations, and 256
structure evaluations with forces for prediction. These exercise realistic
code paths and moderately sized networks, but are not production datasets.
Training uses synthetic energy labels; placeholder force labels are not used
for fitting. No model-accuracy claim is made.

Each backend receives an identical fixed starting network and data within
one thread-limit comparison. Each timed process starts from a fresh copy;
prediction uses the fixed seed model, not the newly trained outputs. One
warm-up is excluded and three measured runs are alternated between backends.
Timings are end-to-end process wall times including input/output.
OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, and VECLIB_MAXIMUM_THREADS are all set
to the indicated limit. This controls requested library limits, not proof of
actual simultaneous worker counts; it does not enable OpenMP in AENET.

| Requested thread limit | Operation | OpenBLAS median (s) | Accelerate median (s) |
| --- | --- | --- | --- |
| 1 | train | 3.527 | 3.696 |
| 1 | predict | 0.498 | 0.518 |
| 4 | train | 3.491 | 3.947 |
| 4 | predict | 0.491 | 0.517 |

Single-thread measured ranges were 3.425-3.592 s versus 3.540-3.949 s for
training and 0.474-0.535 s versus 0.490-0.553 s for prediction
(OpenBLAS versus Accelerate). Timing variation and the limited workload do
not support a general performance ranking. No Accelerate speed advantage
was demonstrated here; larger networks or different LAPACK-heavy training
methods may behave differently and need representative user workloads.

Fixed-model prediction energies agreed at printed precision (eight decimal
places) in both thread-limit comparisons. In the single-thread comparison,
all 20,736 printed force components also agreed (six decimal places).
This establishes agreement at output precision, not bitwise internal equality.
No tolerance changes were made. Training trajectories from different random
initializations are not used as a numerical equivalence check.

The benchmark script, raw timings, and seed model are retained locally under
dev-notes/l4/accelerate/ and the temporary benchmark directories. Seed-model
hash and numeric results are recorded with the benchmark output. Before
claiming a production performance benefit, repeat with representative data
and the intended training method; such a claim is not needed to choose the
simpler macOS runtime packaging.

## Notices, metadata, and release gates

Production packaging installs the AENET MPL-2.0 license, the license supplied
with the bundled L-BFGS-B 3.0 sources, GPLv3 and the GCC Runtime Library
Exception for the bundled GNU runtimes, and the OpenBLAS BSD notice on Linux.
macOS uses system Accelerate and therefore has no OpenBLAS or OpenMP payload or
notice. `metadata.json` records the source revision, AENET and compiler
versions, target platform, build configuration, BLAS/LAPACK provider and
linkage, bundled runtime filenames, installed notices, and the location of the
hashed file inventory. The archive validator enforces this schema and the
platform-specific notice set.

Release gates must cover clean builds/CTest, extracted-file and architecture
checks, transitive dependency resolution, independent-runtime CLI/API tests,
checksums/notices, and issue 8's install instructions. Test actual downloaded
macOS quarantine/security behavior during documentation/publication; local
sandbox execution does not establish it. Choose any Developer ID/notarization
requirements from that evidence rather than claiming they are resolved.

## Reproduction and evidence

The experiment scripts preserve exact commands and failure logs. Essential
configuration for the original OpenBLAS experiment: BUILD_AENET=ON, Release,
USE_MPI=OFF, USE_MKL=OFF, USE_OPENBLAS=ON, and build_all with one build job,
followed by CTest/install. The Accelerate overrides are described above.
On macOS the successful compiler is Homebrew gcc@14's gfortran-14; the linker
flags add the Homebrew libomp library directory and -lomp to executable and
shared-library links. OpenBLAS's prefix is supplied to CMake. Linux uses the
Ubuntu gfortran/libopenblas-dev packages without those macOS linker flags.

After staging, the experiment copies non-system runtime dependencies,
rewrites relative load paths, and moves the staging tree to a new directory.
Linux then runs in an Ubuntu 22.04 container containing Python but no compiler.
macOS runs the C caller and each CLI through sandbox-exec with reads denied
under /opt/homebrew, /usr/local, /Library/Developer, and /Applications/Xcode.app.
The inspected dependency closure is bundled libraries plus system libSystem
and, for the Accelerate variant, the system Accelerate framework. This is evidence of toolchain-independent runtime resolution,
not a pristine macOS installation or a downloaded/quarantined archive test.

Immutable successful evidence:

- [Linux log](https://github.com/atomisticnet/aenet/blob/b0b1db8b2c0d34d367643f41707afe9e3ebdcd17/evidence/ubuntu-22.04.log).
- [macOS log](https://github.com/atomisticnet/aenet/blob/eb61306016163e83f2f82c8bd74174ac0e6ab8c9/evidence/macos-14.log).
- [OpenBLAS experimental scripts](https://github.com/atomisticnet/aenet/tree/e085f99/l4).
- [Accelerate macOS log](https://github.com/atomisticnet/aenet/blob/eb771a152db2aae28a4fdb1a4079d8b34415eefa/evidence/macos-14-accelerate.log).
- [Accelerate experimental scripts](https://github.com/atomisticnet/aenet/tree/53e9fd9/l4).

## Follow-up

The temporary experiment branches were deleted after their durable evidence
was recorded at the immutable links above. Their commits and workflows were
experimental and were not merged. Private raw local logs remain in
dev-notes/l4/.

GitHub currently provides native macOS arm64 and Ubuntu x86_64 runners:
[hosted runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).
Runner labels establish execution availability, not artifact compatibility.

Both platform feasibility experiments and production archive implementations
have successful native evidence. Linux production validation is recorded in
[run 35815467080](https://github.com/atomisticnet/aenet/actions/runs/35815467080),
and macOS arm64 production validation is recorded in
[run 36026156581](https://github.com/atomisticnet/aenet/actions/runs/36026156581).
The completed notice and metadata contract passed for both native candidates
from one revision in [paired run
36027659670](https://github.com/atomisticnet/aenet/actions/runs/36027659670).
The maintained commands under `packaging/` now own build, runtime relocation,
metadata/notices, deterministic archive creation, and extracted-archive
validation. Issue 4 consumes these commands for CI orchestration; issue 8 owns
user instructions, and local publication work verifies actual downloads and
macOS quarantine behavior. No candidate archive has been published.
