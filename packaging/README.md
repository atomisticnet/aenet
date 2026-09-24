# Binary release packaging

These commands are the maintained entry points for local release-candidate
work and release CI. They require Python 3.9 or newer. The build also requires
CMake, CTest, GNU Fortran, and the selected platform's numerical library.

`build` configures a fresh GNU serial Release build, builds all targets with
one build job, runs CTest, and installs into the release archive's top-level
directory. It selects Accelerate on macOS and OpenBLAS on Linux:

```sh
packaging/build \
  --platform macos-arm64 \
  --compiler /path/to/gfortran-14 \
  --build-dir /tmp/aenet-build \
  --stage-parent /tmp/aenet-stage
```

Use `--platform linux-x86_64` with the Linux GNU compiler for the Linux
candidate. Both the build directory and computed stage directory must be new.
This prevents an older build or staged runtime from entering a candidate.

The platform packaging steps add runtime libraries and license material to
the stage. Once it contains `bin/`, `tools/`, `lib/`, `include/`, and
`licenses/`, `package` creates the archive, `manifest.txt`, and SHA-256
sidecar:

```sh
packaging/package \
  --stage /tmp/aenet-stage/aenet-2.0.4-macos-arm64-gnu-serial \
  --output /tmp/aenet-dist \
  --epoch 1789920000
```

The epoch normalizes every archive timestamp. CI should pass the source
commit timestamp explicitly; `SOURCE_DATE_EPOCH` is used when `--epoch` is
omitted, with zero as the deterministic fallback. Archive ownership is
normalized to numeric user and group zero. File permission bits and relative
symbolic links are preserved. Absolute links, links escaping the archive,
special files, and names outside the two approved platform contracts fail
packaging.

`validate` verifies the archive name and root, SHA-256 sidecar, required
directories, safe paths and links, file modes, content hashes, and manifest:

```sh
packaging/validate \
  /tmp/aenet-dist/aenet-2.0.4-macos-arm64-gnu-serial.tar.gz
```

Platform-specific dependency inspection and independent-runtime smoke tests
extend this validation under issue 7. Run the common archive regression tests
with:

```sh
python3 packaging/tests/test_archive.py
```

## Linux x86_64

After `build` produces the Linux install tree, bundle the runtime libraries
reported by the same GNU compiler and rewrite every ELF runtime search path:

```sh
packaging/linux/prepare \
  --stage /tmp/aenet-stage/aenet-2.0.4-linux-x86_64-gnu-serial \
  --compiler /usr/bin/gfortran
```

This step requires `patchelf` and `readelf`. It copies
`libgfortran.so.5`, `libquadmath.so.0`, and `libgcc_s.so.1`, then verifies the
full shipped ELF closure. It rejects shared OpenBLAS, BLAS, or LAPACK and any
dependency outside the bundled runtimes and the Ubuntu 22.04 baseline system
libraries. Add the required license material, then create the archive with
the common `package` command.

Build the C API test before entering the independent runtime environment:

```sh
packaging/linux/build_api_smoke \
  --output /tmp/aenet-validation/api-smoke
```

The test is linked only to the system dynamic-loading API. It receives the
extracted `libaenet.so` path at runtime and checks exported API symbols,
initialization, atom-type conversion, and finalization.

The release gate runs the structural validator, extracts the actual archive,
and checks it in an Ubuntu 22.04 container:

```sh
packaging/linux/validate_container \
  /tmp/aenet-dist/aenet-2.0.4-linux-x86_64-gnu-serial.tar.gz \
  --api-smoke /tmp/aenet-validation/api-smoke
```

The container installs Python and binary-inspection tools but no compiler. It
requires an empty `LD_LIBRARY_PATH`, verifies x86_64 ELF metadata and relative
RUNPATHs, checks that bundled dependencies resolve from the extracted tree,
runs exact CLI version checks and the prebuilt C API test, and executes the
generate/train/predict numerical smoke workflow. Docker and network access to
the Ubuntu package repositories are prerequisites for this validation entry
point.

## macOS arm64

The first macOS candidate uses GNU Fortran 14 and system Accelerate. The common
build command sets `CMAKE_OSX_DEPLOYMENT_TARGET=14.0`; runtime preparation
rejects other GNU major versions:

```sh
packaging/macos/prepare \
  --stage /tmp/aenet-stage/aenet-2.0.4-macos-arm64-gnu-serial \
  --compiler /path/to/gfortran-14
```

This command recursively copies only the GNU Fortran, quadmath, and GCC
support libraries. Accelerate, libSystem, and other `/System/Library` or
`/usr/lib` dependencies remain operating-system dependencies. OpenBLAS and
OpenMP runtimes are rejected. The command rewrites bundled dependencies to
`@loader_path`, assigns relative dynamic-library IDs, verifies that every
Mach-O object is arm64 with a deployment target no newer than macOS 14.0, and
applies and verifies ad-hoc signatures after all load-command changes.

After adding the required license material, use the common `package` command.
Build the separately mounted API test and validate the final archive with:

```sh
packaging/macos/build_api_smoke \
  --output /tmp/aenet-validation/api-smoke
packaging/macos/validate_archive \
  /tmp/aenet-dist/aenet-2.0.4-macos-arm64-gnu-serial.tar.gz \
  --api-smoke /tmp/aenet-validation/api-smoke
```

Validation checks the compressed archive and a newly extracted tree. It
repeats the architecture, deployment-target, dependency, and signature checks,
then runs exact CLI version checks, the prebuilt C API test, and the numerical
generate/train/predict workflow through `sandbox-exec`. The default profile
denies reads from common Homebrew and Apple developer-toolchain locations so
the checks cannot obtain the bundled runtimes from the build installation.
