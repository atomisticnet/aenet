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
