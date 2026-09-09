#!/bin/bash
# Experimental AENET feasibility runner; MPL-2.0, see src/license-header.txt.
set -euo pipefail
uname -a
git rev-parse HEAD
if [[ "$RUNNER_OS" == macOS ]]; then
  brew install gcc openblas
  fc="$(brew --prefix gcc)/bin/gfortran"
  args=(-DCMAKE_PREFIX_PATH="$(brew --prefix openblas)" -DCMAKE_EXE_LINKER_FLAGS=-lgomp -DCMAKE_SHARED_LINKER_FLAGS=-lgomp)
else
  sudo apt-get update -qq
  sudo apt-get install -y gfortran libopenblas-dev cmake patchelf
  fc=gfortran
  args=()
fi
"$fc" --version
cmake --version
cmake -S . -B build-l4 -DBUILD_AENET=ON -DCMAKE_Fortran_COMPILER="$fc" -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=OFF -DUSE_MKL=OFF -DUSE_OPENBLAS=ON "${args[@]}"
cmake --build build-l4 --target build_all --parallel 1
ctest --test-dir build-l4 --output-on-failure
cmake --install build-l4 --prefix "$RUNNER_TEMP/stage"
cc l4/harness.c -Isrc -L"$RUNNER_TEMP/stage/lib" -laenet -Wl,-rpath,"$RUNNER_TEMP/stage/lib" -o "$RUNNER_TEMP/stage/bin/api-smoke"
if [[ "$RUNNER_OS" == macOS ]]; then
  install_name_tool -add_rpath @loader_path/../lib "$RUNNER_TEMP/stage/bin/api-smoke"
  install_name_tool -delete_rpath "$RUNNER_TEMP/stage/lib" "$RUNNER_TEMP/stage/bin/api-smoke"
fi
python3 l4/package.py "$RUNNER_TEMP/stage"
mv "$RUNNER_TEMP/stage" "$RUNNER_TEMP/relocated"
if [[ "$RUNNER_OS" == macOS ]]; then
  otool -l "$RUNNER_TEMP/relocated/bin/predict.x_openblas" | grep -A6 LC_BUILD_VERSION
  # Block reads of non-system toolchains; this is not older-OS validation.
  profile='(version 1)(allow default)(deny file-read* (subpath "/opt/homebrew") (subpath "/usr/local") (subpath "/Library/Developer") (subpath "/Applications/Xcode.app"))'
  sandbox-exec -p "$profile" "$RUNNER_TEMP/relocated/bin/api-smoke"
  # System Python may itself need developer paths; generate inputs outside sandbox,
  # then repeat all three executable checks inside the restricted sandbox.
  python3 l4/smoke.py "$RUNNER_TEMP/relocated" "$RUNNER_TEMP/smoke"
  cd "$RUNNER_TEMP/smoke"
  rm smoke.train smoke.train.scaled Cu.nn
  for name in generate train predict; do
    args=("$name.in"); [[ "$name" != predict ]] || args+=(3.xsf)
    sandbox-exec -p "$profile" "$RUNNER_TEMP/relocated/bin/$name.x_openblas" "${args[@]}" > "$name-restricted.log" 2>&1
    if grep -E 'Error:|runtime error|Library not loaded' "$name-restricted.log"; then exit 1; fi
  done
  grep 'Total energy' predict-restricted.log
else
  docker run --rm -v "$RUNNER_TEMP/relocated:/bundle:ro" -v "$PWD/l4:/checks:ro" ubuntu:22.04 bash -ec 'apt-get update -qq; apt-get install -y python3; /bundle/bin/api-smoke; python3 /checks/smoke.py /bundle /tmp/smoke; ! command -v gfortran'
fi
printf '\nL4 runtime experiment passed\n'
