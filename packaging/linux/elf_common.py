#!/usr/bin/env python3
# This file is part of the AENET package.
#
# Copyright (C) 2012-2026 Nongnuch Artrith and Alexander Urban
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.

from pathlib import Path
import os
import re
import subprocess


BUNDLED_RUNTIMES = ("libgfortran.so.5", "libquadmath.so.0", "libgcc_s.so.1")
SYSTEM_LIBRARIES = {
    "ld-linux-x86-64.so.2",
    "libc.so.6",
    "libdl.so.2",
    "libm.so.6",
    "libpthread.so.0",
    "librt.so.1",
}
DYNAMIC_BLAS_PREFIXES = ("libopenblas.so", "libblas.so", "liblapack.so")


def command_output(command):
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    return subprocess.check_output(
        [str(item) for item in command], text=True, stderr=subprocess.STDOUT,
        env=environment,
    )


def parse_dynamic(output):
    needed = set(re.findall(r"\(NEEDED\).*?\[([^]]+)\]", output))
    paths = re.findall(r"\((?:RPATH|RUNPATH)\).*?\[([^]]*)\]", output)
    if len(paths) > 1:
        raise ValueError("ELF object contains multiple runtime search paths")
    return needed, paths[0] if paths else None


def parse_machine(output):
    elf_class = re.search(r"^\s*Class:\s*(\S+)", output, re.MULTILINE)
    machine = re.search(r"^\s*Machine:\s*(.+?)\s*$", output, re.MULTILINE)
    if (elf_class is None or elf_class.group(1) != "ELF64" or
            machine is None or machine.group(1) !=
            "Advanced Micro Devices X86-64"):
        raise ValueError("release object is not Linux x86_64 ELF64")
    return "x86_64"


def check_dependency_policy(needed, bundled):
    dynamic_blas = sorted(
        name for name in needed if name.startswith(DYNAMIC_BLAS_PREFIXES)
    )
    if dynamic_blas:
        raise ValueError("dynamic BLAS/LAPACK dependency: " +
                         ", ".join(dynamic_blas))
    unresolved = sorted(needed - bundled - SYSTEM_LIBRARIES)
    if unresolved:
        raise ValueError("unresolved dependency: " + ", ".join(unresolved))


def parse_ldd(output):
    resolved = {}
    for line in output.splitlines():
        match = re.match(r"\s*(\S+)\s+=>\s+(\S+)", line)
        if match:
            resolved[match.group(1)] = match.group(2)
    if "not found" in output:
        raise ValueError("runtime dependency is not found")
    return resolved


def check_runtime_resolution(resolved, bundled, library_directory):
    expected_directory = library_directory.resolve()
    for name in sorted(bundled):
        if name not in resolved:
            raise ValueError(f"bundled dependency missing from ldd: {name}")
        path = Path(resolved[name])
        if not path.is_absolute() or path.resolve().parent != expected_directory:
            raise ValueError(f"{name} resolved outside the archive: {path}")


def elf_files(root, readelf="readelf"):
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    files = []
    for directory in ("bin", "tools", "lib"):
        for path in sorted((root / directory).iterdir()):
            if path.is_symlink() or not path.is_file() or path.suffix == ".a":
                continue
            result = subprocess.run(
                [readelf, "-h", str(path)], text=True,
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                env=environment,
            )
            if result.returncode == 0:
                parse_machine(result.stdout)
                files.append(path)
    return files


def inspect_elf(path, readelf="readelf"):
    parse_machine(command_output([readelf, "-h", path]))
    return parse_dynamic(command_output([readelf, "-d", path]))
