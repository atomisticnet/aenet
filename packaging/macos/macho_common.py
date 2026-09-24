#!/usr/bin/env python3
# This file is part of the AENET package.
#
# Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.

import os
from pathlib import Path
import re
import subprocess


SYSTEM_PREFIXES = ("/System/Library/", "/usr/lib/")
GCC_RUNTIME_PATTERNS = (
    re.compile(r"libgfortran(?:\.[0-9]+)*\.dylib$"),
    re.compile(r"libquadmath(?:\.[0-9]+)*\.dylib$"),
    re.compile(r"libgcc_s(?:\.[0-9]+)*\.dylib$"),
)
FORBIDDEN_RUNTIME_NAMES = ("openblas", "libomp", "libgomp")


def command_output(command):
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    return subprocess.check_output(
        [str(item) for item in command], text=True, stderr=subprocess.STDOUT,
        env=environment,
    )


def parse_dependencies(output):
    dependencies = []
    for line in output.splitlines()[1:]:
        line = line.strip()
        if line:
            dependencies.append(line.split(" (", 1)[0])
    return dependencies


def parse_architectures(output):
    architectures = set(output.split())
    if architectures != {"arm64"}:
        raise ValueError("release object must contain only arm64 architecture")
    return architectures


def version_tuple(value):
    parts = value.split(".")
    return tuple(int(part) for part in parts[:2]) + (0,) * (2 - len(parts))


def parse_build_versions(output):
    versions = []
    legacy_block = False
    for line in output.splitlines():
        stripped = line.strip()
        if stripped == "cmd LC_VERSION_MIN_MACOSX":
            legacy_block = True
        elif stripped.startswith("cmd "):
            legacy_block = False
        match = re.match(r"minos\s+(\d+(?:\.\d+)*)$", stripped)
        if match:
            versions.append(version_tuple(match.group(1)))
        elif legacy_block:
            match = re.match(r"version\s+(\d+(?:\.\d+)*)$", stripped)
            if match:
                versions.append(version_tuple(match.group(1)))
                legacy_block = False
    if not versions:
        raise ValueError("Mach-O deployment target is missing")
    if any(version > (14, 0) for version in versions):
        raise ValueError("Mach-O deployment target is newer than macOS 14.0")
    return versions


def is_gcc_runtime(name):
    return any(pattern.fullmatch(name) for pattern in GCC_RUNTIME_PATTERNS)


def check_dependencies(dependencies, bundled):
    for dependency in dependencies:
        name = Path(dependency).name
        lowered = name.lower()
        if any(token in lowered for token in FORBIDDEN_RUNTIME_NAMES):
            raise ValueError(f"OpenBLAS/OpenMP dependency is not allowed: {name}")
        if dependency.startswith(SYSTEM_PREFIXES):
            continue
        if dependency.startswith("@loader_path/") and name in bundled:
            continue
        if dependency.startswith("/"):
            raise ValueError(f"toolchain path remains in dependency: {dependency}")
        raise ValueError(f"unresolved Mach-O dependency: {dependency}")


def verify_adhoc_signature(path, codesign="codesign"):
    result = subprocess.run(
        [codesign, "-d", "--verbose=4", str(path)], text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    if result.returncode != 0 or "Signature=adhoc" not in result.stdout:
        raise ValueError(f"Mach-O object lacks an ad-hoc signature: {path}")


def macho_files(root, lipo="lipo", include_archives=False):
    files = []
    for directory in ("bin", "tools", "lib"):
        for path in sorted((root / directory).iterdir()):
            if path.is_symlink() or not path.is_file():
                continue
            if path.suffix == ".a" and not include_archives:
                continue
            result = subprocess.run(
                [lipo, "-archs", str(path)], text=True,
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            )
            if result.returncode == 0:
                parse_architectures(result.stdout)
                files.append(path)
    return files


def dependencies(path, otool="otool"):
    return parse_dependencies(command_output([otool, "-L", path]))


def dylib_id(path, otool="otool"):
    result = subprocess.run(
        [otool, "-D", str(path)], text=True,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    lines = result.stdout.splitlines()
    return lines[1].strip() if result.returncode == 0 and len(lines) > 1 else None
