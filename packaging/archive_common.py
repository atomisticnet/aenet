#!/usr/bin/env python3
# This file is part of the AENET package.
#
# Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.

import hashlib
import json
from pathlib import Path, PurePosixPath
import re


ROOT_PATTERN = re.compile(
    r"aenet-(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)-"
    r"(macos-arm64|linux-x86_64)-gnu-serial"
)
REQUIRED_DIRECTORIES = ("bin", "tools", "lib", "include", "licenses")
COMMON_NOTICES = {
    "AENET-MPL-2.0.txt",
    "GCC-GPL-3.0.txt",
    "GCC-RUNTIME-LIBRARY-EXCEPTION.txt",
    "L-BFGS-B.txt",
}


def check_root_name(name):
    if ROOT_PATTERN.fullmatch(name) is None:
        raise ValueError(
            "archive root must be aenet-<version>-macos-arm64-gnu-serial "
            "or aenet-<version>-linux-x86_64-gnu-serial"
        )


def check_member_path(name, root):
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"unsafe archive path: {name}")
    if not path.parts or path.parts[0] != root:
        raise ValueError(f"archive member is outside {root}: {name}")


def check_link_target(member_path, target, root):
    target_path = PurePosixPath(target)
    if target_path.is_absolute():
        raise ValueError(f"unsafe symbolic link: {member_path} -> {target}")
    destination = PurePosixPath(member_path).parent.joinpath(target_path)
    stack = []
    for part in destination.parts:
        if part in ("", "."):
            continue
        if part == "..":
            if not stack:
                raise ValueError(
                    f"unsafe symbolic link: {member_path} -> {target}"
                )
            stack.pop()
        else:
            stack.append(part)
    if not stack or stack[0] != root:
        raise ValueError(f"unsafe symbolic link: {member_path} -> {target}")


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def required_notices(platform):
    notices = set(COMMON_NOTICES)
    if platform == "linux-x86_64":
        notices.add("OpenBLAS.txt")
    return notices


def validate_metadata(metadata, root, license_names):
    match = ROOT_PATTERN.fullmatch(root)
    if match is None:
        check_root_name(root)
    version = ".".join(match.group(index) for index in (1, 2, 3))
    platform = match.group(4)
    expected_keys = {
        "schema_version", "aenet_version", "archive_root", "platform",
        "build_configuration", "source_revision", "compiler", "blas",
        "bundled_runtime_libraries", "licenses", "file_inventory",
    }
    if not isinstance(metadata, dict) or set(metadata) != expected_keys:
        raise ValueError("metadata fields do not match schema version 1")
    if metadata["schema_version"] != 1:
        raise ValueError("unsupported metadata schema version")
    if metadata["aenet_version"] != version or metadata["archive_root"] != root:
        raise ValueError("metadata version or archive root is inconsistent")
    if metadata["platform"] != platform:
        raise ValueError("metadata platform is inconsistent")
    if metadata["build_configuration"] != "Release GNU serial":
        raise ValueError("metadata build configuration is inconsistent")
    revision = metadata["source_revision"]
    if (not isinstance(revision, str) or
            re.fullmatch(r"[0-9a-f]{40}", revision) is None):
        raise ValueError("metadata source revision is invalid")
    compiler = metadata["compiler"]
    if not isinstance(compiler, dict) or set(compiler) != {"family", "version"}:
        raise ValueError("metadata compiler fields are invalid")
    if (compiler["family"] != "GNU Fortran" or
            not isinstance(compiler["version"], str) or
            not compiler["version"]):
        raise ValueError("metadata compiler is invalid")
    blas = metadata["blas"]
    if not isinstance(blas, dict) or set(blas) != {"name", "version", "linkage"}:
        raise ValueError("metadata BLAS fields are invalid")
    expected_blas = ("OpenBLAS", "static") if platform == "linux-x86_64" \
        else ("Accelerate", "system")
    if ((blas["name"], blas["linkage"]) != expected_blas or
            not isinstance(blas["version"], str) or not blas["version"]):
        raise ValueError("metadata BLAS dependency is inconsistent")
    runtimes = metadata["bundled_runtime_libraries"]
    if (not isinstance(runtimes, list) or not runtimes or
            not all(isinstance(name, str) and name for name in runtimes) or
            runtimes != sorted(set(runtimes))):
        raise ValueError("metadata runtime inventory is invalid")
    expected_notices = required_notices(platform)
    if set(license_names) != expected_notices:
        raise ValueError("archive license files are incomplete or unexpected")
    if (not isinstance(metadata["licenses"], list) or
            metadata["licenses"] != sorted(expected_notices)):
        raise ValueError("metadata license inventory is inconsistent")
    if metadata["file_inventory"] != "manifest.txt":
        raise ValueError("metadata file inventory reference is invalid")


def read_metadata(path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, UnicodeError) as error:
        raise ValueError(f"invalid metadata.json: {error}") from error
