#!/usr/bin/env python3
# This file is part of the AENET package.
#
# Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.

import hashlib
from pathlib import Path, PurePosixPath
import re


ROOT_PATTERN = re.compile(
    r"aenet-(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)-"
    r"(macos-arm64|linux-x86_64)-gnu-serial"
)
REQUIRED_DIRECTORIES = ("bin", "tools", "lib", "include", "licenses")


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
