#!/usr/bin/env python3
# This file is part of the AENET package.
#
# Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.

from pathlib import Path
import sys
import unittest


PACKAGING = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGING / "linux"))

from elf_common import (check_dependency_policy, check_runtime_resolution,
                        parse_dynamic, parse_ldd, parse_machine)  # noqa: E402


class ElfPolicyTests(unittest.TestCase):
    def test_parses_needed_sonames_and_runpath(self):
        dynamic = """
 0x0000000000000001 (NEEDED) Shared library: [libgfortran.so.5]
 0x0000000000000001 (NEEDED) Shared library: [libm.so.6]
 0x000000000000001d (RUNPATH) Library runpath: [$ORIGIN/../lib]
"""
        needed, runpath = parse_dynamic(dynamic)
        self.assertEqual(needed, {"libgfortran.so.5", "libm.so.6"})
        self.assertEqual(runpath, "$ORIGIN/../lib")

    def test_rejects_dynamic_blas_and_unapproved_dependencies(self):
        with self.assertRaisesRegex(ValueError, "dynamic BLAS"):
            check_dependency_policy({"libopenblas.so.0"}, set())
        with self.assertRaisesRegex(ValueError, "unresolved dependency"):
            check_dependency_policy({"libunexpected.so.1"}, set())

    def test_accepts_bundled_gcc_runtimes_and_baseline_system_libraries(self):
        needed = {"libgfortran.so.5", "libm.so.6", "libc.so.6"}
        check_dependency_policy(needed, {"libgfortran.so.5"})

    def test_requires_linux_x86_64_elf_objects(self):
        header = "  Class:                             ELF64\n" \
                 "  Machine:                           Advanced Micro Devices X86-64\n"
        self.assertEqual(parse_machine(header), "x86_64")
        with self.assertRaisesRegex(ValueError, "x86_64"):
            parse_machine("Class: ELF64\nMachine: AArch64\n")

    def test_runtime_dependencies_must_resolve_inside_the_archive(self):
        output = """
 libgfortran.so.5 => /bundle/lib/libgfortran.so.5 (0x1234)
 libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x5678)
"""
        resolved = parse_ldd(output)
        check_runtime_resolution(
            resolved, {"libgfortran.so.5"}, Path("/bundle/lib")
        )
        with self.assertRaisesRegex(ValueError, "outside the archive"):
            check_runtime_resolution(
                resolved, {"libgfortran.so.5"}, Path("/other/lib")
            )
        with self.assertRaisesRegex(ValueError, "missing from ldd"):
            check_runtime_resolution(
                {}, {"libgfortran.so.5"}, Path("/bundle/lib")
            )

    def test_rejects_missing_ldd_dependencies(self):
        with self.assertRaisesRegex(ValueError, "not found"):
            parse_ldd("libgfortran.so.5 => not found\n")


if __name__ == "__main__":
    unittest.main()
