#!/usr/bin/env python3
# This file is part of the AENET package.
#
# Copyright (C) 2012-2026 Nongnuch Artrith and Alexander Urban
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.

from pathlib import Path
import sys
import unittest


PACKAGING = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGING / "macos"))

from macho_common import (check_dependencies, parse_build_versions,  # noqa: E402
                          parse_dependencies, parse_architectures)


class MachOPolicyTests(unittest.TestCase):
    def test_parses_otool_dependencies(self):
        output = (
            "/tmp/predict.x:\n"
            "\t/opt/homebrew/lib/gcc/14/libgfortran.5.dylib "
            "(compatibility version 6.0.0)\n"
            "\t/System/Library/Frameworks/Accelerate.framework/Versions/A/"
            "Accelerate (compatibility version 1.0.0)\n"
        )
        self.assertEqual(
            parse_dependencies(output),
            ["/opt/homebrew/lib/gcc/14/libgfortran.5.dylib",
             "/System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate"],
        )

    def test_accepts_only_arm64(self):
        self.assertEqual(parse_architectures("arm64\n"), {"arm64"})
        with self.assertRaisesRegex(ValueError, "arm64"):
            parse_architectures("x86_64 arm64\n")

    def test_requires_deployment_target_no_newer_than_macos_15(self):
        output = """      cmd LC_BUILD_VERSION
    minos 15.0
      sdk 15.4
      cmd LC_BUILD_VERSION
    minos 14.6
      sdk 15.4
"""
        self.assertEqual(parse_build_versions(output), [(15, 0), (14, 6)])
        with self.assertRaisesRegex(ValueError, "deployment target"):
            parse_build_versions("minos 15.1\n")

    def test_rejects_toolchain_and_openblas_dependencies(self):
        with self.assertRaisesRegex(ValueError, "toolchain path"):
            check_dependencies(["/opt/homebrew/lib/libgfortran.5.dylib"], set())
        with self.assertRaisesRegex(ValueError, "OpenBLAS"):
            check_dependencies(["@loader_path/libopenblas.dylib"],
                               {"libopenblas.dylib"})
        with self.assertRaisesRegex(ValueError, "unresolved"):
            check_dependencies(["@rpath/libgfortran.5.dylib"],
                               {"libgfortran.5.dylib"})

    def test_accepts_relative_gcc_runtimes_and_system_accelerate(self):
        dependencies = [
            "@loader_path/../lib/libgfortran.5.dylib",
            "/System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate",
            "/usr/lib/libSystem.B.dylib",
        ]
        check_dependencies(dependencies, {"libgfortran.5.dylib"})


if __name__ == "__main__":
    unittest.main()
