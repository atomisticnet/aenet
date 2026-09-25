#!/usr/bin/env python3
# This file is part of the AENET package.
#
# Copyright (C) 2012-2026 Nongnuch Artrith and Alexander Urban
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.

import json
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
FINALIZE = ROOT / "packaging" / "finalize"
VERSION = (ROOT / "src" / "VERSION").read_text(encoding="ascii").strip()


class MetadataTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.work = Path(self.temporary.name)
        self.compiler = self.work / "gfortran"
        self.compiler.write_text(
            "#!/bin/sh\n"
            "test \"$1\" = -dumpfullversion\n"
            "printf '14.3.0\\n'\n",
            encoding="utf-8",
        )
        self.compiler.chmod(0o755)

    def tearDown(self):
        self.temporary.cleanup()

    def stage(self, platform):
        stage = self.work / f"aenet-{VERSION}-{platform}-gnu-serial"
        for directory in ("bin", "tools", "lib", "include"):
            (stage / directory).mkdir(parents=True, exist_ok=True)
        suffixes = {
            "linux-x86_64": ("libgfortran.so.5", "libquadmath.so.0",
                             "libgcc_s.so.1"),
            "macos-arm64": ("libgfortran.5.dylib", "libquadmath.0.dylib",
                            "libgcc_s.1.1.dylib"),
        }
        for name in suffixes[platform]:
            (stage / "lib" / name).write_bytes(b"runtime\n")
        return stage

    def finalize(self, stage, *arguments, check=True):
        return subprocess.run(
            [str(FINALIZE), "--stage", str(stage), "--compiler",
             str(self.compiler), "--source-revision", "a" * 40, *arguments],
            check=check, capture_output=True, text=True,
        )

    def test_linux_metadata_and_notices(self):
        stage = self.stage("linux-x86_64")
        self.finalize(stage, "--openblas-version", "0.3.20")
        metadata = json.loads(
            (stage / "metadata.json").read_text(encoding="utf-8")
        )
        self.assertEqual(metadata["aenet_version"], VERSION)
        self.assertEqual(metadata["source_revision"], "a" * 40)
        self.assertEqual(metadata["platform"], "linux-x86_64")
        self.assertEqual(metadata["compiler"],
                         {"family": "GNU Fortran", "version": "14.3.0"})
        self.assertEqual(
            metadata["blas"],
            {"linkage": "static", "name": "OpenBLAS", "version": "0.3.20"},
        )
        self.assertEqual(metadata["file_inventory"], "manifest.txt")
        self.assertEqual(
            metadata["bundled_runtime_libraries"],
            ["libgcc_s.so.1", "libgfortran.so.5", "libquadmath.so.0"],
        )
        self.assertEqual(
            sorted(path.name for path in (stage / "licenses").iterdir()),
            ["AENET-MPL-2.0.txt", "GCC-GPL-3.0.txt",
             "GCC-RUNTIME-LIBRARY-EXCEPTION.txt", "L-BFGS-B.txt",
             "OpenBLAS.txt"],
        )

    def test_macos_metadata_uses_accelerate_without_openblas_notice(self):
        stage = self.stage("macos-arm64")
        self.finalize(stage)
        metadata = json.loads(
            (stage / "metadata.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            metadata["blas"],
            {"linkage": "system", "name": "Accelerate", "version": "system"},
        )
        self.assertFalse((stage / "licenses" / "OpenBLAS.txt").exists())

    def test_linux_requires_an_openblas_version(self):
        result = self.finalize(self.stage("linux-x86_64"), check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--openblas-version", result.stderr)


if __name__ == "__main__":
    unittest.main()
