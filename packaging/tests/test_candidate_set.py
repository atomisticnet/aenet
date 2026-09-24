#!/usr/bin/env python3
# This file is part of the AENET package.
#
# Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban
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
PACKAGE = ROOT / "packaging" / "package"
VERIFY = ROOT / "packaging" / "verify_candidates"
VERSION = (ROOT / "src" / "VERSION").read_text(encoding="ascii").strip()
REVISION = "a" * 40
COMMON_NOTICES = [
    "AENET-MPL-2.0.txt",
    "GCC-GPL-3.0.txt",
    "GCC-RUNTIME-LIBRARY-EXCEPTION.txt",
    "L-BFGS-B.txt",
]


class CandidateSetTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.work = Path(self.temporary.name)
        self.candidates = self.work / "candidates"
        self.create_candidate("linux-x86_64", REVISION)
        self.create_candidate("macos-arm64", REVISION)

    def tearDown(self):
        self.temporary.cleanup()

    def create_candidate(self, platform, revision):
        root = f"aenet-{VERSION}-{platform}-gnu-serial"
        stage = self.work / root
        for directory in ("bin", "tools", "lib", "include", "licenses"):
            (stage / directory).mkdir(parents=True, exist_ok=True)
        program = stage / "bin" / "generate.x"
        program.write_text("#!/bin/sh\n", encoding="utf-8")
        program.chmod(0o755)
        notices = list(COMMON_NOTICES)
        if platform == "linux-x86_64":
            notices.append("OpenBLAS.txt")
            blas = {"name": "OpenBLAS", "version": "0.3.20",
                    "linkage": "static"}
            runtime = "libgfortran.so.5"
        else:
            blas = {"name": "Accelerate", "version": "system",
                    "linkage": "system"}
            runtime = "libgfortran.5.dylib"
        (stage / "lib" / runtime).write_bytes(b"runtime\n")
        for name in notices:
            (stage / "licenses" / name).write_text(
                "notice\n", encoding="utf-8"
            )
        metadata = {
            "schema_version": 1,
            "aenet_version": VERSION,
            "archive_root": root,
            "platform": platform,
            "build_configuration": "Release GNU serial",
            "source_revision": revision,
            "compiler": {"family": "GNU Fortran", "version": "14.3.0"},
            "blas": blas,
            "bundled_runtime_libraries": [runtime],
            "licenses": sorted(notices),
            "file_inventory": "manifest.txt",
        }
        (stage / "metadata.json").write_text(
            json.dumps(metadata) + "\n", encoding="utf-8"
        )
        subprocess.run(
            [PACKAGE, "--stage", stage, "--output", self.candidates,
             "--epoch", "1234567890"], check=True,
        )

    def run_verify(self, *arguments, check=True):
        return subprocess.run(
            [VERIFY, "--directory", self.candidates,
             "--source-revision", REVISION, *arguments],
            check=check, capture_output=True, text=True,
        )

    def test_accepts_complete_matching_candidate_set(self):
        result = self.run_verify("--tag", f"v{VERSION}")
        self.assertIn(f"validated candidate set for {VERSION}", result.stdout)

    def test_rejects_missing_platform_candidate(self):
        macos = self.candidates / (
            f"aenet-{VERSION}-macos-arm64-gnu-serial.tar.gz"
        )
        macos.unlink()
        Path(str(macos) + ".sha256").unlink()
        result = self.run_verify(check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("candidate files", result.stderr)

    def test_rejects_mismatched_source_revision(self):
        for path in self.candidates.iterdir():
            path.unlink()
        self.create_candidate("linux-x86_64", REVISION)
        self.create_candidate("macos-arm64", "b" * 40)
        result = self.run_verify(check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("source revision", result.stderr)

    def test_rejects_tag_that_does_not_match_version(self):
        result = self.run_verify("--tag", "v9.9.9", check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("tag does not match", result.stderr)


if __name__ == "__main__":
    unittest.main()
