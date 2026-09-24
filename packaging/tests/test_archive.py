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
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "packaging" / "build"
PACKAGE = ROOT / "packaging" / "package"
VALIDATE = ROOT / "packaging" / "validate"
VERSION = (ROOT / "src" / "VERSION").read_text(encoding="ascii").strip()
ARCHIVE_ROOT = f"aenet-{VERSION}-linux-x86_64-gnu-serial"


class ArchiveTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.work = Path(self.temporary.name)
        self.stage = self.work / ARCHIVE_ROOT
        (self.stage / "bin").mkdir(parents=True)
        (self.stage / "lib").mkdir()
        (self.stage / "tools").mkdir()
        (self.stage / "include").mkdir()
        (self.stage / "licenses").mkdir()
        program = self.stage / "bin" / "generate.x"
        program.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        program.chmod(0o755)
        library = self.stage / "lib" / "libaenet.so.2.0.4"
        library.write_bytes(b"library\n")
        os.symlink(library.name, self.stage / "lib" / "libaenet.so.2")
        notices = ["AENET-MPL-2.0.txt", "GCC-GPL-3.0.txt",
                   "GCC-RUNTIME-LIBRARY-EXCEPTION.txt", "L-BFGS-B.txt",
                   "OpenBLAS.txt"]
        for name in notices:
            (self.stage / "licenses" / name).write_text(
                "license\n", encoding="utf-8"
            )
        metadata = {
            "schema_version": 1,
            "aenet_version": VERSION,
            "archive_root": ARCHIVE_ROOT,
            "platform": "linux-x86_64",
            "build_configuration": "Release GNU serial",
            "source_revision": "a" * 40,
            "compiler": {"family": "GNU Fortran", "version": "11.4.0"},
            "blas": {"name": "OpenBLAS", "version": "0.3.20",
                     "linkage": "static"},
            "bundled_runtime_libraries": ["libgfortran.so.5"],
            "licenses": notices,
            "file_inventory": "manifest.txt",
        }
        (self.stage / "metadata.json").write_text(
            json.dumps(metadata) + "\n", encoding="utf-8"
        )

    def tearDown(self):
        self.temporary.cleanup()

    def run_package(self, output, stage=None, check=True):
        return subprocess.run(
            [str(PACKAGE), "--stage", str(stage or self.stage),
             "--output", str(output), "--epoch", "1234567890"],
            check=check, capture_output=True, text=True,
        )

    def test_archive_is_reproducible_and_self_describing(self):
        first = self.work / "first"
        second = self.work / "second"
        self.run_package(first)
        self.run_package(second)
        archive_name = ARCHIVE_ROOT + ".tar.gz"
        archive1 = first / archive_name
        archive2 = second / archive_name
        self.assertEqual(archive1.read_bytes(), archive2.read_bytes())

        checksum = first / (archive_name + ".sha256")
        expected = hashlib.sha256(archive1.read_bytes()).hexdigest()
        self.assertEqual(checksum.read_text(encoding="ascii"),
                         f"{expected}  {archive_name}\n")

        with tarfile.open(archive1, "r:gz") as tar:
            names = tar.getnames()
            self.assertTrue(all(name == ARCHIVE_ROOT or
                                name.startswith(ARCHIVE_ROOT + "/")
                                for name in names))
            executable = tar.getmember(ARCHIVE_ROOT + "/bin/generate.x")
            self.assertEqual(executable.mode, 0o755)
            link = tar.getmember(ARCHIVE_ROOT + "/lib/libaenet.so.2")
            self.assertTrue(link.issym())
            self.assertEqual(link.linkname, "libaenet.so.2.0.4")
            manifest = tar.extractfile(ARCHIVE_ROOT + "/manifest.txt")
            self.assertIsNotNone(manifest)
            contents = manifest.read().decode("utf-8")
            self.assertIn("bin/generate.x", contents)
            self.assertIn("lib/libaenet.so.2 -> libaenet.so.2.0.4", contents)
            self.assertIn("metadata.json", contents)
            self.assertIn("licenses/OpenBLAS.txt", contents)

        subprocess.run([str(VALIDATE), str(archive1)], check=True)

    def test_rejects_an_unapproved_archive_name(self):
        invalid = self.work / "aenet-2.0.4-macos-x86_64-gnu-serial"
        invalid.mkdir()
        result = self.run_package(self.work / "output", invalid, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("archive root", result.stderr)

    def test_rejects_a_symlink_that_escapes_the_archive(self):
        os.symlink("../../outside", self.stage / "lib" / "unsafe")
        result = self.run_package(self.work / "output", check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unsafe symbolic link", result.stderr)

    def test_rejects_incomplete_redistribution_notices(self):
        (self.stage / "licenses" / "OpenBLAS.txt").unlink()
        result = self.run_package(self.work / "output", check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("license files", result.stderr)

    def test_validation_rejects_a_changed_archive(self):
        output = self.work / "output"
        self.run_package(output)
        archive = output / (ARCHIVE_ROOT + ".tar.gz")
        archive.write_bytes(archive.read_bytes() + b"changed")
        result = subprocess.run(
            [str(VALIDATE), str(archive)], check=False,
            capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("checksum does not match", result.stderr)


class BuildEntryPointTests(unittest.TestCase):
    def fake_tool(self, work):
        log = work / "commands.txt"
        tool = work / "tool"
        tool.write_text(
            "#!/bin/sh\nprintf '%s\\n' \"$*\" >> \"$COMMAND_LOG\"\n",
            encoding="utf-8",
        )
        tool.chmod(0o755)
        environment = os.environ.copy()
        environment["COMMAND_LOG"] = str(log)
        return tool, log, environment

    def test_macos_build_uses_accelerate_and_serial_execution(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            tool, log, environment = self.fake_tool(work)
            result = subprocess.run(
                [str(BUILD), "--platform", "macos-arm64",
                 "--compiler", "/toolchain/gfortran-14",
                 "--build-dir", str(work / "build-tree"),
                 "--stage-parent", str(work / "stage"),
                 "--cmake", str(tool), "--ctest", str(tool)],
                check=True, capture_output=True, text=True, env=environment,
            )
            commands = log.read_text(encoding="utf-8").splitlines()
            self.assertIn("-DBLA_VENDOR=Apple", commands[0])
            self.assertIn("-DUSE_OPENBLAS=OFF", commands[0])
            self.assertIn("-DCMAKE_OSX_DEPLOYMENT_TARGET=15.0", commands[0])
            self.assertIn("--parallel 1", commands[1])
            self.assertIn("--output-on-failure", commands[2])
            root = f"aenet-{VERSION}-macos-arm64-gnu-serial"
            self.assertIn(root, commands[3])
            self.assertIn(root, result.stdout)

    def test_linux_build_selects_openblas(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            tool, log, environment = self.fake_tool(work)
            subprocess.run(
                [str(BUILD), "--platform", "linux-x86_64",
                 "--compiler", "/usr/bin/gfortran",
                 "--build-dir", str(work / "build-tree"),
                 "--stage-parent", str(work / "stage"),
                 "--cmake", str(tool), "--ctest", str(tool)],
                check=True, capture_output=True, text=True, env=environment,
            )
            configure = log.read_text(encoding="utf-8").splitlines()[0]
            self.assertIn("-DUSE_OPENBLAS=ON", configure)
            self.assertNotIn("-DBLA_VENDOR=Apple", configure)


if __name__ == "__main__":
    unittest.main()
