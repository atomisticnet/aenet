#!/usr/bin/env python3
# This file is part of the AENET package.
#
# Copyright (C) 2012-2026 Nongnuch Artrith and Alexander Urban
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.

import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
VERIFY = ROOT / "packaging" / "verify_release_tag"


class ReleaseTagTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.repository = Path(self.temporary.name)
        subprocess.run(["git", "init", "-q", self.repository], check=True)
        subprocess.run(
            ["git", "-C", self.repository, "config", "user.name", "Test"],
            check=True,
        )
        subprocess.run(
            ["git", "-C", self.repository, "config", "user.email",
             "test@example.com"],
            check=True,
        )
        (self.repository / "src").mkdir()
        (self.repository / "src" / "VERSION").write_text("2.1.0\n")
        subprocess.run(
            ["git", "-C", self.repository, "add", "src/VERSION"], check=True
        )
        subprocess.run(
            ["git", "-C", self.repository, "commit", "-qm", "initial"],
            check=True,
        )

    def tearDown(self):
        self.temporary.cleanup()

    def verify(self, tag):
        return subprocess.run(
            [VERIFY, "--repository", self.repository, tag],
            text=True,
            capture_output=True,
        )

    def test_accepts_annotated_tag_matching_target_version(self):
        subprocess.run(
            ["git", "-C", self.repository, "tag", "-am", "release", "v2.1.0"],
            check=True,
        )
        result = self.verify("v2.1.0")
        expected = subprocess.check_output(
            ["git", "-C", self.repository, "rev-parse", "HEAD"], text=True
        ).strip()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), expected)

    def test_rejects_lightweight_tag(self):
        subprocess.run(
            ["git", "-C", self.repository, "tag", "v2.1.0"], check=True
        )
        result = self.verify("v2.1.0")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("annotated", result.stderr)

    def test_rejects_version_mismatch(self):
        subprocess.run(
            ["git", "-C", self.repository, "tag", "-am", "release", "v2.0.9"],
            check=True,
        )
        result = self.verify("v2.0.9")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("src/VERSION", result.stderr)

    def test_rejects_missing_tag(self):
        result = self.verify("v2.1.0")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("does not exist", result.stderr)

    def test_rejects_non_release_tag_name(self):
        subprocess.run(
            ["git", "-C", self.repository, "tag", "-am", "test", "testing"],
            check=True,
        )
        result = self.verify("testing")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("vMAJOR.MINOR.PATCH", result.stderr)


if __name__ == "__main__":
    unittest.main()
