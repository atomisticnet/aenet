#!/usr/bin/env python3
# This file is part of the AENET package.
#
# Copyright (C) 2012-2026 Nongnuch Artrith and Alexander Urban
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.

import argparse
import math
import os
from pathlib import Path
import re
import subprocess


def run(root, work, name, arguments, command_prefix):
    executable = root / "bin" / f"{name}.x"
    result = subprocess.run(
        [*command_prefix, str(executable), *arguments], cwd=work, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90,
    )
    (work / f"{name}.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0 or "Error:" in result.stdout:
        raise RuntimeError(f"{name} failed:\n{result.stdout[-2000:]}")
    return result.stdout


def main():
    parser = argparse.ArgumentParser(
        description="Run the binary release numerical CLI smoke workflow."
    )
    parser.add_argument("root", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument(
        "--sandbox-profile",
        help="run each backend executable through sandbox-exec with this profile",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    command_prefix = []
    if args.sandbox_profile:
        command_prefix = ["sandbox-exec", "-p", args.sandbox_profile]
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    (work / "Cu.stp").write_text(
        "DESCR\nSmoke\nEND DESCR\nATOM Cu\nENV 1\nCu\nRMIN 0.5\n"
        "BASIS type=chebyshev\n"
        "radial_Rc=5 radial_N=3 angular_Rc=5 angular_N=2\n",
        encoding="utf-8",
    )
    for index in range(8):
        distance = 2.0 + index * 0.15
        energy = (distance - 2.5) ** 2 - 1.0
        force = 2.0 * (distance - 2.5)
        (work / f"{index}.xsf").write_text(
            f"# total energy = {energy} eV\nATOMS\n"
            f"Cu 0 0 0 {force} 0 0\n"
            f"Cu {distance} 0 0 {-force} 0 0\n",
            encoding="utf-8",
        )
    files = "".join(f"{index}.xsf\n" for index in range(8))
    (work / "generate.in").write_text(
        "OUTPUT smoke.train\nTYPES\n1\nCu 0\nSETUPS\nCu Cu.stp\n"
        f"FILES\n8\n{files}", encoding="utf-8",
    )
    run(root, work, "generate", ["generate.in"], command_prefix)
    if (work / "smoke.train").stat().st_size == 0:
        raise RuntimeError("generate produced an empty training set")

    (work / "train.in").write_text(
        "TRAININGSET smoke.train\nTESTPERCENT 0\nITERATIONS 2\n"
        "METHOD\nbfgs\nNETWORKS\nCu Cu.nn 1 3:tanh\n", encoding="utf-8",
    )
    run(root, work, "train", ["train.in"], command_prefix)
    if (work / "Cu.nn").stat().st_size == 0:
        raise RuntimeError("train produced an empty network")

    (work / "predict.in").write_text(
        "TYPES\n1\nCu\nNETWORKS\nCu Cu.nn\nFORCES\n", encoding="utf-8"
    )
    output = run(root, work, "predict", ["predict.in", "3.xsf"], command_prefix)
    matches = re.findall(r"Total energy\s*:\s*([-+0-9.Ee]+)", output)
    if len(matches) != 1 or not math.isfinite(float(matches[0])):
        raise RuntimeError("predict did not report one finite total energy")
    energy = float(matches[0])
    repeated = run(
        root, work, "predict", ["predict.in", "3.xsf"], command_prefix
    )
    repeated_energy = float(re.findall(
        r"Total energy\s*:\s*([-+0-9.Ee]+)", repeated
    )[0])
    if abs(energy - repeated_energy) >= 1.0e-8:
        raise RuntimeError("repeated prediction changed by 1e-8 eV or more")
    print(f"generate/train/predict passed; repeated energy: {energy}")


if __name__ == "__main__":
    main()
