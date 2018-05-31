#!/usr/bin/env python
"""
An example script using aenetLib.

This script requires the Atomistic Simulation Environment (ASE)
(obtainable from: https://wiki.fysik.dtu.dk/ase).

The input file that points to the ANN potential files is in JSON format.

Example input file:

{
    "potentials" : {
        "Ti" : Ti.10t-10t.nn,
        "O"  : O.10t-10t.nn,
        "H"  : H.10t-10t.nn
    }
}
"""

from __future__ import print_function, division

__author__ = "Alexander Urban, Nongnuch Artrith"
__email__ = "alexurba@mit.edu, nartrith@mit.edu"
__date__   = "2014-09-15"
__version__ = "0.1"

import sys
import argparse
import json
import numpy as np

try:
    import ase.io
    from ase.optimize import BFGS
except ImportError:
    sys.stderr.write("Error: unable to import ASE.")
    sys.exit()

from aenet.ase_calculator import ANNCalculator

#----------------------------------------------------------------------#

def get_ANN_energy(input_file, structure_files, structure_format,
                   calc_forces, do_relax):
    if do_relax is not None:
        calc_forces = True
    with open(input_file, 'r') as fp:
        inp = json.load(fp)
    calc = ANNCalculator(inp['potentials'])
    for strucfile in structure_files:
        atoms = ase.io.read(strucfile, format=structure_format)
        atoms.set_calculator(calc)
        if calc_forces:
            forces = atoms.get_forces()
        energy = atoms.get_potential_energy()
        symbol = atoms.get_chemical_symbols()
        print("Structure file: {}".format(strucfile))
        if calc_forces:
            print("Atomic positions and forces:")
        else:
            print("Atomic positions:")
        for i in range(len(atoms)):
            print(("{:2s}" + 3*"{:15.8f} ").format(
                symbol[i], *atoms.positions[i]), end="")
            if calc_forces:
                print((3*"{:15.8f} ").format(*forces[i]))
            else:
                print("")
        if calc_forces:
            print("Average force (must be zero): "
                  + (3*"{:15.8f} ").format(*np.sum(forces, axis=0)))
        print("Total energy = {} eV".format(energy))
        if do_relax is not None:
            relax = BFGS(atoms)
            relax.run(fmax=do_relax)
            ext = strucfile.split(".")[-1]
            outfile = ".".join(strucfile.split(".")[:-1]) + "-relaxed." + ext
            ase.io.write(outfile, atoms, format=structure_format)

#----------------------------------------------------------------------#

if (__name__ == "__main__"):

    parser = argparse.ArgumentParser(
        description     = __doc__+"\n{} {}".format(__date__,__author__),
        formatter_class = argparse.RawDescriptionHelpFormatter )

    parser.add_argument(
        "input_file",
        help    = "Path to the principle input file with parameters.")

    parser.add_argument(
        "structure_files",
        help    = "Path(s) to atomic structure file(s).",
        nargs   = "+")

    parser.add_argument(
        "--format",
        help    = "Structure format known by ASE (default: auto-detect).",
        type    = str,
        default = None)

    parser.add_argument(
        "--forces", "-f",
        help    = "Calculate atomic force.",
        action  = "store_true")

    parser.add_argument(
        "--relax",
        help    = "Perform geometry optimization to given force criterion.",
        type    = float,
        default = None)

    args = parser.parse_args()

    get_ANN_energy(input_file=args.input_file,
                   structure_files=args.structure_files,
                   structure_format=args.format,
                   calc_forces=args.forces,
                   do_relax=args.relax)
