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

try:
    import ase.io
except ImportError:
    sys.stderr.write("Error: unable to import ASE.")
    sys.exit()

from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units

from aenet.ase_calculator import ANNCalculator


def printenergy(istep, a):
    """ print the potential, kinetic and total energy """
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print("{:5d} {:15.8e} {:15.8e} {:7.1f} {:15.8e}".format(
        istep, epot, ekin, ekin/(1.5*units.kB), epot+ekin))
    sys.stdout.flush()


def parse_input(input_file):
    # defaults
    trajectory_file = "md.traj"
    # read parameters from input file
    with open(input_file, 'r') as fp:
        inp = json.load(fp)
    structure_file = str(inp['structure_file'])
    T = float(inp['temperature'])
    dt = float(inp['time_step'])
    md_steps = int(inp['md_steps'])
    print_steps = int(inp['print_steps'])
    potentials = inp['potentials']
    if "trajectory_file" in inp:
        trajectory_file = str(inp['trajectory_file'])
    return (structure_file, T, dt, md_steps, print_steps,
            trajectory_file, potentials)


def get_ANN_energy(input_file):
    (structure_file, T, dt, md_steps, print_steps,
     trajectory_file, potentials) = parse_input(input_file)
    # atomic structure
    atoms = ase.io.read(structure_file, format='vasp')
    # ANN calculator
    calc = ANNCalculator(potentials)
    atoms.set_calculator(calc)
    # initialize velocities
    MaxwellBoltzmannDistribution(atoms, temp=T*units.kB)
    # initialize MD
    md = VelocityVerlet(atoms, dt*units.fs, trajectory=trajectory_file)
    print("# {:5s} {:15s} {:15s} {:7s} {:15s}".format(
        "step", "E_pot", "E_kin", "T", "E_tot"))
    printenergy(0, atoms)
    istep = 0
    for i in range(int(md_steps/print_steps)):
        md.run(steps=print_steps)
        istep += print_steps
        printenergy(istep, atoms)


if (__name__ == "__main__"):

    parser = argparse.ArgumentParser(
        description=__doc__+"\n{} {}".format(__date__, __author__),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "input_file",
        help="Path to the principle input file with parameters.")

    args = parser.parse_args()

    get_ANN_energy(input_file=args.input_file)
