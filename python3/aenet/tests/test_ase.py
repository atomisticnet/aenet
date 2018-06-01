"""Unit tests for the ASE Calculator
"""

from __future__ import print_function, division

__author__ = "Alexander Urban"
__email__ = "alexurba@mit.edu"
__date__ = "2014-09-19"
__version__ = "0.1"

import unittest
import os
import numpy as np

try:
    import ase
    HAS_ASE = True
except ImportError:
    HAS_ASE = False

from aenet.ase_calculator import ANNCalculator

POTENTIAL_Ti = os.path.join(os.path.dirname(__file__), "Ti.15t-15t.nn")
POTENTIAL_O = os.path.join(os.path.dirname(__file__), "O.15t-15t.nn")

potential_files = {"Ti": POTENTIAL_Ti, "O": POTENTIAL_O}


class ASETest(unittest.TestCase):

    def get_atoms(self):
        positions = [(0.42651082, 2.82851745, 1.22027956), (1.76480257, -0.94296447,
                                                            3.66083922),
                     (0.14225662, 0.94278546, 0.79509432), (2.04905661, 0.94276737,
                                                            4.08602446),
                     (-0.14223328, -0.94276538,
                      3.23565398), (2.33354646, 2.82831821, 1.64546480)]
        cell = [(3.81382872, 0.00000000, 0.00000000),
                (0.56886601, 3.77116440, 0.00000000), (-2.19138173, -1.88561177,
                                                       4.88111931)]
        atoms = ase.Atoms(
            "Ti2O4", positions=positions, cell=cell, pbc=[True, True, True])
        return atoms

    @unittest.skipIf(not HAS_ASE, "ASE not available")
    def test_init(self):
        atoms = self.get_atoms()
        calc = ANNCalculator(potential_files)
        atoms.set_calculator(calc)
        calc.release()

    @unittest.skipIf(not HAS_ASE, "ASE not available")
    def test_energy(self):
        atoms = self.get_atoms()
        calc = ANNCalculator(potential_files)
        atoms.set_calculator(calc)
        energy = atoms.get_potential_energy()
        self.assertAlmostEqual(energy, -4990.310783539421)
        calc.release()

    @unittest.skipIf(not HAS_ASE, "ASE not available")
    def test_forces(self):
        atoms = self.get_atoms()
        calc = ANNCalculator(potential_files)
        atoms.set_calculator(calc)
        forces = atoms.get_forces()
        num_forces = self._get_num_forces(atoms)
        for iat in range(len(atoms)):
            for i in range(3):
                self.assertAlmostEqual(
                    forces[iat, i], num_forces[iat, i], places=2)
        calc.release()

    def _get_num_forces(self, atoms):
        d = 0.01
        dd = np.identity(3) * d
        forces = np.zeros((len(atoms), 3))
        for iat in range(len(atoms)):
            for i in range(3):
                atoms.positions[iat] -= dd[i]
                E1 = atoms.get_potential_energy()
                atoms.positions[iat] += 2 * dd[i]
                E2 = atoms.get_potential_energy()
                atoms.positions[iat] -= dd[i]
                forces[iat, i] = -(E2 - E1) / (2 * d)
        return forces


if __name__ == "__main__":
    unittest.main()
