"""
Tests the aenet interface in core.py.
"""

from __future__ import print_function, division

__author__ = "Nonguch Artrith, Alexander Urban"
__email__ = "nartrith@mit.edu, alexurba@mit.edu"
__date__   = "2014-09-11"
__version__ = "0.1"

import unittest
import os

from aenet.core import ANNPotentials

POTENTIAL_Ti = os.path.join(os.path.dirname(__file__), 'Ti.10t-10t.nn')
POTENTIAL_O  = os.path.join(os.path.dirname(__file__), 'O.10t-10t.nn')
POTENTIAL_H  = os.path.join(os.path.dirname(__file__), 'H.10t-10t.nn')
potential_files = {
    'Ti' : POTENTIAL_Ti,
    'O'  : POTENTIAL_O,
    'H'  : POTENTIAL_H,
}

class AenetCoreTest(unittest.TestCase):

    def test_load_potentials(self):
        pot = ANNPotentials(potential_files)
        self.assertEqual(pot.Rc_min, 0.75)
        self.assertEqual(pot.Rc_max, 6.5)
        del pot

    def test_atomic_energy(self):
        pot = ANNPotentials(potential_files)
        self.assertEqual(pot.Rc_min, 0.75)
        self.assertEqual(pot.Rc_max, 6.5)
        del pot

if __name__=="__main__":
    unittest.main()
