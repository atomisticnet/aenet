"""
To be used with the Atomic Simulation Environment (ASE).
See also: https://wiki.fysik.dtu.dk/ase
"""

from __future__ import print_function, division

__author__ = "Alexander Urban, Nongnuch Artrith"
__email__ = "aurban@atomistic.net, nartrith@atomistic.net"
__date__   = "2014-09-19"
__version__ = "0.1"

import sys
import numpy as np
try:
    from ase.calculators.calculator import all_changes
    from ase.calculators.calculator import Calculator
    from ase.neighborlist import NeighborList
except ImportError:
    sys.stderr.write("Error: could not import ASE "
                     "(https://wiki.fysik.dtu.dk/ase).")
    sys.exit()

from aenet.core import ANNPotentials

class ANNCalculator(Calculator):

    def __init__(self, potentials, **kwargs):
        Calculator.__init__(self, **kwargs)
        self.ann = ANNPotentials(potentials)
        self.cutoff = self.ann.Rc_max
        self.cutoff2 = self.cutoff**2
        self.neighbors = None
        self.implemented_properties = {
            'energy' : self.calculate_energy,
            'forces' : self.calculate_energy_and_forces,
        }
        self.results = {}

    def __del__(self):
        try:
            del self.ann
        except (AttributeError, NameError):
            pass

    def release(self):
        del self.ann

    def update(self, atoms):
        if ((self.neighbors is None) or 
            (len(self.neighbors.nl.cutoffs) != len(atoms))):
            cutoffs = self.cutoff*np.ones(len(atoms))
            self.neighbors = NeighborList(cutoffs,
                                          self_interaction=False,
                                          bothways=True)
        self.neighbors.update(atoms)

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):

        # Calculator essentially does: self.atoms = atoms
        Calculator.calculate(self, atoms, properties, system_changes)

        has_results = [p in self.results for p in properties]
        if (len(system_changes) > 0) or (not np.all(has_results)):
            self.update(self.atoms)
            if ('energy' in properties) and ('forces' in properties):
                # Forces evaluation requires energy. No need to compute
                # energy twice.
                del properties[properties.index('energy')]
            for p in properties:
                if p in self.implemented_properties:
                    self.implemented_properties[p](self.atoms)
                else:
                    raise NotImplementedError(
                        "Property not implemented: {}".format(p))

    def calculate_energy(self, atoms):
        energy = 0.0
        atom_types = atoms.get_chemical_symbols()
        for i in range(len(atoms)):
            indices, offsets = self.neighbors.get_neighbors(i)
            type_i = atom_types[i]
            coords_i = atoms.positions[i]
            coords_j = np.empty((len(indices),3), dtype=np.double)
            types_j = []
            inb = 0
            for j, offset in zip(indices, offsets):
                coo = (atoms.positions[j] + np.dot(offset, atoms.get_cell()))
                d2 = np.sum((coo - coords_i)**2)
                if d2 <= self.cutoff2:
                    coords_j[inb] = coo
                    types_j.append(atom_types[j])
                    inb += 1
            energy += self.ann.atomic_energy(coords_i, type_i, coords_j[:inb], types_j)
        self.results['energy'] = energy

    def calculate_energy_and_forces(self, atoms):
        energy = 0.0
        atom_types = atoms.get_chemical_symbols()
        forces = np.zeros((len(atoms), 3), dtype=np.double)
        for i in range(len(atoms)):
            indices, offsets = self.neighbors.get_neighbors(i)
            type_i = atom_types[i]
            index_i = i + 1
            index_j = np.empty(indices.size, dtype=np.intc)
            coords_i = atoms.positions[i]
            coords_j = np.empty((indices.size,3), dtype=np.double)
            types_j = []
            inb = 0
            for j, offset in zip(indices, offsets):
                coo = (atoms.positions[j] + np.dot(offset, atoms.get_cell()))
                d2 = np.sum((coo - coords_i)**2)
                if d2 <= self.cutoff2:
                    coords_j[inb] = coo
                    index_j[inb] = j + 1
                    types_j.append(atom_types[j])
                    inb += 1
            energy += self.ann.atomic_energy_and_forces(
                coords_i, type_i, index_i, coords_j[:inb], types_j,
                index_j[:inb], forces)
        self.results['energy'] = energy
        self.results['forces'] = forces
