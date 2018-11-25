"""Utilities for interacting with aenet xsf files.

aenet uses xsf (http://www.xcrysden.org/doc/XSF.html) for structure data. This
module provides two utilities to read aenet-compatible xsf files into ase.Atoms
objects, and to write ase.Atoms objects into aenet-compatible xsf files.

The minor variation of xsf that aenet uses is that the energy is stored in a
comment in the file.

"""

import os
import re
import ase.io


def read_xsf(xsfile):
    """Return an atoms object for the xsfile.

    Parameters
    ----------
    xsfile : string
        Filename to read xsf data from

    Returns
    -------
    atoms : ase.Atoms
        The energy and forces are stored in a SinglePointCalculator on
        the returned atoms.

    """
    atoms = ase.io.read(xsfile)

    with open(xsfile, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith('# total energy'):
                m = re.findall(r'[-+]?\d*\.\d+|\d+', line)
                energy = float(m[0])
                atoms.calc.results['energy'] = energy
                break

    return atoms


def write_xsf(xsfile, atoms):
    """Create an aenet compatible xsf from an atoms object.

    Parameters
    ----------
    xsfile : string
        Filename to write xsf file to

    atoms : ase.Atoms
        an ase atoms object with an attached calculator to get energy and
        forces.

    Returns
    -------
    output : string
        The string written to the file.
    """
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    xsf = ['# total energy = {} eV'.format(energy), '']

    if True in atoms.pbc:
        xsf += ['CRYSTAL', 'PRIMVEC']
        for v in atoms.get_cell():
            xsf += ['{} {} {}'.format(*v)]
        xsf += ['PRIMCOORD', '{} 1'.format(len(atoms))]

    else:
        xsf += ['ATOMS']

    S = ('{atom.symbol:<3s} {atom.x: .12f} {atom.y: .12f} {atom.z: .12f}'
         ' {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}')
    xsf += [S.format(atom=atom, f=forces[i]) for i, atom in enumerate(atoms)]

    output = '\n'.join(xsf)

    base, fname = os.path.split(xsfile)

    if not os.path.isdir(base):
        os.makedirs(base)

    with open(xsfile, 'w') as f:
        f.write(output)

    return output
