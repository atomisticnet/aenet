"""Utility functions for interacting with aenet output files.


"""

import re
import numpy as np
from io import StringIO
import hashlib


def get_training_details(outfile='train.out'):
    """Parse the training output file to get the convergence data.

    Parameters
    ----------
    outfile : string, default=train.out

    Returns
    -------
    epochs, train_mae, train_rmse, test_mae, test_rmse
    """
    with open(outfile, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if ' Number of iterations' in line:
            nepochs = int(line.replace(" Number of iterations    :", ""))

        if ' epoch             MAE          <RMSE>' in line:
            break

    # Here we break out the section containing the data, and then use numpy to
    # read it in. In case the file is not done yet, we make sure to only go to
    # the end.
    end = min(i + nepochs + 2, len(lines))
    s = StringIO("\n".join(lines[i + 1: end]))
    (epochs, train_mae, train_rmse,
     test_mae, test_rmse, _) = np.loadtxt(s, dtype={'names': ('epoch',
                                                              'train-mae',
                                                              'train-rmse',
                                                              'test-mae',
                                                              'test-rmse',
                                                              'c'),
                                                    'formats': ('i', 'f',
                                                                'f', 'f',
                                                                'f', 'S1')},
                                          unpack=True)

    return epochs, train_mae, train_rmse, test_mae, test_rmse


def read_ascii_train_file(ascii_train='train.ascii'):
    """Parse the ascii train file into a dictionary.

    This contains information about the dataset, and each configuration in the
    dataset, including the fingerprints. The atomic numbers, positions and
    forces for each configuration is also in the output. There are no unit cells
    for the configurations though, that does not appear to be in the train file.

    Parameters
    ----------
    ascii_train : string, default='train.ascii'
        filename to the ascii train file. This should be created with
        trnset2ASCII.x

    Returns
    -------
    a dictionary containing the data.
    {'header': {information about the computation},
     'configurations': [{information about each configuration},]}

    """
    f = open(ascii_train, 'r')
    fname = f.readline()
    scaled = f.readline().strip()
    scale = float(f.readline())
    shift = float(f.readline())

    nelements = int(f.readline())
    elements = [f.readline().strip() for i in range(nelements)]
    atomic_energies = [float(x) for x in f.readline().split()]
    total_atoms = float(f.readline())
    nconfigurations = int(f.readline())
    emin, emax, e_avg = [float(x) for x in f.readline().split()]

    data = {
        'header':
            dict(
                fname=fname,
                scaled=scaled,
                scale=scale,
                shift=shift,
                nelements=nelements,
                elements=elements,
                atomic_energies=atomic_energies,
                total_atoms=total_atoms,
                nconfigurations=nconfigurations,
                emin=emin,
                emax=emax,
                e_avg=e_avg),
        'configurations': []
    }

    for i in range(nconfigurations):
        path = f.readline().strip()
        natoms, nspecies = [int(x) for x in f.readline().split()]
        energy = float(f.readline())

        configuration = dict(
            path=path, energy=energy, natoms=natoms, nspecies=nspecies)

        positions, forces = [], []
        FP = []
        numbers = []
        for j in range(natoms):
            species_number = int(f.readline())
            numbers += [species_number]
            x, y, z = [float(pos) for pos in f.readline().split()]
            positions += [[x, y, z]]

            sfx, sfy, sfz = [float(pos) for pos in f.readline().split()]
            forces += [[sfx, sfy, sfz]]
            nfingerprints = int(f.readline())
            fp = [float(x) for x in f.readline().split()]
            if len(fp) != nfingerprints:
                raise Exception('Found wrong # of fingerprints')
            FP += [fp]

        configuration['fingerprint'] = FP
        configuration['numbers'] = numbers
        configuration['positions'] = positions
        configuration['forces'] = forces
        data['configurations'] += [configuration]

    return data


def read_predict_out(predictfile='predict.out'):
    """Read the predict.out results.

    This file should be the output from predict.x.

    Parameters
    ----------
    predictfile : string, default='predict.out'


    Returns
    -------
     a zip iterator of (filename, natoms, energy)
    """
    files, natoms, energies = [], [], []
    energy_evaluation = False

    with open(predictfile) as f:
        # Get to the energy evaluation section
        for line in f:
            if not energy_evaluation and 'Energy evaluation' in line:
                energy_evaluation = True

            if energy_evaluation:
                if 'File name' in line:
                    files.append(line.replace(
                        ' File name         : ', '').strip())
                if 'Number of atoms' in line:
                    natoms.append(
                        int(line.replace(' Number of atoms   :', '')))
                if ' Total energy               :' in line:
                    m = re.findall(r'[-+]?\d*\.\d+|\d+', line)
                    energies.append(float(m[0]))
                if 'Atomic Energy Network done.' in line:
                    break

    return zip(files, natoms, energies)


def aenet_hash(atoms):
    """Get an sha1 hash for the atoms.

    This is a string that is practically unique for each atoms. It is a hash of
    a string containing the atomic numbers, positions, cell and periodic
    boundary conditions.

    Parameters
    ----------

    atoms : ase.Atoms

    Returns
    -------
    The hash

    """

    s = f'{atoms.numbers}{atoms.positions}{atoms.cell}{atoms.pbc}'
    return hashlib.sha1(s.encode('utf-8')).hexdigest()


def aenet_write_db(db, atoms, **kwargs):
    """Write atoms to a database if it has not been written before.
    The hash is automatically added for you.

    Parameters
    ----------

    atoms : ase.Atoms

    **kwargs : extra keywords to pass to db.write.

    Returns
    -------

    """
    hsh = aenet_hash(atoms)
    try:
        row = db.get(hash=hsh)
        return row
    except KeyError:
        return db.write(atoms, hash=hsh, **kwargs)
