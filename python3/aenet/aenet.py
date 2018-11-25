'''aenet python interface.

This module helps automate creation of the input files.
'''

import glob
import json
import os
from subprocess import Popen, PIPE
from aenet.ase_calculator import ANNCalculator
import numpy as np
from io import StringIO
from itertools import combinations_with_replacement

GENERATEX = 'generate.x'
TRAINX = 'train.x'
PREDICTX = 'predict.x'


def aenet_save(aenet, fname=None):
    """Save the aenet dictionary to disk.

    Parameters
    ----------

    aenet : dict
        This is a parameter dictionary.

    fname : str, optional, default=None
        The filename to save to. If None, defaults to label.json

    Returns
    -------
        None
    """

    if fname is None:
        fname = aenet['label'] + '.json'
    with open(fname, 'w') as f:
        f.write(json.dumps(aenet, indent=4))


def aenet_load(fname=None):
    """load the aenet dictionary from disk.

    Parameters
    ----------

    fname : str optional, default=None
        The file name to load.

    Returns
    -------
        The dictionary loaded from the file.
    """
    with open(fname) as f:
        return json.loads(f.read())


def aenet(label='aenet', elements=(), xsfiles=()):
    """Initialization function for aenet.

    Creates and saves a json file containing information that is needed for
    setting up and running aenet computations.

    Parameters
    ----------

    label : string
        A label for the simulation

    elements : list
        list of tuples containing (symbol, atomic_energy) for each element

    xsfiles : list
         a list of paths to xsfiles to be used for training.

    Returns
    -------
    A dictionary that is the first argument for all the aenet functions.

    """

    p = {'label': label,
         'elements': elements,
         'xsfiles': xsfiles,
         'setups': {},
         'generate': {},
         'train': {},
         'predict': {}}
    aenet_save(p)
    return p


def symmfunc_centered(aenet, Rc=6.5, g2_etas=(),
                      g4_etas=(), g4_lambdas=(), g4_zetas=()):
    elements = [sym for sym, _ae in aenet['elements']]

    radial_fingerprints = [[('G', 2), ('type2', sym),
                            ('eta', eta), ('Rs',  0.0000), ('Rc', Rc)]
                           for eta in g2_etas for sym in elements]

    # angular fingerprints
    angular_fingerprints = []

    for zeta in g4_zetas:
        for _lambda in g4_lambdas:
            for eta in g4_etas:
                for (t2, t3) in set(combinations_with_replacement(elements, 2)):
                    fp = [('G', 4), ('type2', t2), ('type3', t3),
                          ('eta',  eta), ('lambda', _lambda),
                          ('zeta', zeta), ('Rc', Rc)]
            angular_fingerprints += [fp]

    fingerprints = radial_fingerprints + angular_fingerprints
    s = '\n'.join(['  '.join(['{}={}'.format(key, val) for key, val in fp])
                   for fp in fingerprints])

    sf = f'''SYMMFUNC type=Behler2011
{len(fingerprints)}
{s}\n'''
    return sf


def chebyshev_basis(aenet, radial_rc=4.0, radial_n=6,
                    angular_rc=4.0, angular_n=2):
    """Returns the Chebshev basis section of a fingerprint file.

    Parameters
    ----------

    aenet : dict
        The aenet state dictionary.

    radial_rc : float optional, default=4.0
        Cutoff radius for the radial part.

    radial_n : int optional, default=6
        Number of radial Chebyshev functions.

    angular_rc : float optional, default=4.0
        Cutoff radius for the angular part.

    angular_n : int optional, default=2
        Number of angular Chebyshev functions

    Returns
    -------
        String of the section.
    """
    return f'''BASIS type=Chebyshev
radial_Rc = {radial_rc} radial_N = {radial_n} angular_Rc = {angular_rc} angular_N = {angular_n}'''


def generate_setup(aenet,
                   element,
                   setupfile=None,
                   description=None,
                   rmin=0.75,
                   fingerprints=None):
    """Generate an element setup file.

    This will create the setup, and write it to a file. The aenet dictionary
    will be updated and saved.

    Parameters
    ----------

    aenet : dict
        The state dictionary.

    element : string
        string for chemical symbol to make a setup for

    setupfile : string
        filename to write setup to. Defaults to {element}.stp

    description : string, optional, default=None.
        Short description of the fingerprint.

    rmin : float
        Minimum distance to consider.

    fingerprints : string
        This should be a string representation of the symmetry functions or
        Chebyshev basis. Usually this should be the output of a function that
        generates these.

    Returns
    -------
        Returns the string representing the setup.

    """

    if setupfile is None:
        setupfile = element + '.stp'

    if description is None:
        description = 'Setup for ' + element

    environment = [sym for sym, atomic_energy in aenet['elements']]

    p = {'setupfile': setupfile,
         'description': description,
         'environment': environment,
         'rmin': rmin,
         'fingerprints': fingerprints}

    aenet['setups'][element] = p
    aenet_save(aenet)

    out = ['DESCR',
           description,
           'END DESCR',
           '',
           'ATOM {}'.format(element),
           '',
           'ENV {}'.format(len(environment))]
    for sym in environment:
        out += [sym]

    out += ['',
            'RMIN {}d0'.format(float(rmin)),
            '',
            fingerprints]

    output = '\n'.join(out)
    with open(setupfile, 'w') as f:
        f.write(output)

    return output


def generate_in(aenet,
                inputfile='generate.in',
                trainfile=None,
                debug=None,
                timing=None):
    """Create the input file for generation.

    http://ann.atomistic.net/Documentation/#training-set-generation-with-generate-x

    Parameters
    ----------

    inputfile : string, defaults to generate.in
        The name of the file to output

    trainfile : string, defaults to {aenet['label']}.train.
        The name of the file to save the train data in.

    debug : boolean
        if True, put DEBUG in the input file.

    timint : boolean
        if True, put TIMING in the input file.

    Returns
    -------
    The contents of the generate.in file
    """

    if trainfile is None:
        trainfile = '{}.train'.format(aenet['label'])

    p = {'inputfile': inputfile,
         'trainfile': trainfile,
         'debug': debug,
         'timing': timing}

    aenet['generate'] = p
    aenet_save(aenet)

    types = aenet['elements']
    setups = []
    for sym, ae in types:
        setups += [(sym, aenet['setups'][sym]['setupfile'])]

    xsf = aenet['xsfiles']

    out = ['OUTPUT {}'.format(trainfile),
           '']

    if debug is not None:
        out += ['DEBUG']

    if timing is not None:
        out += ['TIMING']

    out += ['TYPES',
            '{}'.format(len(types))]

    for t, e in types:
        out += ['{0:3s} {1}'.format(t, e)]

    out += ['',
            'SETUPS']

    for t, f in setups:
        out += ['{0:3s} {1}'.format(t, f)]

    out += ['',
            'FILES',
            '{}'.format(len(xsf))]
    for p in xsf:
        out += [p]

    output = '\n'.join(out)
    with open(inputfile, 'w') as f:
        f.write(output)

    return output


def generate_x(aenet,
               inputfile=None,
               txtfile=None):
    """Run generate.x on the input file.

    Parameters
    ----------

    inputfile : string, defaults to {aenet['generate']['inputfile']}

    txtfile : string, defaults to {inputfile + '.out'}.

    Returns
    -------
        None
    """

    if inputfile is None:
        inputfile = aenet['generate']['inputfile']
    if txtfile is None:
        txtfile = inputfile + '.out'

    aenet['generate']['txtfile'] = txtfile
    aenet_save(aenet)
    cmd = [GENERATEX, inputfile]

    if os.path.exists(aenet['generate']['trainfile']):
        os.unlink(aenet['generate']['trainfile'])

    print('Running {} to generate {}'.format(
        ' '.join(cmd), aenet['generate']['trainfile']))

    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    output, stderr = process.communicate()

    with open(txtfile, 'wb') as f:
        f.write(output)

    print('Finished. Now you can setup the training file.')


def train_in(aenet,
             inputfile='train.in',
             trainfile=None,
             testpercent=10,
             method='bfgs',
             iterations=10,
             networks=(),
             debug=None,
             maxenergy=0.0,
             saveenergies=True,
             timing=None):
    """Generate the train.in file used by train.x.

    http://ann.atomistic.net/Documentation/#ann-potential-training-with-train-x

    Parameters
    ----------

    aenet : type
        The state dictionary.

             inputfile='train.in',
             trainfile=None,
             testpercent=10,
             method='bfgs',
             iterations=10,
             networks=(),
             debug=None,
             maxenergy=0.0,
             saveenergies=None,
             timing=None,


    inputfile: string for filename to write the input to.

    testpercent: int, percent of training data to use for testing

    method is a string including parameters.

    iterations: int, number of epochs to train

    networks: list of (symbol, ((#neurons, act), ...))


    Returns
    -------
    return
    """

    if trainfile is None:
        trainfile = aenet['generate']['trainfile']

    out = ['TRAININGSET {}'.format(trainfile),
           'TESTPERCENT {}'.format(testpercent),
           'ITERATIONS {}'.format(iterations),
           '']

    if maxenergy is not None:
        out += ['MAXENERGY {}'.format(maxenergy), '']

    if saveenergies is not None:
        out += ['SAVE_ENERGIES', '']

    if timing is not None:
        out += ['TIMING']

    if debug is not None:
        out += ['DEBUG']

    out += ['',
            'METHOD',
            method,
            '']

    out += ['NETWORKS',
            '! atom   network         hidden',
            '! types  file-name       layers  nodes:activation']

    syms = [x[0] for x in networks]
    nns = [x[1] for x in networks]
    paths = []

    nn_string = '{:8s} {:<15s} {:<7d} {}'
    for sym, nn in zip(syms, nns):
        path = f'{sym}.nn'
        paths.append(path)
        out += [nn_string.format(sym,
                                 path,
                                 len(nn),
                                 ' '.join(['{}:{}'.format(n, act) for n,
                                           act in nn]))]

    output = '\n'.join(out)

    p = {'inputfile': inputfile,
         'testpercent': testpercent,
         'method': method,
         'iterations': iterations,
         'networks': list(zip(syms, paths, nns)),
         'debug': debug,
         'maxenergy': maxenergy,
         'saveenergies': saveenergies,
         'timing': timing}

    aenet['train'] = p
    aenet_save(aenet)

    with open(inputfile, 'w') as f:
        f.write(output)

    return output


def train_x(aenet,
            inputfile=None,
            txtfile=None,
            save=False):
    """Run train.x"""

    if inputfile is None:
        inputfile = aenet['train']['inputfile']

    if txtfile is None:
        txtfile = inputfile + '.out'

    aenet['train']['outputfile'] = txtfile
    aenet_save(aenet)

    cmd = [TRAINX, inputfile, '|', 'tee', txtfile]

    print('Starting the training process.')
    for f in ['TRAIN.0', 'TEST.0']:
        if os.path.exists(f):
            os.unlink(f)

    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    output, stderr = process.communicate()

    with open(txtfile, 'wb') as f:
        f.write(output)

    if not save:
        for nn in glob.glob('*.nn-*'):
            os.unlink(nn)

    print('Done with training.')


def get_ann(aenet):
    """Return an ase calculator that uses the trained network."""
    potentials = {}
    for sym, path, nn in aenet['train']['networks']:
        potentials[sym] = path

    return ANNCalculator(potentials)


def predict_in(aenet,
               inputfile='predict.in',
               forces=None,
               relax=None,
               timing=None,
               debug=None,
               xsfiles=()):
    """Creates the input file for predict.x

    xsfiles is a list of xsf files.
    If empty, pass a file at the command line.

    http://ann.atomistic.net/Documentation/#prediction-of-structural-energies-and-atomic-forces-with-predict-x

    """

    p = {'inputfile': inputfile,
         'forces': forces,
         'relax': relax,
         'timing': timing,
         'debug': debug,
         'xsfiles': xsfiles}

    aenet['predict'] = p
    aenet_save(aenet)

    types = [sym for sym in aenet['elements']]
    networks = aenet['train']['networks']

    out = ['TYPES',
           '{}'.format(len(types))]

    for sym, ae in types:
        out += [sym]

    out += ['',
            'NETWORKS']

    for sym, path, nodes in networks:
        out += ['{:3s} {}'.format(sym, path)]

    out += ['']
    if forces is not None:
        out += ['FORCES']

    if relax is not None:
        out += ['RELAX']
        out += [relax]

    if timing is not None:
        out += ['TIMING']

    if xsfiles:
        out += ['',
                'FILES',
                '{}'.format(len(xsfiles))]

        for xsf in xsfiles:
            out += [xsf]

    output = '\n'.join(out)

    with open(inputfile, 'w') as f:
        f.write(output)

    return output


def predict_x(aenet,
              outputfile=None):
    """Run predict.x.

    outputfile is where the output will go.
    """
    if outputfile is None:
        outputfile = aenet['predict']['inputfile'] + '.out'

    cmd = [PREDICTX, aenet['predict']['inputfile'], '>', outputfile]
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    output, stderr = process.communicate()

    with open(outputfile, 'w') as f:
        f.write(output)

    aenet['predict']['outputfile'] = outputfile
    aenet_save(aenet)


def read_ascii_train_file(aenet, raw=True):
    """Parse the ascii train file into a dictionary.

    This contains information about the dataset, and each configuration in the
    dataset, including the fingerprints. The atomic numbers, positions and
    forces for each configuration is also in the output. There are no unit cells
    for the configurations though, that does not appear to be in the train file.

    Note there is a similar function in aenet.utils. This version is a little
    more automated and creates the ascii files as needed.

    Parameters
    ----------
    aenet : string
        the aenet state file.

    raw: boolean
        if True, use unscaled fingerprints.

    Returns
    -------
    a dictionary containing the data.
    {'header': {information about the computation},
     'configurations': [{information about each configuration},]}

    """
    trainfile = aenet['generate']['trainfile']
    if raw:
        rawopt = '--raw'
        rawext = '.raw'
    else:
        rawopt = ''
        rawext = ''

    asciifile = trainfile + '.ascii' + rawext

    if not os.path.exists(asciifile):
        cmd = ['trnset2ASCII.x', rawopt, trainfile, asciifile]
        process = Popen(cmd, stdout=PIPE, stderr=PIPE)
        process.communicate()

    f = open(asciifile, 'r')
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


def get_training_details(aenet):
    """Parse the training output file to get the convergence data.

    Parameters
    ----------
    aenet : dict
        The aenet state dictionary

    Returns
    -------
    epochs, train_mae, train_rmse, test_mae, test_rmse
    """
    outfile = aenet['train']['outputfile']
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
