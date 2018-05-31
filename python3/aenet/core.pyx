"""
Cython interface to aenetLib
"""

__author__ = "Alexander Urban, Nongnuch Artrith"
__date__   = "2014-09-03"

cimport aenet.core as core
from aenet.timing import Tick
ticker=Tick()

from libc.stdlib cimport malloc, free
from cython.view cimport array as cvarray
from cpython cimport bool

cimport numpy as np
import numpy as np

import sys

ctypedef unsigned char char_type

cdef char_type[:] _chars(s):
    if isinstance(s, unicode):
        # encode to the specific encoding used inside of the module
        s = (<unicode>s).encode('utf8')
    return s

class AenetError(Exception):
    pass

class InitError(Exception):
    pass

class SpeciesError(Exception):
    pass

cdef class ANNPotentials:
    """
    Artificial Neural Network (ANN) Potential
    """

    cdef:
        int ntypes

    cdef readonly:
        list atom_types
        dict potential_files
        double Rc_min
        double Rc_max
        int nnb_max

    def __init__(self, dict potential_files):
        cdef:
            char** cstring_array = <char**>malloc(len(potential_files)*sizeof(char*))
            int stat
            str filename
            str typename

        self.atom_types = list(potential_files.keys())
        self.potential_files = potential_files
        self.ntypes = len(self.atom_types)
        self.Rc_min = 0.0
        self.Rc_max = 0.0

        bytes_types = []
        bytes_files = []
        for i in range(self.ntypes):
            if isinstance(self.atom_types[i], unicode):
                bytes_types.append((<unicode>self.atom_types[i]).encode('UTF-8'))
            else:
                bytes_types.append(self.atom_types[i].encode('UTF-8'))
            if isinstance(potential_files[self.atom_types[i]], unicode):
                bytes_files.append((<unicode>potential_files[self.atom_types[i]]
                                ).encode('UTF-8'))
            else:
                bytes_files.append(potential_files[self.atom_types[i]
                                               ].encode('UTF-8'))
            cstring_array[i] = bytes_types[i]
        self.atom_types = [x.decode('utf-8') for x in bytes_types]

        # initialize aenetLib
        core.aenet_init(self.ntypes, cstring_array, &stat)
        self._check_aenet_return_status(stat)

        # load potentials
        for i in range(self.ntypes):
            typename = self.atom_types[i]
            filename = bytes_files[i].decode('utf-8')
            self._load_potential(typename, filename)
        assert(core.aenet_all_loaded() == True)

        self.Rc_min = core.aenet_Rc_min
        self.Rc_max = core.aenet_Rc_max
        self.nnb_max = core.aenet_nnb_max

    def __dealloc__(self):
        cdef int stat
        core.aenet_final(&stat)
        if stat == core.AENET_ERR_MALLOC:
            sys.stderr.write("Warning: memory deallocation failed (aenetLib).")

    def __str__(self):
        out  = " ANN Potentials\n"
        out += "   Atomic species: "
        out += (self.ntypes*"%s ") % tuple(self.atom_types)
        out += "\n"
        return out

    def _load_potential(self, str atom_type, str filename):
        cdef:
            char* cfilename
            int type_id
            int stat
        bytestring = filename.encode('UTF-8')
        cfilename = bytestring
        type_id = self.get_type_id(atom_type)

        core.aenet_load_potential(type_id, cfilename, &stat)
        if stat == core.AENET_ERR_TYPE:
            raise SpeciesError("Unknown type ID: %d" % (type_id,))
        elif stat == core.AENET_ERR_IO:
            raise IOError("File not found: %s" % (filename,))
        self._check_aenet_return_status(stat)

    def print_info_to_stdout(self):
        core.aenet_print_info()

    cpdef double atomic_energy(self, double[:] coo_i, str type_i,
                               double[:,::1] coo_j, list type_j) except? -1:
        cdef:
            int itype_i
            int n_j = len(type_j)
            int[:] itype_j = cvarray(shape=(n_j,), itemsize=sizeof(int), format="i")
            int stat
            double E_i
            str t
            int i

        itype_i = self.get_type_id(type_i)
        for i,t in enumerate(type_j):
            itype_j[i] = self.get_type_id(str(t))
        core.aenet_atomic_energy(&coo_i[0], itype_i, n_j,
                                 &coo_j[0][0], &itype_j[0],
                                 &E_i, &stat)
        self._check_aenet_return_status(stat)

        return E_i

    cpdef double atomic_energy_and_forces(self, double[:] coo_i, str type_i,
                                          int index_i, double[:,::1] coo_j,
                                          list type_j, int[:] index_j,
                                          double[:,::1] F) except? -1:
        cdef:
            int itype_i
            int n_j = len(type_j)
            int natoms = F.shape[0]
            int[:] itype_j = cvarray(shape=(n_j,), itemsize=sizeof(int),
                                     format="i")
            int stat
            double E_i
            str t
            int i

        itype_i = self.get_type_id(type_i)
        for i,t in enumerate(type_j):
            itype_j[i] = self.get_type_id(str(t))
        core.aenet_atomic_energy_and_forces(&coo_i[0], itype_i, index_i,
                                            n_j, &coo_j[0][0], &itype_j[0],
                                            &index_j[0], natoms, &E_i,
                                            &F[0][0], &stat)
        self._check_aenet_return_status(stat)

        return E_i

    cpdef int get_type_id(self, str atom_type) except? -1:
        cdef:
            int i
            str t
            int type_id
        type_id = 0
        for (i,t) in enumerate(self.atom_types):
            if t.strip() == atom_type.strip():
                type_id = i + 1
                break
        if type_id == 0:
            raise SpeciesError("Unknown atomic species: {}".format(atom_type))
        return type_id

    cdef int _check_aenet_return_status(self, stat) except? -1:
        """Check the return status and raise corresponding exceptions."""
        if stat == core.AENET_ERR_MALLOC:
            raise MemoryError("aenetLib: Memory allocation failed.")
        elif stat == core.AENET_ERR_INIT:
            raise InitError("aenetLib: Module not initialized "
                            "(or attempt to re-initialize module).")
        elif stat == core.AENET_ERR_IO:
            raise IOError("aenetLib: I/O error encountered.")
        elif stat == core.AENET_ERR_TYPE:
            raise SpeciesError("aenetLib: Invalid atomic species.")
        elif stat != 0:
            raise AenetError("aenetLib: Unspecified error encountered.")
