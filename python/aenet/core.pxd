cdef extern from 'aenet.h':
   void aenet_init(int ntypes, char* atom_types[], int* stat)
   void aenet_final(int* stat)
   void aenet_print_info()
   void aenet_load_potential(int type_id, char* filename, int* stat)
   bint aenet_all_loaded()

   double aenet_free_atom_energy(int type_id)

   void aenet_atomic_energy(double coo_i[3], int type_i, int n_j,
                            double coo_j[], int type_j[], double* E_i,
                            int* stat)

   void aenet_atomic_energy_and_forces(double coo_i[3], int type_i, int index_i,
                                       int n_j, double coo_j[], int type_j[],
                                       int index_j[], int natoms, double* E_i,
                                       double F[], int* stat)

   void aenet_convert_atom_types(int ntypes_in, char* atom_types[],
                                 int natoms_in, int type_id_in[],
                                 int type_id_out[], int* stat)
   int AENET_OK
   int AENET_ERR_INIT
   int AENET_ERR_MALLOC
   int AENET_ERR_IO
   int AENET_ERR_TYPE
   int AENET_TYPELEN
   int AENET_PATHLEN

   int aenet_nsf_max
   int aenet_nnb_max
   double aenet_Rc_min
   double aenet_Rc_max
