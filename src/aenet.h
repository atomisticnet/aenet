#ifndef AENET_H_INCLUDED
#define AENET_H_INCLUDED

void aenet_init(int ntypes, char* atom_types[], int* stat);
void aenet_final(int* stat);
void aenet_print_info(void);
void aenet_load_potential(int type_id, char* filename, int* stat);
_Bool aenet_all_loaded(void);

double aenet_free_atom_energy(int type_id);

void aenet_atomic_energy(double coo_i[3], int type_i, int n_j,
                         double coo_j[], int type_j[], double* E_i,
                         int* stat);

void aenet_atomic_energy_and_forces(double coo_i[3], int type_i, int index_i,
                                    int n_j, double coo_j[], int type_j[],
                                    int index_j[], int natoms, double* E_i,
                                    double F[], int* stat);

void aenet_convert_atom_types(int ntypes_in, char* atom_types[],
                              int natoms_in, int type_id_in[],
                              int type_id_out[], int* stat);

void aenet_sfb_init(int ntypes, char* atom_types[], int radial_order,
                    int angular_order, double radial_Rc, double angular_Rc,
                    int* stat);

void aenet_sfb_final(int* stat);

int aenet_sfb_nvalues(void);

void aenet_sfb_eval(int itype0, double coo0[3], int nat, int itype1[],
                    double coo1[], int nv, double values[], int* stat);

void aenet_sfb_reconstruct_radial(int nv, double values[], int nx,
                                  double x[], double y[], int* stat);

extern int AENET_OK;
extern int AENET_ERR_INIT;
extern int AENET_ERR_MALLOC;
extern int AENET_ERR_IO;
extern int AENET_ERR_TYPE;
extern int AENET_TYPELEN;
extern int AENET_PATHLEN;

extern int aenet_nsf_max;
extern int aenet_nnb_max;
extern double aenet_Rc_min;
extern double aenet_Rc_max;

#endif
