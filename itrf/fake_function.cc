/*
 * fake the missing functions
 */

#define REALTYPE double

typedef struct pdata{
    REALTYPE F[21];
    REALTYPE U[6][3];
    REALTYPE twozeta_a;
    REALTYPE twozeta_b;
    REALTYPE twozeta_c;
    REALTYPE twozeta_d;
    REALTYPE oo2z;
    REALTYPE oo2n;
    REALTYPE oo2zn;
    REALTYPE poz;
    REALTYPE pon;
    REALTYPE oo2p;
    REALTYPE ss_r12_ss;
} prim_data;

typedef struct {
    REALTYPE *int_stack;
    prim_data *PrimQuartet;
    REALTYPE AB[3];
    REALTYPE CD[3];
    REALTYPE *vrr_classes[11][11];
    REALTYPE *vrr_stack;
} Libint_t;

typedef struct {
    double *int_stack;
    prim_data *PrimQuartet;
    double *zero_stack;
    double *ABCD[12+144];
    double AB[3];
    double CD[3];
    double *deriv_classes[9][9][12];
    double *deriv2_classes[9][9][144];
    double *dvrr_classes[9][9];
    double *dvrr_stack;
} Libderiv_t;


extern "C" {
double *(*build_eri[1][1][1][1])(Libint_t *, int) = {0};
void (*build_deriv1_eri[1][1][1][1])(Libderiv_t *, int) = {0};
void (*build_deriv12_eri[1][1][1][1])(Libderiv_t *, int) = {0};

int init_libint(Libint_t *libint, int max_am, int max_num_prim_quartets) {}
int init_libint_base() {}
void erd__memory_eri_batch_() {}
void erd__gener_eri_batch_() {}
int free_libint(Libint_t *libint) {}

int init_libderiv1(Libderiv_t *libderiv, int a1, int a2, int a3) {}
int init_libderiv12(Libderiv_t *libderiv, int a1, int a2, int a3) {}
int init_libderiv_base() {}
int free_libderiv(Libderiv_t *libderiv) {}
}

