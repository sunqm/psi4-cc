/*
 * read_den_mat.cc
 */

#include <iostream>
#include <assert.h>
#include <libparallel/parallel.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libciomr/libciomr.h>
#define EXTERN
#include "psi4itrf.h"

#define LET(VAR,X0,X1,X2,X3)    VAR[0] = X0; \
                                VAR[1] = X1; \
                                VAR[2] = X2; \
                                VAR[3] = X3;

using namespace psi;

namespace psi { namespace ccdensity {
    void init_io();
    void get_moinfo_light();
    void get_params(Options& options);
    void get_rho_params(Options& options);
    int **cacheprep_rhf(int level, int *cachefiles);
    void cachedone_rhf(int **cachelist);
    void cleanup_light();
    void exit_io();
    }
    namespace ccenergy {
        int **cacheprep_rhf(int level, int *cachefiles);
        int **cacheprep_uhf(int level, int *cachefiles);
    }
}
void init_dpd_rhf(int nmo, std::vector<int *> spaces,
                  int *cachefiles, int **cachelist);
void init_dpd_uhf(int nmo, std::vector<int *> spaces,
                  int *cachefiles, int **cachelist);

static void copy_rdm1(double *rdm1, const int nmo, int start[2], int count[2],
                      const int pnum, const int qnum, const char lbl[])
{
    const int irrep = 0;
    int p, q, p0, q0, pq;
    dpdfile2 D;
    global_dpd_->file2_init(&D, PSIF_CC_OEI, irrep, pnum, qnum, lbl);
    global_dpd_->file2_mat_init(&D);
    global_dpd_->file2_mat_rd(&D);

    for (p = start[0], p0 = 0; p0 < count[0]; p++, p0++) {
    for (q = start[1], q0 = 0; q0 < count[1]; q++, q0++) {
        rdm1[p*nmo+q] = D.matrix[0][p0][q0];
    } }

    global_dpd_->file2_mat_close(&D);
    global_dpd_->file2_close(&D);
}

static void copy_rdm2(double *rdm2, int nmo, int start[4], int count[4],
                      int pqnum, int rsnum, int fpqnum, int frsnum,
                      const char lbl[])
{
    int p, q, r, s, pq, rs;
    dpdbuf4 G;
    global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, pqnum, rsnum, fpqnum, frsnum, 0, lbl);
    global_dpd_->buf4_mat_irrep_init(&G, 0);
    global_dpd_->buf4_mat_irrep_rd(&G, 0);

    pq = 0;
    for (p = start[0]; p < start[0]+count[0]; p++) {
    for (q = start[1]; q < start[1]+count[1]; q++) {
        rs = 0;
        for (r = start[2]; r < start[2]+count[2]; r++) {
        for (s = start[3]; s < start[3]+count[3]; s++) {
            rdm2[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s] = G.matrix[0][pq][rs];
            rs++;
        } }
        pq++;
    } }

    global_dpd_->buf4_mat_irrep_close(&G, 0);
    global_dpd_->buf4_close(&G);
}

// transpose a block
static void rdm1_dagger(double *rdm1, const int nmo,
                        int start[2], int count[2])
{
    int p, q, pi, po;

    for (p = start[0]; p < start[0]+count[0]; p++) {
    for (q = start[1]; q < start[1]+count[1]; q++) {
        pi = p * nmo + q;
        po = q * nmo + p;
        rdm1[po] = rdm1[pi];
    } }
}

static void rdm2_swap_e1e2(double *rdm2, const int nmo,
                           int start[4], int count[4])
{
    int p, q, r, s, pi, po;

    for (p = start[0]; p < start[0]+count[0]; p++) {
    for (q = start[1]; q < start[1]+count[1]; q++) {
    for (r = start[2]; r < start[2]+count[2]; r++) {
    for (s = start[3]; s < start[3]+count[3]; s++) {
        pi = p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s;
        po = q*nmo*nmo*nmo+p*nmo*nmo+s*nmo+r;
        rdm2[po] = rdm2[pi];
    } } } }
}
static void rdm2_copyswap_e1e2(double *rdm2i, double *rdm2o,
                               const int nmo, int start[4], int count[4])
{
    int p, q, r, s, pi, po;

    for (p = start[0]; p < start[0]+count[0]; p++) {
    for (q = start[1]; q < start[1]+count[1]; q++) {
    for (r = start[2]; r < start[2]+count[2]; r++) {
    for (s = start[3]; s < start[3]+count[3]; s++) {
        pi = p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s;
        po = q*nmo*nmo*nmo+p*nmo*nmo+s*nmo+r;
        rdm2o[po] = rdm2i[pi];
    } } } }
}


static void rdm2_dagger(double *rdm2, const int nmo,
                        int start[4], int count[4])
{
    int p, q, r, s, pi, po;

    for (p = start[0]; p < start[0]+count[0]; p++) {
    for (q = start[1]; q < start[1]+count[1]; q++) {
    for (r = start[2]; r < start[2]+count[2]; r++) {
    for (s = start[3]; s < start[3]+count[3]; s++) {
        pi = p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s;
        po = r*nmo*nmo*nmo+s*nmo*nmo+p*nmo+q;
        rdm2[po] = rdm2[pi];
    } } } }
}

// rdm2 is anti-symmetric
/* there are three kinds of density matrix defined in psi4.
 * they are normal CC-density, Fock-adjusted density, Mulliken density
 * this one based on Mulliken density, see comments in deanti.cc */
void read_cc_dm_rhf(double *rdm1, double *rdm2, int nmo)
{
    using namespace ccdensity;
    int nocc, nvir;
    int start[4], count[4];

    init_io();
    chkpt_init(PSIO_OPEN_OLD);
    std::vector<int*> spaces;
    int *cachefiles = init_int_array(PSIO_MAXUNIT);
    int **cachelist = ccenergy::cacheprep_rhf(CACHELEVEL, cachefiles);
    init_dpd_rhf(nmo, spaces, cachefiles, cachelist);
    int nirreps = chkpt_rd_nirreps();
    int *occpi = chkpt_rd_clsdpi();
    nocc = 0;
    for (int i = 0; i < nirreps; i++) {
        nocc += occpi[i];
    }
    free(occpi);

    nvir = nmo - nocc;

    start[0] = 0; count[0] = nocc;
    start[1] = 0; count[1] = nocc;
    copy_rdm1(rdm1, nmo, start, count, 0, 0, "DIJ -1");

    start[0] = nocc; count[0] = nvir;
    start[1] = nocc; count[1] = nvir;
    copy_rdm1(rdm1, nmo, start, count, 1, 1, "DAB -1");

    start[0] = 0   ; count[0] = nocc;
    start[1] = nocc; count[1] = nvir;
    copy_rdm1(rdm1, nmo, start, count, 0, 1, "DAI -1");
    rdm1_dagger(rdm1, nmo, start, count);

    start[0] = 0   ; count[0] = nocc;
    start[1] = nocc; count[1] = nvir;
    copy_rdm1(rdm1, nmo, start, count, 0, 1, "DIA -1");

    LET(start, 0   , 0   , 0   , 0   );
    LET(count, nocc, nocc, nocc, nocc);
    copy_rdm2(rdm2, nmo, start, count, 0, 0, 0, 0, "2 Gijkl - Gijlk");
    
    LET(start, 0   , 0   , 0   , nocc);
    LET(count, nocc, nocc, nocc, nvir);
    copy_rdm2(rdm2, nmo, start, count, 0, 10, 0, 10, "2 Gijka - Gjika");
    rdm2_swap_e1e2(rdm2, nmo, start, count); // oovo
    rdm2_dagger(rdm2, nmo, start, count); // ovoo
    LET(start, 0   , 0   , nocc, 0   );
    LET(count, nocc, nocc, nvir, nocc);
    rdm2_dagger(rdm2, nmo, start, count); // vooo
    
#if defined ccdensity_deanti
    // we cannot use GIjAb because it is ruined by deanti_RHF in PSI4.0b5
    dpdbuf4 G1, G2, G3;
    int i, j, a, b, pq, rs;
    global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
    global_dpd_->buf4_mat_irrep_init(&G1, 0);
    global_dpd_->buf4_mat_irrep_rd(&G1, 0);
    global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
    global_dpd_->buf4_mat_irrep_init(&G2, 0);
    global_dpd_->buf4_mat_irrep_rd(&G2, 0);
    global_dpd_->buf4_init(&G3, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "2 Gijab - Gijba");
    global_dpd_->buf4_mat_irrep_init(&G3, 0);
    global_dpd_->buf4_mat_irrep_rd(&G3, 0);
    for (i = 0, pq = 0; i < nocc; i++) {
    for (j = 0; j < nocc; j++, pq++) {
    for (a = nocc, rs = 0; a < nmo; a++) {
    for (b = nocc; b < nmo; b++, rs++) {
        rdm2[i*nmo*nmo*nmo+j*nmo*nmo+a*nmo+b] = G3.matrix[0][pq][rs] \
                + 0*G1.matrix[0][i*nvir+b-nocc][j*nvir+a-nocc] \
                + 0*G2.matrix[0][i*nvir+b-nocc][j*nvir+a-nocc];
    } } } } // oovv fixme due to the change in fold.cc
    LET(start, 0   , 0   , nocc, nocc);
    LET(count, nocc, nocc, nvir, nvir);
    rdm2_dagger(rdm2, nmo, start, count); // vvoo

    for (i = 0, pq = 0; i < nocc; i++) {
    for (b = nocc; b < nmo; b++, pq++) {
    for (j = 0, rs = 0; j < nocc; j++) {
    for (a = nocc; a < nmo; a++, rs++) {
        rdm2[i*nmo*nmo*nmo+b*nmo*nmo+a*nmo+j] =-G1.matrix[0][pq][rs] \
                                              - G2.matrix[0][pq][rs];
    } } } } // ovvo
    global_dpd_->buf4_mat_irrep_close(&G1, 0);
    global_dpd_->buf4_mat_irrep_close(&G2, 0);
    global_dpd_->buf4_mat_irrep_close(&G3, 0);
    global_dpd_->buf4_close(&G1);
    global_dpd_->buf4_close(&G2);
    global_dpd_->buf4_close(&G3);
    
    LET(start, 0   , nocc, nocc, 0   );
    LET(count, nocc, nvir, nvir, nocc);
    rdm2_dagger(rdm2, nmo, start, count); // voov
    
    LET(start, 0   , nocc, 0   , nocc);
    LET(count, nocc, nvir, nocc, nvir);
    copy_rdm2(rdm2, nmo, start, count, 10, 10, 10, 10, "GIbJa");
    rdm2_swap_e1e2(rdm2, nmo, start, count); // vovo
#else
    dpdbuf4 G1, G2;
    int i, j, a, b, pq, rs;
    global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
    global_dpd_->buf4_mat_irrep_init(&G1, 0);
    global_dpd_->buf4_mat_irrep_rd(&G1, 0);
    for (i = 0, pq = 0; i < nocc; i++) {
    for (j = 0; j < nocc; j++, pq++) {
    for (a = nocc, rs = 0; a < nmo; a++) {
    for (b = nocc; b < nmo; b++, rs++) {
        rdm2[i*nmo*nmo*nmo+j*nmo*nmo+a*nmo+b] = \
                + G1.matrix[0][pq][rs] * 2 \
                - G1.matrix[0][pq][(b-nocc)*nvir+a-nocc];
    } } } } // oovv
    global_dpd_->buf4_mat_irrep_close(&G1, 0);
    global_dpd_->buf4_close(&G1);
    LET(start, 0   , 0   , nocc, nocc);
    LET(count, nocc, nocc, nvir, nvir);
    rdm2_dagger(rdm2, nmo, start, count); // vvoo

    global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
    global_dpd_->buf4_mat_irrep_init(&G1, 0);
    global_dpd_->buf4_mat_irrep_rd(&G1, 0);
    global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
    global_dpd_->buf4_mat_irrep_init(&G2, 0);
    global_dpd_->buf4_mat_irrep_rd(&G2, 0);
    for (i = 0, pq = 0; i < nocc; i++) {
    for (b = nocc; b < nmo; b++, pq++) {
    for (j = 0, rs = 0; j < nocc; j++) {
    for (a = nocc; a < nmo; a++, rs++) {
        rdm2[i*nmo*nmo*nmo+b*nmo*nmo+a*nmo+j] =-G1.matrix[0][pq][rs] \
                                              - G2.matrix[0][pq][rs];
    } } } } // ovvo
    global_dpd_->buf4_mat_irrep_close(&G1, 0);
    global_dpd_->buf4_mat_irrep_close(&G2, 0);
    global_dpd_->buf4_close(&G1);
    global_dpd_->buf4_close(&G2);
    LET(start, 0   , nocc, nocc, 0   );
    LET(count, nocc, nvir, nvir, nocc);
    rdm2_dagger(rdm2, nmo, start, count); // voov
    
    global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
    global_dpd_->buf4_mat_irrep_init(&G1, 0);
    global_dpd_->buf4_mat_irrep_rd(&G1, 0);
    global_dpd_->buf4_init(&G2, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
    global_dpd_->buf4_mat_irrep_init(&G2, 0);
    global_dpd_->buf4_mat_irrep_rd(&G2, 0);
    for (i = 0, pq = 0; i < nocc; i++) {
    for (b = nocc; b < nmo; b++, pq++) {
    for (j = 0, rs = 0; j < nocc; j++) {
    for (a = nocc; a < nmo; a++, rs++) {
        rdm2[i*nmo*nmo*nmo+b*nmo*nmo+j*nmo+a] = G1.matrix[0][pq][rs] \
                                              + G2.matrix[0][pq][rs];
    } } } } // ovvo
    global_dpd_->buf4_mat_irrep_close(&G1, 0);
    global_dpd_->buf4_mat_irrep_close(&G2, 0);
    global_dpd_->buf4_close(&G1);
    global_dpd_->buf4_close(&G2);
    LET(start, 0   , nocc, 0   , nocc);
    LET(count, nocc, nvir, nocc, nvir);
    rdm2_swap_e1e2(rdm2, nmo, start, count); // vovo
#endif
    
    LET(start, nocc, 0   , nocc, nocc);
    LET(count, nvir, nocc, nvir, nvir);
    copy_rdm2(rdm2, nmo, start, count, 11, 5, 11, 5, "2 Gciab - Gciba");
    rdm2_swap_e1e2(rdm2, nmo, start, count); // ovvv
    rdm2_dagger(rdm2, nmo, start, count); // vvvo
    LET(start, 0   , nocc, nocc, nocc);
    LET(count, nocc, nvir, nvir, nvir);
    rdm2_dagger(rdm2, nmo, start, count); // vvov

    LET(start, nocc, nocc, nocc, nocc);
    LET(count, nvir, nvir, nvir, nvir);
    copy_rdm2(rdm2, nmo, start, count, 5, 5, 5, 5, "2 Gabcd - Gabdc");

    for (i = 0; i < nmo*nmo*nmo*nmo; i++) {
        rdm2[i] *= 2;
    }

    dpd_close(0);
    chkpt_close();
    exit_io();
    for (int i = 0; i < spaces.size(); i++) {
        free(spaces[i]);
    }
    free(cachefiles);
    free_int_matrix(cachelist);
}

/*
 * rdm2[ia,bj] = -rdm2[ia,jb]
 */
static void copy_uhf_iabj(double *rdm2, const int nmo,
                          int start[4], int count[4])
{
    int p, q, r, s, pi, po;

    for (p = start[0]; p < start[0]+count[0]; p++) {
    for (q = start[1]; q < start[1]+count[1]; q++) {
    for (r = start[2]; r < start[2]+count[2]; r++) {
    for (s = start[3]; s < start[3]+count[3]; s++) {
        pi = p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s;
        po = p*nmo*nmo*nmo+q*nmo*nmo+s*nmo+r;
        rdm2[po] = -rdm2[pi];
    } } } }
}
/*
 * rdm2_ab[pq,rs]: [p=e1_alpha,q=e2_beta; r=e1_alpha,s=e2_beta]
 * rdm2_ba[pq,rs]: [p=e1_beta,q=e2_alpha; r=e1_beta,s=e2_alpha]
 */
void read_cc_dm_uhf(double *rdm1_a, double *rdm1_b,
                    double *rdm2_aa, double *rdm2_bb,
                    double *rdm2_ab, double *rdm2_ba, int nmo)
{
    using namespace ccdensity;
    int start[4], count[4];

    init_io();
    chkpt_init(PSIO_OPEN_OLD);
    std::vector<int*> spaces;
    int *cachefiles = init_int_array(PSIO_MAXUNIT);
    int **cachelist = ccenergy::cacheprep_uhf(CACHELEVEL, cachefiles);
    init_dpd_uhf(nmo, spaces, cachefiles, cachelist);
    int *clsdpi = chkpt_rd_clsdpi();
    int *openpi = chkpt_rd_openpi();
    int nocc_a = clsdpi[0] + openpi[0];
    int nocc_b = clsdpi[0];
    int nvir_a = nmo - nocc_a;
    int nvir_b = nmo - nocc_b;
    free(clsdpi);
    free(openpi);

    start[0] = 0; count[0] = nocc_a;
    start[1] = 0; count[1] = nocc_a;
    copy_rdm1(rdm1_a, nmo, start, count, 0, 0, "DIJ -1");
    start[0] = 0; count[0] = nocc_b;
    start[1] = 0; count[1] = nocc_b;
    copy_rdm1(rdm1_b, nmo, start, count, 2, 2, "Dij -1");

    start[0] = nocc_a; count[0] = nvir_a;
    start[1] = nocc_a; count[1] = nvir_a;
    copy_rdm1(rdm1_a, nmo, start, count, 1, 1, "DAB -1");
    start[0] = nocc_b; count[0] = nvir_b;
    start[1] = nocc_b; count[1] = nvir_b;
    copy_rdm1(rdm1_b, nmo, start, count, 3, 3, "Dab -1");

    start[0] = 0     ; count[0] = nocc_a;
    start[1] = nocc_a; count[1] = nvir_a;
    copy_rdm1(rdm1_a, nmo, start, count, 0, 1, "DAI -1");
    rdm1_dagger(rdm1_a, nmo, start, count);
    start[0] = 0     ; count[0] = nocc_b;
    start[1] = nocc_b; count[1] = nvir_b;
    copy_rdm1(rdm1_b, nmo, start, count, 2, 3, "Dai -1");
    rdm1_dagger(rdm1_b, nmo, start, count);

    start[0] = 0     ; count[0] = nocc_a;
    start[1] = nocc_a; count[1] = nvir_a;
    copy_rdm1(rdm1_a, nmo, start, count, 0, 1, "DIA -1");
    start[0] = 0     ; count[0] = nocc_b;
    start[1] = nocc_b; count[1] = nvir_b;
    copy_rdm1(rdm1_b, nmo, start, count, 2, 3, "Dia -1");

    LET(start, 0     , 0     , 0     , 0     );
    LET(count, nocc_a, nocc_a, nocc_a, nocc_a);
    copy_rdm2(rdm2_aa, nmo, start, count, 0, 0, 2, 2, "GIJKL");
    LET(start, 0     , 0     , 0     , 0     );
    LET(count, nocc_b, nocc_b, nocc_b, nocc_b);
    copy_rdm2(rdm2_bb, nmo, start, count, 10, 10, 12, 12, "Gijkl");
    LET(start, 0     , 0     , 0     , 0     );
    LET(count, nocc_a, nocc_b, nocc_a, nocc_b);
    copy_rdm2(rdm2_ab, nmo, start, count, 22, 22, 22, 22, "GIjKl");
    rdm2_copyswap_e1e2(rdm2_ab, rdm2_ba, nmo, start, count);
    
    LET(start, 0     , 0     , 0     , nocc_a);
    LET(count, nocc_a, nocc_a, nocc_a, nvir_a);
    copy_rdm2(rdm2_aa, nmo, start, count, 0, 20, 2, 20, "GIJKA");
    rdm2_swap_e1e2(rdm2_aa, nmo, start, count); // OOVO
    rdm2_dagger(rdm2_aa, nmo, start, count); // OVOO
    LET(start, 0     , 0     , nocc_a, 0     );
    LET(count, nocc_a, nocc_a, nvir_a, nocc_a);
    rdm2_dagger(rdm2_aa, nmo, start, count); // VOOO
    LET(start, 0     , 0     , 0     , nocc_b);
    LET(count, nocc_b, nocc_b, nocc_b, nvir_b);
    copy_rdm2(rdm2_bb, nmo, start, count, 10, 30, 12, 30, "Gijka");
    rdm2_swap_e1e2(rdm2_bb, nmo, start, count); // oovo
    rdm2_dagger(rdm2_bb, nmo, start, count); // ovoo
    LET(start, 0     , 0     , nocc_b, 0     );
    LET(count, nocc_b, nocc_b, nvir_b, nocc_b);
    rdm2_dagger(rdm2_bb, nmo, start, count); // vooo

    LET(start, 0     , 0     , 0     , nocc_b);
    LET(count, nocc_a, nocc_b, nocc_a, nvir_b);
    copy_rdm2(rdm2_ab, nmo, start, count, 22, 24, 22, 24, "GIjKa");
    rdm2_copyswap_e1e2(rdm2_ab, rdm2_ba, nmo, start, count); // oOvO
    rdm2_dagger(rdm2_ab, nmo, start, count); // OvOo
    LET(start, 0     , nocc_b, 0     , 0     );
    LET(count, nocc_a, nvir_b, nocc_a, nocc_b);
    rdm2_copyswap_e1e2(rdm2_ab, rdm2_ba, nmo, start, count); // vOoO
    LET(start, 0     , 0     , 0     , nocc_a);
    LET(count, nocc_b, nocc_a, nocc_b, nvir_a);
    copy_rdm2(rdm2_ba, nmo, start, count, 23, 27, 23, 27, "GiJkA");
    rdm2_copyswap_e1e2(rdm2_ba, rdm2_ab, nmo, start, count); // OoVo
    rdm2_dagger(rdm2_ba, nmo, start, count); // oVoO
    LET(start, 0     , nocc_a, 0     , 0     );
    LET(count, nocc_b, nvir_a, nocc_b, nocc_a);
    rdm2_copyswap_e1e2(rdm2_ba, rdm2_ab, nmo, start, count); // VoOo

    LET(start, 0     , 0     , nocc_a, nocc_a);
    LET(count, nocc_a, nocc_a, nvir_a, nvir_a);
    copy_rdm2(rdm2_aa, nmo, start, count, 0, 5, 2, 7, "GIJAB");
    rdm2_dagger(rdm2_aa, nmo, start, count);
    LET(start, 0     , 0     , nocc_b, nocc_b);
    LET(count, nocc_b, nocc_b, nvir_b, nvir_b);
    copy_rdm2(rdm2_bb, nmo, start, count, 10, 15, 12, 17, "Gijab");
    rdm2_dagger(rdm2_bb, nmo, start, count);

    LET(start, 0     , 0     , nocc_a, nocc_b);
    LET(count, nocc_a, nocc_b, nvir_a, nvir_b);
    copy_rdm2(rdm2_ab, nmo, start, count, 22, 28, 22, 28, "GIjAb");
    rdm2_dagger(rdm2_ab, nmo, start, count); // VvOo
    rdm2_copyswap_e1e2(rdm2_ab, rdm2_ba, nmo, start, count); // oOvV
    LET(start, 0     , 0     , nocc_b, nocc_a);
    LET(count, nocc_b, nocc_a, nvir_b, nvir_a);
    rdm2_dagger(rdm2_ba, nmo, start, count); // vVoO

    // G[IbAj] * <Ib|Aj> = -G[IbjA] * <Ij|Ab>
    dpdbuf4 G1;
    int i, j, a, b, pq, rs;
    double val;
    global_dpd_->buf4_init(&G1, PSIF_CC_GAMMA, 0, 24, 27, 24, 27, 0, "GIbjA");
    global_dpd_->buf4_mat_irrep_init(&G1, 0);
    global_dpd_->buf4_mat_irrep_rd(&G1, 0);
    for (i = 0, pq = 0; i < nocc_a; i++) {
    for (b = nocc_b; b < nmo; b++, pq++) {
    for (j = 0, rs = 0; j < nocc_b; j++) {
    for (a = nocc_a; a < nmo; a++, rs++) {
        val = G1.matrix[0][pq][rs];
        rdm2_ab[i*nmo*nmo*nmo+b*nmo*nmo+a*nmo+j] =-val;
        rdm2_ab[a*nmo*nmo*nmo+j*nmo*nmo+i*nmo+b] =-val;
        rdm2_ba[b*nmo*nmo*nmo+i*nmo*nmo+j*nmo+a] =-val;
        rdm2_ba[j*nmo*nmo*nmo+a*nmo*nmo+b*nmo+i] =-val;
    } } } }
    global_dpd_->buf4_mat_irrep_close(&G1, 0);
    global_dpd_->buf4_close(&G1);
    
    LET(start, 0     , nocc_a, 0     , nocc_a);
    LET(count, nocc_a, nvir_a, nocc_a, nvir_a);
    copy_rdm2(rdm2_aa, nmo, start, count, 20, 20, 20, 20, "GIBJA");
    rdm2_swap_e1e2(rdm2_aa, nmo, start, count); // vovo
    copy_uhf_iabj(rdm2_aa, nmo, start, count); // ovvo
    LET(start, 0     , nocc_a, nocc_a, 0     );
    LET(count, nocc_a, nvir_a, nvir_a, nocc_a);
    rdm2_dagger(rdm2_aa, nmo, start, count); // voov

    LET(start, 0     , nocc_b, 0     , nocc_b);
    LET(count, nocc_b, nvir_b, nocc_b, nvir_b);
    copy_rdm2(rdm2_bb, nmo, start, count, 30, 30, 30, 30, "Gibja");
    rdm2_swap_e1e2(rdm2_bb, nmo, start, count);
    copy_uhf_iabj(rdm2_bb, nmo, start, count);
    LET(start, 0     , nocc_b, nocc_b, 0     );
    LET(count, nocc_b, nvir_b, nvir_b, nocc_b);
    rdm2_dagger(rdm2_bb, nmo, start, count);

    LET(start, 0     , nocc_b, 0     , nocc_b);
    LET(count, nocc_a, nvir_b, nocc_a, nvir_b);
    copy_rdm2(rdm2_ab, nmo, start, count, 24, 24, 24, 24, "GIbJa");
    rdm2_copyswap_e1e2(rdm2_ab, rdm2_ba, nmo, start, count);

    LET(start, 0     , nocc_a, 0     , nocc_a);
    LET(count, nocc_b, nvir_a, nocc_b, nvir_a);
    copy_rdm2(rdm2_ba, nmo, start, count, 27, 27, 27, 27, "GiBjA");
    rdm2_copyswap_e1e2(rdm2_ba, rdm2_ab, nmo, start, count);

    LET(start, nocc_a, 0     , nocc_a, nocc_a);
    LET(count, nvir_a, nocc_a, nvir_a, nvir_a);
    copy_rdm2(rdm2_aa, nmo, start, count, 21, 5, 21, 7, "GCIAB");
    rdm2_swap_e1e2(rdm2_aa, nmo, start, count); // OVVV
    rdm2_dagger(rdm2_aa, nmo, start, count); // VVVO
    LET(start, 0     , nocc_a, nocc_a, nocc_a);
    LET(count, nocc_a, nvir_a, nvir_a, nvir_a);
    rdm2_dagger(rdm2_aa, nmo, start, count); // VVOV
    LET(start, nocc_b, 0     , nocc_b, nocc_b);
    LET(count, nvir_b, nocc_b, nvir_b, nvir_b);
    copy_rdm2(rdm2_bb, nmo, start, count, 31, 15, 31, 17, "Gciab");
    rdm2_swap_e1e2(rdm2_bb, nmo, start, count); // ovvv
    rdm2_dagger(rdm2_bb, nmo, start, count); // vvvo
    LET(start, 0     , nocc_b, nocc_b, nocc_b);
    LET(count, nocc_b, nvir_b, nvir_b, nvir_b);
    rdm2_dagger(rdm2_bb, nmo, start, count); // vvov

    LET(start, nocc_a, 0     , nocc_a, nocc_b);
    LET(count, nvir_a, nocc_b, nvir_a, nvir_b);
    copy_rdm2(rdm2_ab, nmo, start, count, 26, 28, 26, 28, "GCiAb");
    rdm2_copyswap_e1e2(rdm2_ab, rdm2_ba, nmo, start, count); // oVvV
    rdm2_dagger(rdm2_ab, nmo, start, count); // VvVo
    LET(start, nocc_a, nocc_b, nocc_a, 0     );
    LET(count, nvir_a, nvir_b, nvir_a, nocc_b);
    rdm2_copyswap_e1e2(rdm2_ab, rdm2_ba, nmo, start, count); // vVoV
    LET(start, nocc_b, 0     , nocc_b, nocc_a);
    LET(count, nvir_b, nocc_a, nvir_b, nvir_a);
    copy_rdm2(rdm2_ba, nmo, start, count, 25, 29, 25, 29, "GcIaB");
    rdm2_copyswap_e1e2(rdm2_ba, rdm2_ab, nmo, start, count); // OvVv
    rdm2_dagger(rdm2_ba, nmo, start, count); // vVvO
    LET(start, nocc_b, nocc_a, nocc_b, 0     );
    LET(count, nvir_b, nvir_a, nvir_b, nocc_a);
    rdm2_copyswap_e1e2(rdm2_ba, rdm2_ab, nmo, start, count); // VvOv

    LET(start, nocc_a, nocc_a, nocc_a, nocc_a);
    LET(count, nvir_a, nvir_a, nvir_a, nvir_a);
    copy_rdm2(rdm2_aa, nmo, start, count, 5, 5, 7, 7, "GABCD");
    LET(start, nocc_b, nocc_b, nocc_b, nocc_b);
    LET(count, nvir_b, nvir_b, nvir_b, nvir_b);
    copy_rdm2(rdm2_bb, nmo, start, count, 15, 15, 17, 17, "Gabcd");
    LET(start, nocc_a, nocc_b, nocc_a, nocc_b);
    LET(count, nvir_a, nvir_b, nvir_a, nvir_b);
    copy_rdm2(rdm2_ab, nmo, start, count, 28, 28, 28, 28, "GAbCd");
    rdm2_copyswap_e1e2(rdm2_ab, rdm2_ba, nmo, start, count);

    dpd_close(0);
    chkpt_close();
    exit_io();
    for (int i = 0; i < spaces.size(); i++) {
        free(spaces[i]);
    }
    free(cachefiles);
    free_int_matrix(cachelist);
}
