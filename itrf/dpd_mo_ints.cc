/*
 * dpd_mo_ints.cc
 * use psi4 DPD format to save the 1e/2e integrals in MO representation
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>

#include <psi4-dec.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <liboptions/liboptions.h>
#include <psifiles.h>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include "psi4itrf.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define INDEX(i,j) (i > j ? i * (i + 1) / 2 + j : j + (j + 1) / 2 + i)
#define LOWERTRI_INDEX(I,J)     ((I) > (J) ? ((I)*((I)+1)/2+(J)) : ((J)*((J)+1)/2+(I)))

#define NUM_SUBSPACES_RHF 2
#define NUM_SUBSPACES_UHF 4

using namespace psi;
namespace psi {
    namespace ccenergy {
        int **cacheprep_rhf(int level, int *cachefiles);
        int **cacheprep_uhf(int level, int *cachefiles);
    }
}

void assign_sym(int nirreps, int *orbpi, int *orb_sym)
{
    int h, p, q;
    for(h = 0, q = 0; h < nirreps; h++)
        for(p = 0; p < orbpi[h]; p++) {
            orb_sym[q++] = h;
        }
}

// see transqt2/transqt.cc
void save_mo_ints_1e(const char *label, const double *h1e, int nmo)
{
    int ntri_mo = nmo * (nmo + 1) / 2;
    double *buf = new double[ntri_mo];
    int i, j, n;
    n = 0;
    for (i = 0; i < nmo; i++) {
        for (j = 0; j <= i; j++) {
            buf[n] = h1e[i*nmo+j];
            n++;
        }
    }
    iwl_wrtone(PSIF_OEI, label, ntri_mo, buf);
    delete[] buf;
}

void save_mo_ints_2e(int filenum, const double *eri, int nmo)
{
    int n_pairs = nmo * (nmo + 1) / 2;

    struct iwlbuf MBuff;
    double tolerance = 1e-14; // = Process::environment.options.get_double("INTS_TOLERANCE");
    iwl_buf_init(&MBuff, filenum, tolerance, 0, 0);

    int p, q, r, s, pq, rs;

    for(p = 0; p < nmo; p++) {
        for(q = 0; q <= p; q++) {
            pq = LOWERTRI_INDEX(p, q);
            for(r = 0; r < nmo; r++) {
                for(s = 0; s <= r; s++) {
                    rs = LOWERTRI_INDEX(r, s);
                    iwl_buf_wrt_val(&MBuff, p, q, r, s, eri[pq*n_pairs+rs],
                                    0, "outfile", 0);
                }
            }
        }
    }

    iwl_buf_flush(&MBuff, 1);
    iwl_buf_close(&MBuff, 1);
}

void init_dpd_rhf(int nmo, std::vector<int *> spaces,
                  int *cachefiles, int **cachelist)
{
    int nirreps = chkpt_rd_nirreps();
    int *mopi = chkpt_rd_orbspi();
    int h, p, q, count, i;

    //int *frdocc = chkpt_rd_frzcpi();
    //int *fruocc = chkpt_rd_frzvpi();
    //int *actpi = init_int_array(nirreps);
    //int nactive = 0;
    //for(h=0; h < nirreps; h++) {
    //    actpi[h] = mopi[h] - frdocc[h] - fruocc[h];
    //    nactive += actpi[h];
    //}
    //int *actsym = init_int_array(nactive);
    //for(h=0,count=0; h < nirreps; h++)
    //    for(i=0; i < actpi[h]; i++,count++)
    //        actsym[count] = h;

    int *occpi = chkpt_rd_clsdpi();
    int *virpi = (int *)malloc(sizeof(int) * nirreps);
    for (h = 0; h < nirreps; h++) {
        virpi[h] = mopi[h] - occpi[h];
    }
    int *occ_sym = (int *)malloc(sizeof(int)*nmo*nmo);
    int *vir_sym = (int *)malloc(sizeof(int)*nmo*nmo);
    assign_sym(nirreps, occpi, occ_sym);
    assign_sym(nirreps, virpi, vir_sym);
    spaces.push_back(occpi);
    spaces.push_back(occ_sym);
    spaces.push_back(virpi);
    spaces.push_back(vir_sym);

    dpd_file4_cache_entry *priority = NULL;
    dpd_init(0, nirreps, Process::environment.get_memory(),
             0, cachefiles, cachelist, priority, NUM_SUBSPACES_RHF,
             spaces);

    free(mopi);
}

void init_dpd_uhf(int nmo, std::vector<int *> spaces,
                  int *cachefiles, int **cachelist)
{
    int nirreps = chkpt_rd_nirreps();
    int *mopi = chkpt_rd_orbspi();
    int h, p, q, count, i;
    int *clsdpi = chkpt_rd_clsdpi();
    int *openpi = chkpt_rd_openpi();
    int *aoccpi = (int *)malloc(sizeof(int) * nirreps);
    int *avirpi = (int *)malloc(sizeof(int) * nirreps);
    int *boccpi = (int *)malloc(sizeof(int) * nirreps);
    int *bvirpi = (int *)malloc(sizeof(int) * nirreps);
    for (h = 0; h < nirreps; h++) {
        aoccpi[h] = clsdpi[h] + openpi[h];
        avirpi[h] = mopi[h] - aoccpi[h];
        boccpi[h] = clsdpi[h];
        bvirpi[h] = mopi[h] - boccpi[h];
    }
    int *aocc_sym = (int *)malloc(sizeof(int)*nmo*nmo);
    int *avir_sym = (int *)malloc(sizeof(int)*nmo*nmo);
    int *bocc_sym = (int *)malloc(sizeof(int)*nmo*nmo);
    int *bvir_sym = (int *)malloc(sizeof(int)*nmo*nmo);
    assign_sym(nirreps, aoccpi, aocc_sym);
    assign_sym(nirreps, avirpi, avir_sym);
    assign_sym(nirreps, boccpi, bocc_sym);
    assign_sym(nirreps, bvirpi, bvir_sym);
    // see libdpd/init.cc. spaces[2*i] saves the number of orbitals per irrep
    // space[2*i+1] saves the irrep for each orbital
    spaces.push_back(aoccpi); spaces.push_back(aocc_sym);
    spaces.push_back(avirpi); spaces.push_back(avir_sym);
    spaces.push_back(boccpi); spaces.push_back(bocc_sym);
    spaces.push_back(bvirpi); spaces.push_back(bvir_sym);

    dpd_file4_cache_entry *priority = NULL;
    dpd_init(0, nirreps, Process::environment.get_memory(),
             0, cachefiles, cachelist, priority, NUM_SUBSPACES_UHF,
             spaces);

    free(mopi);
}


// simulate transqt2/transqt.cc
/* h1e is hcore
 * eri carries 4-fold symmetry ijkl: i>j,k>l */
void dpd_mo_ints_rhf(const double *h1e, const double *eri, int nmo)
{
    psio_open(PSIF_CC_INFO, PSIO_OPEN_NEW);
    chkpt_init(PSIO_OPEN_OLD);
    std::vector<int*> spaces;
    int *cachefiles = init_int_array(PSIO_MAXUNIT);
    int **cachelist = ccenergy::cacheprep_rhf(CACHELEVEL, cachefiles);
    init_dpd_rhf(nmo, spaces, cachefiles, cachelist);

    save_mo_ints_1e(PSIF_MO_FZC, h1e, nmo);
    save_mo_ints_2e(PSIF_MO_TEI, eri, nmo);

    dpd_close(0);
    chkpt_close();
    psio_close(PSIF_CC_INFO,1);
    for (int i = 0; i < spaces.size(); i++) {
        free(spaces[i]);
    }
    free(cachefiles);
    free_int_matrix(cachelist);
}

void dpd_mo_ints_uhf(double *hcore_a, double *hcore_b,
                     double *eri_aa, double *eri_bb, double *eri_ab,
                     int nmo)
{
    psio_open(PSIF_CC_INFO, PSIO_OPEN_NEW);
    chkpt_init(PSIO_OPEN_OLD);
    std::vector<int*> spaces;
    int *cachefiles = init_int_array(PSIO_MAXUNIT);
    int **cachelist = ccenergy::cacheprep_uhf(CACHELEVEL, cachefiles);
    init_dpd_uhf(nmo, spaces, cachefiles, cachelist);

    save_mo_ints_1e(PSIF_MO_A_FZC, hcore_a, nmo);
    save_mo_ints_1e(PSIF_MO_B_FZC, hcore_b, nmo);
    save_mo_ints_2e(PSIF_MO_AA_TEI, eri_aa, nmo);
    save_mo_ints_2e(PSIF_MO_BB_TEI, eri_bb, nmo);
    save_mo_ints_2e(PSIF_MO_AB_TEI, eri_ab, nmo);

    dpd_close(0);
    chkpt_close();
    psio_close(PSIF_CC_INFO,1);
    for (int i = 0; i < spaces.size(); i++) {
        free(spaces[i]);
    }
    free(cachefiles);
    free_int_matrix(cachelist);
}
