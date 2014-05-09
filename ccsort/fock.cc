#include <cstdio>
#include <cstdlib>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include "ccsort/Params.h"
#include "ccsort/MOInfo.h"
#define EXTERN
#include "ccsort/globals.h"

#include <libmints/wavefunction.h>
#include <libtrans/mospace.h>
#include <libmints/matrix.h>

#define LOWERTRI_INDEX(I,J)     ((I) > (J) ? ((I)*((I)+1)/2+(J)) : ((J)*((J)+1)/2+(I)))

namespace psi { namespace ccsort {

void fock_rhf(void);
void fock_uhf(void);

void fock(void)
{
  if(params.ref == 2) fock_uhf();
  else fock_rhf();
}

void fock_uhf(void)
{
  int h, nirreps;
  int i, j, i1, j1;
  int a, b, a1, b1;
  int *aoccpi, *boccpi;
  int *avirtpi, *bvirtpi;
  int *frdocc;
  dpdfile2 fIJ, fij, fAB, fab, fIA, fia;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi;
  boccpi = moinfo.boccpi;
  avirtpi = moinfo.avirtpi;
  bvirtpi = moinfo.bvirtpi;
  frdocc = moinfo.frdocc;

  // Take Fock matrix from the wavefunction and transform it to the MO basis
  chkpt_init(PSIO_OPEN_OLD);
  double *ea = chkpt_rd_alpha_evals();
  double *eb = chkpt_rd_beta_evals();
  chkpt_close();

  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
  global_dpd_->file2_mat_init(&fIJ);
  global_dpd_->file2_mat_init(&fij);

  for(h=0; h < nirreps; h++) {
      for(i=0; i < aoccpi[h]; i++) {
          fIJ.matrix[h][i][i] = ea[i+frdocc[h]];
      }
      for(i=0; i < boccpi[h]; i++) {
          fij.matrix[h][i][i] = eb[i+frdocc[h]];
      }
  }
  global_dpd_->file2_mat_wrt(&fIJ);
  global_dpd_->file2_mat_wrt(&fij);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fij);

  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");

  global_dpd_->file2_mat_init(&fAB);
  global_dpd_->file2_mat_init(&fab);

  for(h=0; h < nirreps; h++) {
      for(a=0; a < avirtpi[h]; a++) {
          fAB.matrix[h][a][a] = ea[a+frdocc[h]+aoccpi[h]];
      }
      for(a=0; a < bvirtpi[h]; a++) {
          fab.matrix[h][a][a] = eb[a+frdocc[h]+boccpi[h]];
      }
  }

  global_dpd_->file2_mat_wrt(&fAB);
  global_dpd_->file2_mat_wrt(&fab);
  global_dpd_->file2_close(&fAB);
  global_dpd_->file2_close(&fab);

  /* Prepare the alpha and beta occ-vir Fock matrix files */
  global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
  global_dpd_->file2_mat_init(&fIA);
  global_dpd_->file2_mat_init(&fia);

  /* Close the alpha and beta occ-vir Fock matrix files */
  global_dpd_->file2_mat_wrt(&fIA);
  global_dpd_->file2_mat_wrt(&fia);
  global_dpd_->file2_mat_close(&fIA);
  global_dpd_->file2_mat_close(&fia);
  global_dpd_->file2_close(&fIA);
  global_dpd_->file2_close(&fia);

  free(ea);
  free(eb);
}

void fock_rhf(void)
{
  int h;
  int i,j,a,b, a1, b1;
  int nirreps;
  int *occpi, *virtpi;
  int *openpi;
  int *frdocc;
  int *fruocc;
  dpdfile2 fIJ, fij, fAB, fab, fIA, fia;

  // Take Fock matrix from the wavefunction and transform it to the MO basis

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;
  frdocc = moinfo.frdocc;
  fruocc = moinfo.fruocc;

  chkpt_init(PSIO_OPEN_OLD);
  double *evals = chkpt_rd_evals();
  chkpt_close();

  /* Prepare the alpha and beta occ-occ Fock matrix files */
  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
  global_dpd_->file2_mat_init(&fIJ);
  global_dpd_->file2_mat_init(&fij);

  /* One-electron (frozen-core) contributions */
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
          fIJ.matrix[h][i][i] = evals[i+frdocc[h]];
      }
      for(i=0; i < (occpi[h]-openpi[h]); i++) {
          fij.matrix[h][i][i] = evals[i+frdocc[h]];
      }
  }

  /* Close the alpha and beta occ-occ Fock matrix files */
  global_dpd_->file2_mat_wrt(&fIJ);
  global_dpd_->file2_mat_wrt(&fij);
  global_dpd_->file2_mat_close(&fIJ);
  global_dpd_->file2_mat_close(&fij);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fij);

  /* Prepare the alpha and beta vir-vir Fock matrix files */
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
  global_dpd_->file2_mat_init(&fAB);
  global_dpd_->file2_mat_init(&fab);

  /* One-electron (frozen-core) contributions */
  for(h=0; h < nirreps; h++) {

      for(a=0; a < (virtpi[h] - openpi[h]); a++) {
          fAB.matrix[h][a][a] = evals[a+frdocc[h]+occpi[h]];
      }

      for(a=0; a < (virtpi[h] - openpi[h]); a++) {
          fab.matrix[h][a][a] = evals[a+frdocc[h]+occpi[h]];
      }

      for(a=0; a < openpi[h]; a++) {
          a1 = virtpi[h] - openpi[h] + a;
          b1 = frdocc[h] + occpi[h] - openpi[h] + a;
          fab.matrix[h][a1][a1] = evals[b1];
      }

//      for(a=0; a < (virtpi[h] - openpi[h]); a++)
//          for(b=0; b < openpi[h]; b++) {
//              a1 = a + frdocc[h] + occpi[h];
//              b1 = b + frdocc[h] + occpi[h] - openpi[h];
//              fab.matrix[h][a][virtpi[h] - openpi[h] + b] =
//                  trif[LOWERTRI_INDEX(a1,b1)];
//          }
//
//      for(a=0; a < openpi[h]; a++)
//          for(b=0; b < (virtpi[h] - openpi[h]); b++) {
//              a1 = a + frdocc[h] + occpi[h] - openpi[h];
//              b1 = b + frdocc[h] + occpi[h];
//              fab.matrix[h][virtpi[h] - openpi[h] + a][b] =
//                  trif[LOWERTRI_INDEX(a1,b1)];
//          }
  }

  /* Close the alpha and beta vir-vir Fock matrix files */
  global_dpd_->file2_mat_wrt(&fAB);
  global_dpd_->file2_mat_wrt(&fab);
  global_dpd_->file2_mat_close(&fAB);
  global_dpd_->file2_mat_close(&fab);
  global_dpd_->file2_close(&fAB);
  global_dpd_->file2_close(&fab);

  /* Prepare the alpha and beta occ-vir Fock matrix files */
  global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
  global_dpd_->file2_mat_init(&fIA);
  global_dpd_->file2_mat_init(&fia);

//  /* One-electron (frozen-core) contributions */
//  for(h=0; h < nirreps; h++) {
//
//      for(i=0; i < occpi[h]; i++)
//          for(a=0; a < (virtpi[h] - openpi[h]); a++) {
//              i1 = i + frdocc[h];
//              a1 = a + frdocc[h] + occpi[h];
//              fIA.matrix[h][i][a] = trif[LOWERTRI_INDEX(i1,a1)];
//          }
//
//      for(i=0; i < (occpi[h] - openpi[h]); i++)
//          for(a=0; a < (virtpi[h] - openpi[h]); a++) {
//              i1 = i + frdocc[h];
//              a1 = a + frdocc[h] + occpi[h];
//              fia.matrix[h][i][a] = trif[LOWERTRI_INDEX(i1,a1)];
//          }
//
//      for(i=0; i < (occpi[h] - openpi[h]); i++)
//          for(a=0; a < openpi[h]; a++) {
//              i1 = i + frdocc[h];
//              a1 = a + frdocc[h] + occpi[h] - openpi[h];
//              fia.matrix[h][i][virtpi[h] - openpi[h] + a] =
//                  trif[LOWERTRI_INDEX(i1,a1)];
//          }
//  }

  /* Close the alpha and beta occ-vir Fock matrix files */
  global_dpd_->file2_mat_wrt(&fIA);
  global_dpd_->file2_mat_wrt(&fia);
  global_dpd_->file2_mat_close(&fIA);
  global_dpd_->file2_mat_close(&fia);
  global_dpd_->file2_close(&fIA);
  global_dpd_->file2_close(&fia);

  free(evals);
}

}} // namespace psi::ccsort
