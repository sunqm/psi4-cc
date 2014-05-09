
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include <libmints/wavefunction.h>
#include "cclambda/MOInfo.h"
#include "cclambda/Params.h"
#include "cclambda/Local.h"
#define EXTERN
#include "cclambda/globals.h"
#include "cclambda/cclambda.h"

namespace psi { namespace cclambda {

void init_io(void);
void title(void);
void get_moinfo_light(void);
void get_params(Options& options);
void cleanup_light(void);
void init_amps(struct L_Params L_params);
double pseudoenergy(struct L_Params L_params);
void exit_io(void);
void G_build(int L_irr);
void L1_build(struct L_Params L_params);
void L2_build(struct L_Params L_params);
void sort_amps(int L_irr);
void Lsave(int L_irr);
void Lnorm(struct L_Params L_params);
void Lmag(void);
void update(void);
int converged(int L_irr);
void diis(int iter, int L_irr);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void denom(struct L_Params);
void overlap(int L_irr);
void overlap_LAMPS(struct L_Params L_params);
void Lsave_index(struct L_Params L_params);
void Lamp_write(struct L_Params L_params);
void check_ortho(struct L_Params *pL_params);
void projections(struct L_Params *pL_params);
void L_zero(int irrep);
void c_clean(dpdfile2 *LIA, dpdfile2 *Lia, dpdbuf4 *LIJAB, dpdbuf4 *Lijab, dpdbuf4 *LIjAb);
void L_clean(struct L_Params pL_params);
void zeta_norm(struct L_Params pL_params);
void spinad_amps(void);
void status(const char *, FILE *);
void hbar_extra(void);
void ortho_Rs(struct L_Params *pL_params, int current_L);

void cc2_L1_build(struct L_Params L_params);
void cc2_L2_build(struct L_Params L_params);
void cc2_Gai_build(int L_irr);
void cc2_hbar_extra(void);

void cc3_t3z(void);
void cc3_t3x(void);
void cc3_l3l2(void);
void cc3_l3l1(void);

void local_init(void);
void local_done(void);
}} //namespace psi::cclambda

namespace psi { namespace cclambda {

PsiReturnType cclambda_light(Options& options)
{
  int done=0, i, root_L_irr;
  int **cachelist, *cachefiles;
  dpdfile2 L1;

  init_io();
  title();
  moinfo.iter=0;
  get_moinfo_light();
  get_params(options);

  /* throw any existing CC_LAMBDA, CC_DENOM away */
  /* Do this only if we're not running an analytic gradient on the
     ground state. Keeping the files around should allow us to
     restart from old Lambda amplitudes. -TDC, 11/2007 */
  if(!(params.dertype==1 && !cc_excited(params.wfn))) {
    fprintf(outfile, "\tDeleting old CC_LAMBDA data.\n");
    psio_close(PSIF_CC_LAMBDA,0);
    psio_open(PSIF_CC_LAMBDA,PSIO_OPEN_NEW);
    psio_close(PSIF_CC_DENOM,0);
    psio_open(PSIF_CC_DENOM,PSIO_OPEN_NEW);
  }

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    std::vector<int*> spaces;
    spaces.push_back(moinfo.occpi);
    spaces.push_back(moinfo.occ_sym);
    spaces.push_back(moinfo.virtpi);
    spaces.push_back(moinfo.vir_sym);
    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 2, spaces);

    if(params.aobasis) { /* Set up new DPD for AO-basis algorithm */
        std::vector<int*> aospaces;
        aospaces.push_back(moinfo.occpi);
        aospaces.push_back(moinfo.occ_sym);
        aospaces.push_back(moinfo.sopi);
        aospaces.push_back(moinfo.sosym);
        dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 2, aospaces);
        dpd_set_default(0);
    }

  }
  else if(params.ref == 2) { /** UHF **/

    cachelist = cacheprep_uhf(params.cachelev, cachefiles);
    std::vector<int*> spaces;
    spaces.push_back(moinfo.aoccpi);
    spaces.push_back(moinfo.aocc_sym);
    spaces.push_back(moinfo.avirtpi);
    spaces.push_back(moinfo.avir_sym);
    spaces.push_back(moinfo.boccpi);
    spaces.push_back(moinfo.bocc_sym);
    spaces.push_back(moinfo.bvirtpi);
    spaces.push_back(moinfo.bvir_sym);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 4, spaces);

    if(params.aobasis) { /* Set up new DPD's for AO-basis algorithm */
        std::vector<int*> aospaces;
        aospaces.push_back(moinfo.aoccpi);
        aospaces.push_back(moinfo.aocc_sym);
        aospaces.push_back(moinfo.sopi);
        aospaces.push_back(moinfo.sosym);
        aospaces.push_back(moinfo.boccpi);
        aospaces.push_back(moinfo.bocc_sym);
        aospaces.push_back(moinfo.sopi);
        aospaces.push_back(moinfo.sosym);
        dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 4, aospaces);
        dpd_set_default(0);
    }
  }

  if(params.local) local_init();

  if(params.ref == 0) {
    if (params.wfn == "CC2" || params.wfn == "EOM_CC2")
      cc2_hbar_extra();
    else
      hbar_extra();
  }

  /* CC3: Z-build */
  if(params.wfn == "CC3") cc3_t3z();

  for (i=0; i<params.nstates; ++i) {

    /* delete and reopen intermediate files */
    psio_close(PSIF_CC_TMP,0); psio_close(PSIF_CC_TMP0,0);
    psio_close(PSIF_CC_TMP1,0); psio_close(PSIF_CC_TMP2,0);
    psio_open(PSIF_CC_TMP,0); psio_open(PSIF_CC_TMP0,0);
    psio_open(PSIF_CC_TMP1,0); psio_open(PSIF_CC_TMP2,0);
    /* Keep the old lambda amps if this is a ground-state geomopt */
    if(!(params.dertype==1 && !cc_excited(params.wfn))) {
      psio_close(PSIF_CC_LAMBDA,0);
      psio_open(PSIF_CC_LAMBDA,PSIO_OPEN_NEW);
      psio_close(PSIF_CC_DENOM,0); /* aren't these recomputed anyway - perhaps should always delete? */
      psio_open(PSIF_CC_DENOM,PSIO_OPEN_NEW);
    }

    fprintf(outfile,"\tSymmetry of left-hand state: %s\n",
            moinfo.labels[ moinfo.sym^(pL_params[i].irrep) ]);
    fprintf(outfile,"\tSymmetry of left-hand eigenvector: %s\n",
            moinfo.labels[(pL_params[i].irrep)]);

    denom(pL_params[i]); /* uses L_params.cceom_energy for excited states */
    init_amps(pL_params[i]); /* uses denominators for initial zeta guess */

    fprintf(outfile, "\n\t          Solving Lambda Equations\n");
    fprintf(outfile, "\t          ------------------------\n");
    fprintf(outfile, "\tIter     PseudoEnergy or Norm         RMS  \n");
    fprintf(outfile, "\t----     ---------------------     --------\n");

    moinfo.lcc = pseudoenergy(pL_params[i]);
    update();

    for(moinfo.iter=1 ; moinfo.iter <= params.maxiter; moinfo.iter++) {
      sort_amps(pL_params[i].irrep);

      /* must zero New L before adding RHS */
      L_zero(pL_params[i].irrep);

      if(params.wfn == "CC3") cc3_t3x();

      if(params.wfn == "CC2" || params.wfn == "EOM_CC2") {

    cc2_Gai_build(pL_params[i].irrep);
    cc2_L1_build(pL_params[i]);
    if(params.print & 2) status("L1 amplitudes", outfile);
    cc2_L2_build(pL_params[i]);

      }
      else {
    G_build(pL_params[i].irrep);
    L1_build(pL_params[i]);
    if(params.print & 2) status("L1 amplitudes", outfile);
    L2_build(pL_params[i]);

    if(params.wfn == "CC3") {
      cc3_l3l2();
      cc3_l3l1();
    }
      }

      if (params.ref == 1) L_clean(pL_params[i]);
      if (params.nstates > 2) ortho_Rs(pL_params, i);

      if(converged(pL_params[i].irrep)) {
        done = 1;  /* Boolean for convergence */
        Lsave(pL_params[i].irrep); /* copy "New L" to "L" */
        moinfo.lcc = pseudoenergy(pL_params[i]);
        update();
        if (!pL_params[i].ground && !params.zeta) {
          Lnorm(pL_params[i]); /* normalize against R */
        }
        Lsave_index(pL_params[i]); /* save Ls with indices in LAMPS */
        Lamp_write(pL_params[i]); /* write out largest  Ls */
    /* sort_amps(); to be done by later functions */
        fprintf(outfile, "\n\tIterations converged.\n");
        fflush(outfile);
        moinfo.iter = 0;
        break;
      }

      if(params.diis) diis(moinfo.iter, pL_params[i].irrep);
      Lsave(pL_params[i].irrep);
      moinfo.lcc = pseudoenergy(pL_params[i]);
      update();
    }
    fprintf(outfile, "\n");
    if(!done) {
      fprintf(outfile, "\t ** Lambda not converged to %2.1e ** \n",
          params.convergence);
      fflush(outfile);
      dpd_close(0);
      cleanup_light();
      exit_io();
      throw PsiException("cclambda: error", __FILE__, __LINE__);
    }
    if (pL_params[i].ground)
      overlap(pL_params[i].irrep);
  }

  if (params.zeta) {
    zeta_norm(pL_params[0]);
  }
  else if (params.nstates > 1) { /* some excited states are present */
    check_ortho(pL_params);
    projections(pL_params);
  }

  if(params.local) local_done();

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup_light();
  exit_io();
  return Success;
}

}} // namespace psi::cclambda
