/*
 * modified from ccdensity.cc
 * remove unused calling
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libiwl/iwl.h>
#include <liboptions/liboptions.h>
#include <psi4-dec.h>
#include <cmath>
#include <psifiles.h>
#include "ccdensity/MOInfo.h"
#include "ccdensity/Params.h"
#include "ccdensity/Frozen.h"
#define EXTERN
#include "ccdensity/globals.h"
#include <libqt/qt.h>
namespace psi { namespace ccdensity {

void init_io(void);
void title(void);
void get_moinfo_light(void);
void get_frozen(void);
void get_params(Options& options);
void exit_io(void);
void onepdm(struct RHO_Params);
void sortone(struct RHO_Params);
void twopdm(void);
void energy(struct RHO_Params);
void resort_tei(void);
void resort_gamma(void);
void lag(struct RHO_Params rho_params);
void build_X(void);
void build_A(void);
void build_Z(void);
void relax_I(void);
void relax_D(struct RHO_Params rho_params);
void sortI(void);
void fold_RHF_light(struct RHO_Params rho_params);
void fold_UHF_light(struct RHO_Params rho_params);
void fold(struct RHO_Params rho_params);
void deanti(struct RHO_Params rho_params);
void add_ref_RHF(struct iwlbuf *);
void add_ref_ROHF(struct iwlbuf *);
void add_ref_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *);
void add_core_ROHF(struct iwlbuf *);
void add_core_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *);
void dump_RHF(struct iwlbuf *, struct RHO_Params rho_params);
void dump_ROHF(struct iwlbuf *, struct RHO_Params rho_params);
void dump_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *, struct RHO_Params rho_params);
void kinetic(void);
void dipole(void);
void probable(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void setup_LR(struct RHO_Params);
void G_build(void);
void x_oe_intermediates(struct RHO_Params);
void x_onepdm(struct RHO_Params);
void x_te_intermediates(void);
void x_Gijkl(void);
void x_Gabcd(void);
void x_Gibja(void);
void x_Gijka(void);
void x_Gijab(void);
void x_Gciab(void);
void V_build_x(void);
void x_xi1(void);
void x_xi_zero(void);
void x_xi2(void);
void x_xi_oe_intermediates(void);
void G_norm(void);
void zero_onepdm(struct RHO_Params rho_params);
void zero_twopdm(void);
void get_rho_params(Options& options);
void get_td_params(Options& options);
void td_setup(struct TD_Params S);
void tdensity(struct TD_Params S);
void td_print(void);
void oscillator_strength(struct TD_Params *S);
void rotational_strength(struct TD_Params *S);
void ael(struct RHO_Params *rho_params);
void cleanup_light(void);
void td_cleanup(void);
void x_oe_intermediates_rhf(struct RHO_Params rho_params);
void x_te_intermediates_rhf(void);
void x_xi_intermediates(void);
void V_build(void);
void ex_tdensity(char hand, struct TD_Params S, struct TD_Params U);
void ex_td_setup(struct TD_Params S, struct TD_Params U);
void ex_td_cleanup();
void ex_oscillator_strength(struct TD_Params *S, struct TD_Params *U, struct XTD_Params *xtd_data);
void ex_rotational_strength(struct TD_Params *S, struct TD_Params *U, struct XTD_Params *xtd_data);
void ex_td_print(std::vector<struct XTD_Params>);

PsiReturnType ccdensity_light(Options& options)
{
  int i;
  int **cachelist, *cachefiles;
  struct iwlbuf OutBuf;
  struct iwlbuf OutBuf_AA, OutBuf_BB, OutBuf_AB;
  dpdfile2 D;
  double tval;

  init_io();
  title();
  /*  get_frozen(); */
  get_params( options );
  get_moinfo_light();
  get_rho_params(options);

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    std::vector<int*> spaces;
    spaces.push_back(moinfo.occpi);
    spaces.push_back(moinfo.occ_sym);
    spaces.push_back(moinfo.virtpi);
    spaces.push_back(moinfo.vir_sym);
    delete dpd_list[0];
    dpd_list[0] = new DPD(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
             2, spaces);
    dpd_set_default(0);

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
  }

  for (i=0; i<params.nstates; ++i) {

    /* CC_GLG will contain L, or R0*L + Zeta, if relaxed and zeta is available */
    /* CC_GL will contain L */
    setup_LR(rho_params[i]);

    /* Calculate Xi, put Xi in EOM_XI, and quit */
    if ( params.calc_xi ) {
      /* these intermediates go into EOM_TMP and are used to compute Xi;
         they may be reused to compute the excited-state density matrix */
      if (params.ref == 0) {
        x_oe_intermediates_rhf(rho_params[i]);
        x_te_intermediates_rhf();
      }
      else {
        x_oe_intermediates(rho_params[i]);
        x_te_intermediates();
      }
      x_xi_intermediates(); /*Xi intermediates put in EOM_TMP_XI */
      x_xi_zero(); /* make blank Xi */
      x_xi1();
      x_xi2();
      dpd_close(0);
      if(params.ref == 2) cachedone_uhf(cachelist);
      else cachedone_rhf(cachelist);
      free(cachefiles);
      cleanup_light();
      psio_close(PSIF_EOM_TMP_XI,0); /* delete EOM_TMP_XI */
      psio_open(PSIF_EOM_TMP_XI,PSIO_OPEN_NEW);
      exit_io();
      return Success;
    }

    /* compute ground state parts of onepdm or put zeroes there */
    if ( ((rho_params[i].L_irr == rho_params[i].G_irr) || (params.use_zeta)) ) {
      zero_onepdm(rho_params[i]);
      onepdm(rho_params[i]);
    }
    else
      zero_onepdm(rho_params[i]);

    /* if the one-electron excited-state intermediates are not already on disk (from a Xi
       calculation, compute them.  They are nearly all necessary to compute the excited-state
       onepdm. Then complete excited-state onepdm.*/
    if (!rho_params[i].R_ground) {
      x_oe_intermediates(rho_params[i]); /* change to x_oe_intermediates_rhf() when rho gets spin-adapted */
      x_onepdm(rho_params[i]);
    }

    /* begin construction of twopdm */
    if (!params.onepdm) {

      /* Compute intermediates for construction of ground-state twopdm */
      if ( (params.L_irr == params.G_irr) || (params.use_zeta) ) {
        V_build(); /* uses CC_GLG, writes tau2*L2 to CC_MISC */
        G_build(); /* uses CC_GLG, writes t2*L2 to CC_GLG */
      }

      /* Compute ground-state twopdm or ground-state-like contributions to the excited twodpm */
      if ( (params.L_irr == params.G_irr) || (params.use_zeta) )
        twopdm();
      else
        zero_twopdm();

      /* Compute intermediates for construction of excited-state twopdm */
      if (!params.ground) {
        x_te_intermediates(); /* change to x_te_intermediates_rhf() when rho gets spin-adapted */
        V_build_x(); /* uses CC_GL, writes t2*L2 to EOM_TMP */

        /* add in non-R0 parts of onepdm and twopdm */
        x_Gijkl();
        x_Gabcd();
        x_Gibja();
        x_Gijka();
        x_Gciab();
        x_Gijab();
      }
    }

    //sortone(rho_params[i]); /* puts full 1-pdm into moinfo.opdm */ // check me if this neccessary

    if (!params.onepdm) {
      //if(!params.aobasis) energy(rho_params[i]); // regular CC-DM

//      kinetic(); /* puts kinetic energy integrals into MO basis */

      lag(rho_params[i]); /* builds the orbital lagrangian pieces, I */
//
//      /* dpd_init(1, moinfo.nirreps, params.memory, 2, frozen.occpi, frozen.occ_sym,
//         frozen.virtpi, frozen.vir_sym); */
//
//      /*  if(moinfo.nfzc || moinfo.nfzv) {
//          resort_gamma();
//          resort_tei();
//          } */
//
      build_X(); /* builds orbital rotation gradient X */
      build_A(); /* construct MO Hessian A */
      build_Z(); /* solves the orbital Z-vector equations */

      relax_I(); /* adds orbital response contributions to Lagrangian */

//      if (params.relax_opdm) {
//        relax_D(rho_params[i]); /* adds orbital response contributions to onepdm */
//      }
      sortone(rho_params[i]); /* builds large moinfo.opdm matrix */
      sortI(); /* builds large lagrangian matrix I */

      //fock-adjust DM
      if (params.ref == 0) {
          fold_RHF_light(rho_params[i]);
      } else if (params.ref == 2) {
          fold_UHF_light(rho_params[i]);
      } else {
          fold(rho_params[i]);
      }
#if defined ccdensity_deanti
      deanti(rho_params[i]);
#endif
    }

    /*  dpd_close(0); dpd_close(1); */

    if(params.ref == 0) { /** RHF **/

      iwl_buf_init(&OutBuf, PSIF_MO_TPDM, params.tolerance, 0, 0);

      add_core_ROHF(&OutBuf);
      add_ref_RHF(&OutBuf);

      // ==> One-Electron Properties <== //
//      fprintf(outfile, "  ==> Properties: Root %d <==\n\n", i);
//      dipole();

      if(params.onepdm_grid_dump) dx_write(options, moinfo.opdm);
 
      dump_RHF(&OutBuf, rho_params[i]);

      iwl_buf_flush(&OutBuf, 1);
      iwl_buf_close(&OutBuf, 1);

    }
    else if(params.ref == 1) { /** ROHF **/

      iwl_buf_init(&OutBuf, PSIF_MO_TPDM, params.tolerance, 0, 0);

      add_core_ROHF(&OutBuf);
      add_ref_ROHF(&OutBuf);
      dump_ROHF(&OutBuf, rho_params[i]);

      iwl_buf_flush(&OutBuf, 1);
      iwl_buf_close(&OutBuf, 1);
    }
    else if(params.ref == 2) { /** UHF **/

      iwl_buf_init(&OutBuf_AA, PSIF_MO_AA_TPDM, params.tolerance, 0, 0);
      iwl_buf_init(&OutBuf_BB, PSIF_MO_BB_TPDM, params.tolerance, 0, 0);
      iwl_buf_init(&OutBuf_AB, PSIF_MO_AB_TPDM, params.tolerance, 0, 0);

      /*    add_core_UHF(&OutBuf_AA, &OutBuf_BB, &OutBuf_AB); */
      add_ref_UHF(&OutBuf_AA, &OutBuf_BB, &OutBuf_AB);
      dump_UHF(&OutBuf_AA, &OutBuf_BB, &OutBuf_AB, rho_params[i]);

      iwl_buf_flush(&OutBuf_AA, 1);
      iwl_buf_flush(&OutBuf_BB, 1);
      iwl_buf_flush(&OutBuf_AB, 1);
      iwl_buf_close(&OutBuf_AA, 1);
      iwl_buf_close(&OutBuf_BB, 1);
      iwl_buf_close(&OutBuf_AB, 1);
    }

    psio_close(PSIF_CC_TMP,0);   psio_open(PSIF_CC_TMP,PSIO_OPEN_NEW);
    psio_close(PSIF_EOM_TMP0,0); psio_open(PSIF_EOM_TMP0,PSIO_OPEN_NEW);
    psio_close(PSIF_EOM_TMP1,0); psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);
    psio_close(PSIF_CC_GLG,0);   psio_open(PSIF_CC_GLG,PSIO_OPEN_NEW);
    psio_close(PSIF_CC_GL,0);    psio_open(PSIF_CC_GL,PSIO_OPEN_NEW);
    psio_close(PSIF_CC_GR,0);    psio_open(PSIF_CC_GR,PSIO_OPEN_NEW);
    if (!params.calc_xi) {
      psio_close(PSIF_EOM_TMP,0);
      psio_open(PSIF_EOM_TMP,PSIO_OPEN_NEW);
    }
  }

  free(rho_params);
  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup_light();
  exit_io();
  return Success;
}
}}
