#include <stdio.h>
#include <libdpd/dpd.h>
#include "ccdensity/MOInfo.h"
#include "ccdensity/Params.h"
#include "ccdensity/Frozen.h"
#define EXTERN
#include "ccdensity/globals.h"

namespace psi { namespace ccdensity {

    /* FOLD_RHF(): Fold the RHF Fock matrix contributions to the energy
    ** (or energy derivative) into the two-particle density matrix.  Here
    ** we are trying to convert from an energy expression of the form:
    **
    ** E = sum_pq Dpq fpq + 1/4 sum_pqrs Gpqrs <pq||rs>
    **
    ** to the form:
    **
    ** E = sum_pq Dpq hpq + 1/4 sum_pqrs Gpqrs <pq||rs>
    **
    ** We do this by shifting some one-particle density matrix components
    ** into appropriate two-particle density matrix components:
    **
    ** G'pmrm = Dpr + 4 * Gpmrm
    **
    ** One problem is that we need to make sure the resulting density,
    ** G'pqrs, is still antisymmetric to permutation of p and q or r and
    ** s.  So, for example, for the Gimkm component we compute:
    **
    ** G'pmrm = Dpr + Gpmrm
    ** G'mprm = Dpr - Gmprm
    ** G'pmmr = Dpr - Gpmmr
    ** G'mpmr = Dpr + Gmpmr
    ** */

    void fold_RHF_light(struct RHO_Params rho_params)
    {
      int h, nirreps;
      int i, j, k, l, m, a, b;
      int I, J, K, L, M, A, B;
      int IM, JM, MI, MJ, MK, ML, MA, MB;
      int Gi, Gj, Gk, Gl, Gm, Ga, Gb;
      int *occpi, *virtpi;
      int *occ_off, *vir_off;
      int *occ_sym, *vir_sym;
      int *openpi;
      dpdfile2 D, D1, D2, F;
      dpdbuf4 G, Aints, E, C, DInts, FInts, BInts, G1, G2;
      double one_energy=0.0, two_energy=0.0, total_two_energy=0.0;
      double test_energy = 0.0, tmp;
      double this_energy;

      nirreps = moinfo.nirreps;
      occpi = moinfo.occpi; virtpi = moinfo.virtpi;
      occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
      occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
      openpi = moinfo.openpi;

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
      global_dpd_->file2_mat_init(&D);
      global_dpd_->file2_mat_rd(&D);

      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&G, h);
	global_dpd_->buf4_mat_irrep_rd(&G, h);

	for(Gm=0; Gm < nirreps; Gm++) {
	  Gi = Gj = h^Gm;

	  for(i=0; i < occpi[Gi]; i++) {
	    I = occ_off[Gi] + i;
	    for(j=0; j < occpi[Gj]; j++) {
	      J = occ_off[Gj] + j;
	      for(m=0; m < occpi[Gm]; m++) {
		M = occ_off[Gm] + m;

		IM = G.params->rowidx[I][M];
		JM = G.params->colidx[J][M];
		MI = G.params->rowidx[M][I];
		MJ = G.params->colidx[M][J];

		G.matrix[h][IM][JM] += D.matrix[Gi][i][j];
		G.matrix[h][MI][MJ] += D.matrix[Gi][i][j];
	      }
	    }
	  }
	}

	global_dpd_->buf4_mat_irrep_wrt(&G, h); 
	global_dpd_->buf4_mat_irrep_close(&G, h);
      }

      /* Generate spin-adapted Gijkl just for the energy calculation */
      global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijkl - Gijlk", 2);
      global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 0, 0, "2 Gijkl - Gijlk", -1);
      global_dpd_->buf4_close(&G);
      global_dpd_->file2_mat_close(&D);
      global_dpd_->file2_close(&D);

      global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
      global_dpd_->file2_mat_init(&D1);
      global_dpd_->file2_mat_rd(&D1);
      global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
      global_dpd_->file2_mat_init(&D2);
      global_dpd_->file2_mat_rd(&D2);

      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&G, h);
	global_dpd_->buf4_mat_irrep_rd(&G, h);

	for(Gm=0; Gm < nirreps; Gm++) {
	  Gi = Ga = h^Gm;

	  for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
	    I = occ_off[Gi] + i;
	    for(a=0; a < virtpi[Ga]; a++) {
	      A = vir_off[Ga] + a;
	      for(m=0; m < occpi[Gm]; m++) {
		M = occ_off[Gm] + m;

		MI = G.params->rowidx[M][I];
		MA = G.params->colidx[M][A];

		G.matrix[h][MI][MA] += 0.5 * (D1.matrix[Gi][i][a] +
					      D2.matrix[Gi][i][a]);
	      }
	    }
	  }
	}

	global_dpd_->buf4_mat_irrep_wrt(&G, h);
	global_dpd_->buf4_mat_irrep_close(&G, h);
      }

      /* Generate spin-adapted Gijka just for the energy calculation */
      global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijka - Gjika", 2);
      global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, qprs, 0, 10, "2 Gijka - Gjika", -1);
      global_dpd_->buf4_close(&G);

      global_dpd_->file2_mat_close(&D1);
      global_dpd_->file2_close(&D1);
      global_dpd_->file2_mat_close(&D2);
      global_dpd_->file2_close(&D2);

      /* Generate spin-adapted Gijab jut for energy calculation */
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
      global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gijab - Gijba", 2);
      global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 0, 5, "2 Gijab - Gijba", -1);
      global_dpd_->buf4_close(&G);
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "2 Gijab - Gijba");
      global_dpd_->buf4_close(&G);


      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
      global_dpd_->file2_mat_init(&D);
      global_dpd_->file2_mat_rd(&D);

      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&G, h);
	global_dpd_->buf4_mat_irrep_rd(&G, h);

	for(Gm=0; Gm < nirreps; Gm++) {
	  Ga = Gb = h^Gm;

	  for(b=0; b < (virtpi[Gb] - openpi[Gb]); b++) {
	    B = vir_off[Gb] + b;
	    for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
	      A = vir_off[Ga] + a;
	      for(m=0; m < occpi[Gm]; m++) {
		M = occ_off[Gm] + m;

		MB = G.params->rowidx[M][B];
		MA = G.params->colidx[M][A];

		G.matrix[h][MB][MA] += D.matrix[Ga][a][b];
	      }
	    }
	  }
	}

	global_dpd_->buf4_mat_irrep_wrt(&G, h);
	global_dpd_->buf4_mat_irrep_close(&G, h);
      }

      global_dpd_->buf4_close(&G);
      global_dpd_->file2_mat_close(&D);
      global_dpd_->file2_close(&D);

      global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
      global_dpd_->file2_mat_init(&D);
      global_dpd_->file2_mat_rd(&D);

      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&G, h);
	global_dpd_->buf4_mat_irrep_rd(&G, h);

	for(Gm=0; Gm < nirreps; Gm++) {
	  Ga = Gb = h^Gm;

	  for(b=0; b < virtpi[Gb]; b++) {
	    B = vir_off[Gb] + b;
	    for(a=0; a < virtpi[Ga]; a++) {
	      A = vir_off[Ga] + a;
	      for(m=0; m < occpi[Gm]; m++) {
		M = occ_off[Gm] + m;

		MB = G.params->rowidx[M][B];
		MA = G.params->colidx[M][A];

		G.matrix[h][MB][MA] += D.matrix[Ga][a][b];
	      }
	    }
	  }
	}

	global_dpd_->buf4_mat_irrep_wrt(&G, h);
	global_dpd_->buf4_mat_irrep_close(&G, h);
      }

      global_dpd_->buf4_close(&G);

      global_dpd_->file2_mat_close(&D);
      global_dpd_->file2_close(&D);

      /* Generate spin-adapted Gciab just for energy calculation */
      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
      global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gciab - Gciba", 2);
      global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 11, 5, "2 Gciab - Gciba", -1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
      global_dpd_->buf4_scmcopy(&G, PSIF_CC_GAMMA, "2 Gabcd - Gabdc", 2);
      global_dpd_->buf4_sort_axpy(&G, PSIF_CC_GAMMA, pqsr, 5, 5, "2 Gabcd - Gabdc", -1);
      global_dpd_->buf4_close(&G);
    }

  }} // namespace psi::ccdensity
