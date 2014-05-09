#include <cstdio>
#include <strings.h>
#include <string.h>
#include <libdpd/dpd.h>
#include "ccdensity/MOInfo.h"
#include "ccdensity/Params.h"
#include "ccdensity/Frozen.h"
#define EXTERN
#include "ccdensity/globals.h"

namespace psi { namespace ccdensity {

void fold_UHF_light(struct RHO_Params rho_params)
{
  int h, nirreps;
  int i, j, k, l, m, a, b;
  int I, J, K, L, M, A, B;
  int IM, JM, MI, MJ, MK, ML, MA, MB;
  int Gi, Gj, Gk, Gl, Gm, Ga, Gb;
  int *aoccpi, *avirtpi;
  int *boccpi, *bvirtpi;
  int *aocc_off, *avir_off;
  int *bocc_off, *bvir_off;
  int *aocc_sym, *avir_sym;
  int *bocc_sym, *bvir_sym;
  dpdfile2 D, D1, D2, F;
  dpdbuf4 G, Aints, E, C, DInts, FInts, BInts, G1, G2;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi; avirtpi = moinfo.avirtpi;
  boccpi = moinfo.boccpi; bvirtpi = moinfo.bvirtpi;
  aocc_off = moinfo.aocc_off; avir_off = moinfo.avir_off;
  bocc_off = moinfo.bocc_off; bvir_off = moinfo.bvir_off;
  aocc_sym = moinfo.aocc_sym; avir_sym = moinfo.avir_sym;
  bocc_sym = moinfo.bocc_sym; bvir_sym = moinfo.bvir_sym;

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 2, 2, 0, "GIJKL");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Gj = h^Gm;

      for(i=0; i < aoccpi[Gi]; i++) {
	I = aocc_off[Gi] + i;
	for(j=0; j < aoccpi[Gj]; j++) {
	  J = aocc_off[Gj] + j;
	  for(m=0; m < aoccpi[Gm]; m++) {
	    M = aocc_off[Gm] + m;

	    IM = G.params->rowidx[I][M];
	    JM = G.params->colidx[J][M];
	    MI = G.params->rowidx[M][I];
	    MJ = G.params->colidx[M][J];

	    G.matrix[h][IM][JM] += D.matrix[Gi][i][j];
	    G.matrix[h][IM][MJ] -= D.matrix[Gi][i][j];
	    G.matrix[h][MI][MJ] += D.matrix[Gi][i][j];
	    G.matrix[h][MI][JM] -= D.matrix[Gi][i][j];
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

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 12, 12, 0, "Gijkl");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Gj = h^Gm;

      for(i=0; i < boccpi[Gi]; i++) {
	I = bocc_off[Gi] + i;
	for(j=0; j < boccpi[Gj]; j++) {
	  J = bocc_off[Gj] + j;
	  for(m=0; m < boccpi[Gm]; m++) {
	    M = bocc_off[Gm] + m;

	    IM = G.params->rowidx[I][M];
	    JM = G.params->colidx[J][M];
	    MI = G.params->rowidx[M][I];
	    MJ = G.params->colidx[M][J];

	    G.matrix[h][IM][JM] += D.matrix[Gi][i][j];
	    G.matrix[h][MI][JM] -= D.matrix[Gi][i][j];
	    G.matrix[h][MI][MJ] += D.matrix[Gi][i][j];
	    G.matrix[h][IM][MJ] -= D.matrix[Gi][i][j];
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

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Gj = h^Gm;

      for(i=0; i < aoccpi[Gi]; i++) {
	I = aocc_off[Gi] + i;
	for(j=0; j < aoccpi[Gj]; j++) {
	  J = aocc_off[Gj] + j;
	  for(m=0; m < boccpi[Gm]; m++) {
	    M = bocc_off[Gm] + m;

	    IM = G.params->rowidx[I][M];
	    JM = G.params->colidx[J][M];

	    G.matrix[h][IM][JM] += D.matrix[Gi][i][j];
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


  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 22, 22, 22, 0, "GIjKl");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gk = Gl = h^Gm;

      for(k=0; k < boccpi[Gk]; k++) {
	K = bocc_off[Gk] + k;
	for(l=0; l < boccpi[Gl]; l++) {
	  L = bocc_off[Gl] + l;
	  for(m=0; m < aoccpi[Gm]; m++) {
	    M = aocc_off[Gm] + m;

	    MK = G.params->rowidx[M][K];
	    ML = G.params->colidx[M][L];

	    G.matrix[h][MK][ML] += D.matrix[Gk][k][l];
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

  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
  global_dpd_->file2_mat_init(&D1);
  global_dpd_->file2_mat_rd(&D1);
  global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
  global_dpd_->file2_mat_init(&D2);
  global_dpd_->file2_mat_rd(&D2);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 20, 2, 20, 0, "GIJKA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < aoccpi[Gi]; i++) {
	I = aocc_off[Gi] + i;
	for(a=0; a < avirtpi[Ga]; a++) {
	  A = avir_off[Ga] + a;
	  for(m=0; m < aoccpi[Gm]; m++) {
	    M = aocc_off[Gm] + m;

	    MI = G.params->rowidx[M][I];
	    IM = G.params->rowidx[I][M];
	    MA = G.params->colidx[M][A];

	    G.matrix[h][MI][MA] += 0.5 * (D1.matrix[Gi][i][a] +
					  D2.matrix[Gi][i][a]);
	    G.matrix[h][IM][MA] -= 0.5 * (D1.matrix[Gi][i][a] +
					  D2.matrix[Gi][i][a]);
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);
  global_dpd_->file2_mat_close(&D1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_mat_close(&D2);
  global_dpd_->file2_close(&D2);

  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
  global_dpd_->file2_mat_init(&D1);
  global_dpd_->file2_mat_rd(&D1);
  global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
  global_dpd_->file2_mat_init(&D2);
  global_dpd_->file2_mat_rd(&D2);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 30, 12, 30, 0, "Gijka");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < boccpi[Gi]; i++) {
	I = bocc_off[Gi] + i;
	for(a=0; a < bvirtpi[Ga]; a++) {
	  A = bvir_off[Ga] + a;
	  for(m=0; m < boccpi[Gm]; m++) {
	    M = bocc_off[Gm] + m;

	    MI = G.params->rowidx[M][I];
	    IM = G.params->rowidx[I][M];
	    MA = G.params->colidx[M][A];

	    G.matrix[h][MI][MA] += 0.5 * (D1.matrix[Gi][i][a] +
					  D2.matrix[Gi][i][a]);
	    G.matrix[h][IM][MA] -= 0.5 * (D1.matrix[Gi][i][a] +
					  D2.matrix[Gi][i][a]);
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);
  global_dpd_->file2_mat_close(&D1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_mat_close(&D2);
  global_dpd_->file2_close(&D2);

  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
  global_dpd_->file2_mat_init(&D1);
  global_dpd_->file2_mat_rd(&D1);
  global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
  global_dpd_->file2_mat_init(&D2);
  global_dpd_->file2_mat_rd(&D2);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < aoccpi[Gi]; i++) {
	I = aocc_off[Gi] + i;
	for(a=0; a < avirtpi[Ga]; a++) {
	  A = avir_off[Ga] + a;
	  for(m=0; m < boccpi[Gm]; m++) {
	    M = bocc_off[Gm] + m;

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
  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&D1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_mat_close(&D2);
  global_dpd_->file2_close(&D2);


  global_dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
  global_dpd_->file2_mat_init(&D1);
  global_dpd_->file2_mat_rd(&D1);
  global_dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
  global_dpd_->file2_mat_init(&D2);
  global_dpd_->file2_mat_rd(&D2);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Gi = Ga = h^Gm;

      for(i=0; i < boccpi[Gi]; i++) {
	I = bocc_off[Gi] + i;
	for(a=0; a < bvirtpi[Ga]; a++) {
	  A = bvir_off[Ga] + a;
	  for(m=0; m < aoccpi[Gm]; m++) {
	    M = aocc_off[Gm] + m;

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
  global_dpd_->buf4_close(&G);
  global_dpd_->file2_mat_close(&D1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_mat_close(&D2);
  global_dpd_->file2_close(&D2);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 20, 20, 20, 20, 0, "GIBJA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < avirtpi[Gb]; b++) {
	B = avir_off[Gb] + b;
	for(a=0; a < avirtpi[Ga]; a++) {
	  A = avir_off[Ga] + a;
	  for(m=0; m < aoccpi[Gm]; m++) {
	    M = aocc_off[Gm] + m;

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


  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 30, 30, 30, 30, 0, "Gibja");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < bvirtpi[Gb]; b++) {
	B = bvir_off[Gb] + b;
	for(a=0; a < bvirtpi[Ga]; a++) {
	  A = bvir_off[Ga] + a;
	  for(m=0; m < boccpi[Gm]; m++) {
	    M = bocc_off[Gm] + m;

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


  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 24, 24, 24, 24, 0, "GIbJa");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < bvirtpi[Gb]; b++) {
	B = bvir_off[Gb] + b;
	for(a=0; a < bvirtpi[Ga]; a++) {
	  A = bvir_off[Ga] + a;
	  for(m=0; m < aoccpi[Gm]; m++) {
	    M = aocc_off[Gm] + m;

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

  global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 27, 27, 27, 27, 0, "GiBjA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);

    for(Gm=0; Gm < nirreps; Gm++) {
      Ga = Gb = h^Gm;

      for(b=0; b < avirtpi[Gb]; b++) {
	B = avir_off[Gb] + b;
	for(a=0; a < avirtpi[Ga]; a++) {
	  A = avir_off[Ga] + a;
	  for(m=0; m < boccpi[Gm]; m++) {
	    M = bocc_off[Gm] + m;

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
}

}} // namespace psi::ccdensity
