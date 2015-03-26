#include <cstdio>
#include <cstdlib>
#include <string>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include <psi4-dec.h>
#include "cchbar/MOInfo.h"
#include "cchbar/Params.h"
#define EXTERN
#include "cchbar/globals.h"

namespace psi { namespace cchbar {

void init_io();
void title(void);
void get_moinfo_light(Options &);
void get_params(Options &);
void exit_io(void);
void F_build(void);
void Wmbej_build(void);
void Wmnie_build(void);
void Wmbij_build(void);
void Wabij_build(void);
void Wamef_build(void);
void Wabei_build(void);
void cc2_Zmbej_build(void);
void cc2_Wmbej_build(void);
void cc2_Wmbij_build(void);
void cc2_Wabei_build(void);
void purge(void);
void cleanup_light(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhf(int **cachelist);
void sort_amps(void);
void tau_build(void);
void taut_build(void);
void status(const char *, std::string );
void cc3_HET1(void);
void Fai_build(void);
void reference(void);
void norm_HET1(void);

using namespace psi;

PsiReturnType cchbar_light(Options &options)
{
  int **cachelist, *cachefiles;

  init_io();
  title();
  get_moinfo_light(options);
  get_params(options);

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    std::vector<int*> spaces;
    spaces.push_back(moinfo.occpi);
    spaces.push_back(moinfo.occ_sym);
    spaces.push_back(moinfo.virtpi);
    spaces.push_back(moinfo.vir_sym);
    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 2, spaces);
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

  sort_amps();
  tau_build();
  taut_build();

  if (params.Tamplitude) {
    reference();
    Fai_build();
    Wabij_build();
  }

  F_build();
  if(params.print & 2) status("F elements", "outfile");

  Wamef_build();
  if(params.print & 2) status("Wamef elements", "outfile");
  Wmnie_build();
  if(params.print & 2) status("Wmnie elements", "outfile");

  if(params.wfn == "CC2" || params.wfn == "EOM_CC2") {
    cc2_Wmbej_build();
    if(params.print & 2) status("Wmbej elements", "outfile");
    cc2_Zmbej_build();
    if(params.print & 2) status("Zmbej elements", "outfile");
    cc2_Wmbij_build();
    if(params.print & 2) status("Wmbij elements", "outfile");
    cc2_Wabei_build();
    if(params.print & 2) status("Wabei elements", "outfile");
  }
  else {
    Wabei_build();
    if(params.print & 2) status("Wabei elements", "outfile");
    Wmbej_build();
    if(params.print & 2) status("Wmbej elements", "outfile");
    Wmbij_build();
    if(params.print & 2) status("Wmbij elements", "outfile");

    if( params.wfn == "CC3" || params.wfn == "EOM_CC3" ) {
      /* switch to ROHF to generate all spin cases of He^T1 elements */
      if((params.dertype == 3 || params.dertype == 1) && params.ref == 0) {
        params.ref = 1;
        cc3_HET1(); /* compute remaining Wmbej [H,eT1] */
        norm_HET1();
        params.ref = 0;
      }
      else {
        cc3_HET1(); /* compute remaining Wmbej [H,eT1] */
        norm_HET1();
      }
    }
  }

  if(params.ref == 1) purge(); /** ROHF only **/
  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup_light();
  exit_io();
  return Success;
}

}} // namespace psi::chbar
