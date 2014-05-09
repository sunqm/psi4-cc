
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <liboptions/liboptions.h>
#include <psifiles.h>
#include <psi4-dec.h>
#include "ccsort/MOInfo.h"
#include "ccsort/Params.h"
#include "ccsort/Local.h"
#define EXTERN
#include "ccsort/globals.h"

namespace psi { namespace ccsort {

#define IOFF_MAX 32641

void init_io();
void init_ioff(void);
void get_params(Options & options);
void get_moinfo_light(void);
void sort_oei(void);
void sort_tei(void);
void b_sort(void);
void c_sort(void);
void d_sort(void);
void e_sort(void);
void f_sort(void);
void d_spinad(void);
void e_spinad(void);
void f_spinad(void);
void scf_check(void);
void fock(void);
void denom(void);
void exit_io(void);
void cleanup_light(void);
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhf(int **cachelist);
void local_init(Options & options);
void local_done(void);
void cc_memcheck(void);

int ccsort_light(Options &options)
{
  int i;
  int **cachelist, *cachefiles;
  int bamount,famount; /* Truncated theoretical number of B/F-type ints that
                          could be stored in cache at once                  */

  unsigned long int ia_size, ab_size, ij_size, f_size, t2_size, b_size;

  init_io();
  init_ioff();

  get_params(options);
  // fixme: are they neccessary?
  params.make_aibc = 1;
  params.make_abcd = 1;
  params.make_unpacked_abcd = 1;
  get_moinfo_light();

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 2) { /*** UHF references ***/
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
  else { /*** RHF/ROHF references ***/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    std::vector<int*> spaces;
    spaces.push_back(moinfo.occpi);
    spaces.push_back(moinfo.occ_sym);
    spaces.push_back(moinfo.virtpi);
    spaces.push_back(moinfo.vir_sym);
    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 2, spaces);
  }

  /* run a small computation of memory and disk requirements */
  cc_memcheck();

  sort_oei();
  sort_tei();
  c_sort();
  d_sort();
  e_sort();
  f_sort();
  if(params.ref == 0) {
    d_spinad();
    e_spinad();
/*     f_spinad(); */
  }
//  scf_check(); // moinfo.eref initialized here
//  psio_write_entry(PSIF_CC_INFO, "Reference Energy", (char *) &(moinfo.eref),
//           sizeof(double));
  fock();
  denom();

//  /* CPHF stuff for local correlation tests */
//  if(params.local) {
//    local_init(options); // local_init calls environment.wavefunction()
//    local_done();
//  }

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup_light();
  exit_io();
  return (Success);
}

}} //namespace psi::ccsort

