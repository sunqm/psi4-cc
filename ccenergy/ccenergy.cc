/*
 * modified from ccenergy.cc
 * remove unused calling
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <libint/libint.h>
#include <libmints/wavefunction.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <sys/types.h>
#include <psifiles.h>
#include "ccenergy/Params.h"
#include "ccenergy/MOInfo.h"
#include "ccenergy/Local.h"
#define EXTERN
#include "ccenergy/globals.h"
#include "ccenergy/ccwave.h"

namespace psi { namespace ccenergy {

void init_io();
void title(void);
void init_ioff(void);
void get_moinfo_light(void);
void get_params(Options &);
void init_amps(void);
void tau_build(void);
void taut_build(void);
double energy(void);
double mp2_energy(void);
void sort_amps(void);
void Fae_build(void);
void Fmi_build(void);
void Fme_build(void);
void t1_build(void);
void Wmnij_build(void);
void Z_build(void);
void Wmbej_build(void);
void t2_build(void);
void tsave(void);
int converged(double);
double diagnostic(void);
double d1diag(void);
double new_d1diag(void);
double d2diag(void);
void exit_io(void);
void cleanup_light(void);
void update(void);
void diis(int iter);
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
struct dpd_file4_cache_entry *priority_list(void);
void spinad_amps(void);
void status(const char *, FILE *);
void lmp2(void);
void amp_write(void);
int rotate(void);
void analyze(void);
void cc3_Wmnie(void);
void cc3_Wamef(void);
void cc3_Wmnij(void);
void cc3_Wmbij(void);
void cc3_Wabei(void);
void cc3(void);
void cc2_Wmnij_build(void);
void cc2_Wmbij_build(void);
void cc2_Wabei_build(void);
void cc2_t2_build(void);
void one_step(void);
void denom(void);
void pair_energies(double** epair_aa, double** epair_ab);
void print_pair_energies(double* emp2_aa, double* emp2_ab, double* ecc_aa,
                         double* ecc_ab);
void checkpoint(void);
void form_df_ints(Options &options, int **cachelist, int *cachefiles, dpd_file4_cache_entry *priority);

/* local correlation functions */
void local_init(void);
void local_done(void);
}} //namespace psi::ccenergy

//TODO: namespace psi { namespace cctriples {
//TODO: PsiReturnType cctriples(Options &options);
//TODO: }}

namespace psi { namespace ccenergy {

PsiReturnType ccenergy_light(Options &options)
{
    int done=0, brueckner_done=0;
    int h, i, j, a, b, row, col, natom;
    double **geom, *zvals, value;
    FILE *efile;
    int **cachelist, *cachefiles;
    struct dpd_file4_cache_entry *priority;
    dpdfile2 t1;
    dpdbuf4 t2;
    double *emp2_aa, *emp2_ab, *ecc_aa, *ecc_ab, tval;

    moinfo.iter=0;

    init_io();
    title();
    init_ioff();

    get_moinfo_light();
    get_params(options);

    cachefiles = init_int_array(PSIO_MAXUNIT);

    if(params.ref == 2) { /** UHF **/
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
        delete[] dpd_list[0];
        dpd_list[0] = new DPD(0, moinfo.nirreps, params.memory, 0, cachefiles,
                              cachelist, NULL, 4, spaces);
        dpd_set_default(0);

        if( params.df ){
            //form_df_ints(options, cachelist, cachefiles, priority);
            // don't do form_df_ints since it needs environment.wavefunction()
        }else if( params.aobasis != "NONE" ) { /* Set up new DPD's for AO-basis algorithm */
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
    else { /** RHF or ROHF **/
        cachelist = cacheprep_rhf(params.cachelev, cachefiles);

        priority = priority_list();
        std::vector<int*> spaces;
        spaces.push_back(moinfo.occpi);
        spaces.push_back(moinfo.occ_sym);
        spaces.push_back(moinfo.virtpi);
        spaces.push_back(moinfo.vir_sym);

        dpd_init(0, moinfo.nirreps, params.memory, params.cachetype, cachefiles, cachelist, priority, 2, spaces);

        if( params.df ){
            //form_df_ints(options, cachelist, cachefiles, priority);
            // don't do form_df_ints since it needs environment.wavefunction()
        }else if( params.aobasis != "NONE") { /* Set up new DPD for AO-basis algorithm */
            std::vector<int*> aospaces;
            aospaces.push_back(moinfo.occpi);
            aospaces.push_back(moinfo.occ_sym);
            aospaces.push_back(moinfo.sopi);
            aospaces.push_back(moinfo.sosym);
            dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 2, aospaces);
            dpd_set_default(0);
        }

    }

    if ( (params.just_energy) || (params.just_residuals) ) {
        one_step();
        if(params.ref == 2) cachedone_uhf(cachelist); else cachedone_rhf(cachelist);
        free(cachefiles);
        cleanup_light();
        exit_io();
        return Success;
    }

    if(params.local) {
        local_init();
        if(local.weakp=="MP2") lmp2();
    }

    init_amps();
//
//    /* Compute the MP2 energy while we're here */
//    if(params.ref == 0 || params.ref == 2) {
//        moinfo.emp2 = mp2_energy();
//        psio_write_entry(PSIF_CC_INFO, "MP2 Energy", (char *) &(moinfo.emp2),sizeof(double));
//        Process::environment.globals["MP2 CORRELATION ENERGY"] = moinfo.emp2;
//        Process::environment.globals["MP2 TOTAL ENERGY"] = moinfo.emp2 + moinfo.eref;
//    }
//
//    if(params.print_mp2_amps) amp_write();

    tau_build();
    taut_build();
    fprintf(outfile, "\t            Solving CC Amplitude Equations\n");
    fprintf(outfile, "\t            ------------------------------\n");
    fprintf(outfile, "  Iter             Energy              RMS        T1Diag      D1Diag    New D1Diag    D2Diag\n");
    fprintf(outfile, "  ----     ---------------------    ---------   ----------  ----------  ----------   --------\n");
    moinfo.ecc = energy(); // fixme
    pair_energies(&emp2_aa, &emp2_ab);
    double last_energy;

    moinfo.t1diag = diagnostic();
    moinfo.d1diag = d1diag();
    moinfo.new_d1diag = new_d1diag();
    fflush(outfile);
    moinfo.d2diag = d2diag();
    update();
    checkpoint();
    for(moinfo.iter=1; moinfo.iter <= params.maxiter; moinfo.iter++) {

        sort_amps();
        Fme_build(); Fae_build(); Fmi_build();
        if(params.print & 2) status("F intermediates", outfile);

        t1_build();
        if(params.print & 2) status("T1 amplitudes", outfile);

        if( params.wfn == "CC2"  || params.wfn == "EOM_CC2" ) {

            cc2_Wmnij_build();
            if(params.print & 2) status("Wmnij", outfile);

            cc2_Wmbij_build();
            if(params.print & 2) status("Wmbij", outfile);

            cc2_Wabei_build();
            if(params.print & 2) status("Wabei", outfile);

            cc2_t2_build();
            if(params.print & 2) status("T2 amplitudes", outfile);

        }

        else {

            Wmbej_build();
            if(params.print & 2) status("Wmbej", outfile);

            Z_build();
            if(params.print & 2) status("Z", outfile);
            Wmnij_build();
            if(params.print & 2) status("Wmnij", outfile);

            t2_build();
            if(params.print & 2) status("T2 amplitudes", outfile);

            if( params.wfn == "CC3" || params.wfn == "EOM_CC3" ) {

                /* step1: build cc3 intermediates, Wabei, Wmnie, Wmbij, Wamef */
                cc3_Wmnij();
                cc3_Wmbij();
                cc3_Wmnie();
                cc3_Wamef();
                cc3_Wabei();

                /* step2: loop over T3's and add contributions to T1 and T2 as you go */
                cc3();
            }
        }

        if (!params.just_residuals)
            denom(); /* apply denominators to T1 and T2 */

        if(converged(last_energy - moinfo.ecc)) {
            done = 1;

            tsave();
            tau_build(); taut_build();
            last_energy = moinfo.ecc;
            moinfo.ecc = energy();
            moinfo.t1diag = diagnostic();
            moinfo.d1diag = d1diag();
            moinfo.new_d1diag = new_d1diag();
            moinfo.d2diag = d2diag();
            sort_amps();
            update();
            fprintf(outfile, "\n\tIterations converged.\n");
            fflush(outfile);
            fprintf(outfile, "\n");
            amp_write();
            if (params.analyze != 0) analyze();
            break;
        }
        if(params.diis) diis(moinfo.iter);
        tsave();
        tau_build(); taut_build();
        last_energy = moinfo.ecc;
        moinfo.ecc = energy();
        moinfo.t1diag = diagnostic();
        moinfo.d1diag = d1diag();
        moinfo.new_d1diag = new_d1diag();
        moinfo.d2diag = d2diag();
        update();
        checkpoint();
    }  // end loop over iterations

    if(!done) {
        fprintf(outfile, "\t ** Wave function not converged to %2.1e ** \n",
                params.convergence);
        fflush(outfile);
        if( params.aobasis != "NONE" ) dpd_close(1);
        dpd_close(0);
        cleanup_light();
        exit_io();
        return Failure;
    }

    //Process::environment.globals["SCF TOTAL ENERGY (CHKPT)"] = moinfo.escf;
    //Process::environment.globals["SCF TOTAL ENERGY"] = moinfo.eref;

    /* Write total energy to the checkpoint file */
    chkpt_init(PSIO_OPEN_OLD);
    chkpt_wt_etot(moinfo.ecc+moinfo.eref);
    chkpt_close();

    /* Generate the spin-adapted RHF amplitudes for later codes */
    if(params.ref == 0) spinad_amps();

    /* Compute pair energies */
    if(params.print_pair_energies) {
        pair_energies(&ecc_aa, &ecc_ab);
        print_pair_energies(emp2_aa, emp2_ab, ecc_aa, ecc_ab);
    }

    if( (params.wfn == "CC3" || params.wfn == "EOM_CC3" )
            && (params.dertype == 1 || params.dertype == 3) && params.ref == 0) {
        params.ref = 1;
        /* generate the ROHF versions of the He^T1 intermediates */
        cc3_Wmnij();
        cc3_Wmbij();
        cc3_Wmnie();
        cc3_Wamef();
        cc3_Wabei();
        //    params.ref == 0;
    }

    if(params.local) {
        /*    local_print_T1_norm(); */
        local_done();
    }

    //if(params.brueckner)
    //    Process::environment.globals["BRUECKNER CONVERGED"] = rotate();
    //don't call rotate as it needs environment.wavefunction()

    if( params.aobasis != "NONE" ) dpd_close(1);
    dpd_close(0);

    if(params.ref == 2) cachedone_uhf(cachelist);
    else cachedone_rhf(cachelist);
    free(cachefiles);

    cleanup_light();

    Process::environment.globals["CURRENT ENERGY"] = moinfo.ecc+moinfo.eref;
    Process::environment.globals["CURRENT CORRELATION ENERGY"] = moinfo.ecc;
    //Process::environment.globals["CC TOTAL ENERGY"] = moinfo.ecc+moinfo.eref;
    //Process::environment.globals["CC CORRELATION ENERGY"] = moinfo.ecc;

    exit_io();
    //  if(params.brueckner && brueckner_done)
    //     throw FeatureNotImplemented("CCENERGY", "Brueckner end loop", __FILE__, __LINE__);
    //else
    return Success;
}

}} //namespace psi::ccenergy

