/*
 * File: psi4itrf.cc
 * interface for version PSI4.0b5
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>

#define MAIN
#undef HAVE_MPI
#include <boost/shared_ptr.hpp>
#include <psi4/psi4.h>
#include <psi4-def.h>
#include <libparallel/parallel.h>
#include <libparallel2/Communicator.h>
#include <libqt/qt.h> // for timer
#include <libmints/wavefunction.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <psiconfig.h>
#include <psifiles.h>
#include "fake_wave.h"
#include "psi4itrf.h"

#define DEF_MAXMEM (1<<30)  // 1024 M bytes

namespace psi {
    extern int psi_start(const char *outfile);
    extern int psi_stop(FILE* infile, std::string OutFileRMR, char* psi_file_prefix);
    extern int read_options(const std::string &name, Options & options, bool suppress_printing);

    namespace ccenergy { PsiReturnType ccenergy_light(Options &options); }
    namespace ccsort { int ccsort_light(Options &options); }
    namespace cchbar { PsiReturnType cchbar_light(Options &options); }
    namespace cclambda { PsiReturnType cclambda_light(Options &options); }
    namespace ccdensity { PsiReturnType ccdensity_light(Options &options); }
    namespace cctriples { PsiReturnType cctriples_light(Options &options); }
}

using namespace psi;

void psi4itrf_init_env(const char *outfile, const char *tmpdir, 
                       unsigned long max_memory,
                       int argc4MPI, char **argv4MPI)
{
    char *fake_argv[] = {NULL};
    Process::arguments.initialize(0, fake_argv);
    //Process::environment.initialize();
#if defined HAVE_MPI
    WorldComm = boost::shared_ptr<LibParallel::ParallelEnvironment>(new LibParallel::ParallelEnvironment(argc4MPI, argv4MPI));
#else
    WorldComm = boost::shared_ptr<LibParallel::ParallelEnvironment>(new LibParallel::ParallelEnvironment(0, fake_argv));
#endif
    timer_init();
    Wavefunction::initialize_singletons();
    Process::environment.set_memory(max_memory);

    assert(psi_start(outfile) != PSI_RETURN_FAILURE);

    psio_init();
    //_default_psio_manager_->set_default_path(std::string(tmpdir));
    PSIOManager::shared_object()->set_default_path(std::string(tmpdir));
}

void psi4itrf_del_env()
{
#if defined HAVE_MPI
    // Shut things down:
    WorldComm->sync();
#endif
    // There is only one timer:
    timer_done();

    PSIOManager::shared_object()->psiclean();
    psi_stop(infile, "outfile", psi_file_prefix);

#if defined HAVE_MPI
    WorldComm->sync();
    WorldComm->finalize();
#endif

    Process::environment.wavefunction().reset();

    std::cout << std::flush;
}

/*
 * NOTE: call me whenever the number of basis changes, to make sure
 * there is no confliction for the temporary files due to the change of
 * file size.
 */
void psio_clean()
{
    PSIOManager::shared_object()->psiclean();
}



/*
 * ref_wfn can be one of "RHF", "ROHF", "UHF"
 */
double psi4ccsd_energy(char *ref_wfn)
{
    double energy;

    // ccsort must be carried out before everything else, it initilizes psio
    Process::environment.options.set_current_module("CCSORT");
    read_options("CCSORT", Process::environment.options, false);
    Process::environment.options.set_int("CCSORT", "PRINT", 0);
    Process::environment.options.set_str("CCSORT", "WFN", "CCSD");
    Process::environment.options.set_str("CCSORT", "REFERENCE", ref_wfn);
    Process::environment.options.set_int("CCSORT", "CACHELEVEL", CACHELEVEL);
    Process::environment.options.set_bool("CCSORT", "KEEP_OEIFILE", true);
    Process::environment.options.set_bool("CCSORT", "KEEP_TEIFILE", true);
    ccsort::ccsort_light(Process::environment.options);

    Process::environment.options.set_current_module("CCENERGY");
    read_options("CCENERGY", Process::environment.options, false);
    Process::environment.options.set_int("CCENERGY", "PRINT", 0);
    Process::environment.options.set_str("CCENERGY", "WFN", "CCSD");
    Process::environment.options.set_str("CCENERGY", "REFERENCE", ref_wfn);
    //! options.set_bool("CCENERGY", "SEMICANONICAL", false);// will cause params.ref = 2
    Process::environment.options.set_bool("CCENERGY", "ANALYZE", false); // analyze T2 amplitudes

    Process::environment.options.set_bool("CCENERGY", "RESTART", false);
    Process::environment.options.set_str("CCENERGY", "AO_BASIS", "NONE");
    // for cachelist[pqnum][rsnum]
    // 4 ~ cache everything. decrease me if catch memeory problems
    Process::environment.options.set_int("CCENERGY", "CACHELEVEL", CACHELEVEL);
    Process::environment.options.set_bool("CCENERGY", "LOCAL", false);  // local-CC
    Process::environment.options.set_int("CCENERGY", "MAXITER", 1000);

    PsiReturnType ccsd_return = ccenergy::ccenergy_light(Process::environment.options);
    //if (ccsd_return == Success) {
    //    // Get the total energy of the CCSD wavefunction
    //    energy = Process::environment.globals["CURRENT ENERGY"];
    //}
    energy = Process::environment.globals["CURRENT ENERGY"];

    std::cout << std::flush;
    return energy;
}


void psi4cc_density()
{
    Process::environment.options.set_current_module("CCSORT");
    std::string ccname = Process::environment.options.get_str("WFN");
    std::string ref_wfn = Process::environment.options.get_str("REFERENCE");

    Process::environment.options.set_current_module("CCHBAR");
    read_options("CCHBAR", Process::environment.options, false);
    Process::environment.options.set_int("CCHBAR", "PRINT", 0);
    Process::environment.options.set_str("CCHBAR", "WFN", ccname);
    cchbar::cchbar_light(Process::environment.options);

    Process::environment.options.set_current_module("CCLAMBDA");
    read_options("CCLAMBDA", Process::environment.options, false);
    Process::environment.options.set_int("CCLAMBDA", "PRINT", 0);
    Process::environment.options.set_str("CCLAMBDA", "WFN", ccname);
    Process::environment.options.set_str("CCLAMBDA", "AO_BASIS", "NONE");
    Process::environment.options.set_int("CCLAMBDA", "MAXITER", 1000);
    cclambda::cclambda_light(Process::environment.options);

    Process::environment.options.set_current_module("CCDENSITY");
    read_options("CCDENSITY", Process::environment.options, false);
    Process::environment.options.set_int("CCDENSITY", "PRINT", 0);
    Process::environment.options.set_str("CCDENSITY", "WFN", ccname);
    Process::environment.options.set_str("CCDENSITY", "REFERENCE", ref_wfn);
    // true ~ produce onepdm only
    Process::environment.options.set_bool("CCDENSITY", "ONEPDM", false);
    Process::environment.options.set_bool("CCDENSITY", "XI", false);
    Process::environment.options.set_bool("CCDENSITY", "ZETA", false);
    Process::environment.options.set_bool("CCDENSITY", "OPDM_RELAX", true);
    ccdensity::ccdensity_light(Process::environment.options);

    std::cout << std::flush;
}


// CCSD(T)
double psi4ccsd_t(char *ref_wfn)
{
    double energy;
    psi4ccsd_energy(ref_wfn);
    Process::environment.options.set_str("CCENERGY", "WFN", "CCSD_T");

    //todo: Make sure ccenergy returned Success
    //if (ccsd_return != Success)
    //    throw PSIEXCEPTION("CCEnergyWavefunction: CCSD did not converge, will not proceed to (T) correction.");

    // Run cctriples
    //if (psi::cctriples::cctriples_light(Process::environment.options) == Success) {
    //    energy = Process::environment.globals["CURRENT ENERGY"];
    //}
    energy = Process::environment.globals["CURRENT ENERGY"];

    std::cout << std::flush;
    return energy;
}

// also see doc/progman/sample-codes/df-mp2/main.cc
