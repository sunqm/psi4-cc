/*!
** \file
** \brief Initialize input, output, file prefix, etc.
** \ingroup
* modified by sqm
*/

#include <psi4/psi4.h>
#include <psi4-dec.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <getopt.h>
#include <psifiles.h>
#include <psiconfig.h>
#include <libplugin/plugin.h>
#include <libparallel/parallel.h>
#include <libparallel/ParallelPrinter.h>

using namespace std;

namespace psi {

int psi_start(const char *file6)
{
    std::string ofname  = file6;
    std::string fprefix = PSI_DEFAULT_FILE_PREFIX;
    check_only          = false;
    clean_only          = false;
    verbose             = false;

    Process::environment.set_n_threads(1);
    messy = true;

    if(ofname == "stdout"){
        outfile = boost::shared_ptr<PsiOutStream>(new PsiOutStream());
    } else {
        outfile = boost::shared_ptr<PsiOutStream>(new OutFile(ofname,TRUNCATE));
    }

    /* copy over file prefix, etc. into their appropriate variables */
    psi_file_prefix = strdup(fprefix.c_str());

    outfile_name = ofname;

    return(PSI_RETURN_SUCCESS);
}
}

