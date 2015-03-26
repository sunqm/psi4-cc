#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <psifiles.h>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libparallel/ParallelPrinter.h>

namespace psi {

int psi_stop(FILE* infile, std::string OutFileRMR, char* psi_file_prefix)
{
  delete psi_file_prefix;
 
  outfile->Printf( "\n*** psi4cc successful. Buy Sebastian a beer! ***\n");
  infile = NULL;
  outfile = boost::shared_ptr<OutFile>();
  psi_file_prefix = NULL;

  return(PSI_RETURN_SUCCESS);
}

}

