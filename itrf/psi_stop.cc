#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <psifiles.h>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
namespace psi {

int psi_stop(FILE* infile, FILE* outfile, char* psi_file_prefix)
{
  delete psi_file_prefix;
 
  //fclose(infile);
  if (outfile != stdout) {
      fflush(outfile);
      fclose(outfile);
  }

  infile = NULL;
  outfile = NULL;
  psi_file_prefix = NULL;

  return(PSI_RETURN_SUCCESS);
}

}

