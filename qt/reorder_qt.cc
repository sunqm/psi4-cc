/*

 PSI4: an ab initio quantum chemistry software package

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/

/********************************************************************************
    Modified by Sebastian Wouters from
    https://github.com/psi4/psi4public/blob/master/src/lib/libqt/reorder_qt.cc
    so that no wavefunction object is needed anymore!
********************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libmints/matrix.h>
#include "psi4-dec.h"

namespace psi {

/*!
** reorder_qt_uhf()
**
** Generalization of reorder_qt() for UHF case
**
** \param docc        = doubly occupied orbitals per irrep
** \param socc        = singly occupied orbitals per irrep
** \param frozen_docc = frozen occupied orbitals per irrep
** \param frozen_uocc = frozen unoccupied orbitals per irrep
** \param order_alpha = reordering array for alpha (Pitzer->QT order)
** \param order_beta  = reordering array for beta  (Pitzer->QT order)
** \param nirreps     = number of irreducible representations
**
** \ingroup QT
*/
void reorder_qt_uhf_modified(int *docc, int *socc, int *frozen_docc,
                    int *frozen_uocc, int *order_alpha, int *order_beta,
                    int *orbspi, int nirreps)
{
  int p, nmo;
  int cnt_alpha, cnt_beta, irrep, tmpi;
  int *offset, this_offset;
  int *uocc;

  int * nalphapi = new int[nirreps];
  int * nbetapi  = new int[nirreps];
  for (int h=0; h<nirreps; h++){
      nalphapi[h] = frozen_docc[h] + docc[h] + socc[h];
      nbetapi[h]  = frozen_docc[h] + docc[h];
  }

  offset = init_int_array(nirreps);

  uocc = init_int_array(nirreps);

  /* construct the offset array */
  offset[0] = 0;
  for (irrep=1; irrep<nirreps; irrep++) {
    offset[irrep] = offset[irrep-1] + orbspi[irrep-1];
  }

  /* construct the uocc array */
  nmo = 0;
  for (irrep=0; irrep<nirreps; irrep++) {
    nmo += orbspi[irrep];
    tmpi = frozen_uocc[irrep] + docc[irrep] + socc[irrep];
    if (tmpi > orbspi[irrep]) {
      outfile->Printf( "(reorder_qt_uhf): orbitals don't add up for irrep %d\n",
              irrep);
      return;
    }
    else
      uocc[irrep] = orbspi[irrep] - tmpi;
  }

  cnt_alpha = cnt_beta = 0;

  /* do the frozen core */
  for(irrep=0; irrep<nirreps; irrep++) {
    this_offset = offset[irrep];
    for(p=0; p < frozen_docc[irrep]; p++) {
      order_alpha[this_offset+p] = cnt_alpha++;
      order_beta[this_offset+p] = cnt_beta++;
    }
  }

  /* alpha occupied orbitals */
  for(irrep=0; irrep<nirreps; irrep++) {
    this_offset = offset[irrep] + frozen_docc[irrep];
    for(p=0; p < nalphapi[irrep] - frozen_docc[irrep]; p++) {
      order_alpha[this_offset + p] = cnt_alpha++;
    }
  }

  /* beta occupied orbitals */
  for(irrep=0; irrep<nirreps; irrep++) {
    this_offset = offset[irrep] + frozen_docc[irrep];
    for(p=0; p < nbetapi[irrep] - frozen_docc[irrep]; p++) {
      order_beta[this_offset + p] = cnt_beta++;
    }
  }

  /* alpha unoccupied orbitals */
  for(irrep=0; irrep<nirreps; irrep++) {
    this_offset = offset[irrep] + nalphapi[irrep];
    for(p=0; p < orbspi[irrep] - nalphapi[irrep]; p++) {
      order_alpha[this_offset + p] = cnt_alpha++;
    }
  }

  /* beta unoccupied orbitals */
  for(irrep=0; irrep<nirreps; irrep++) {
    this_offset = offset[irrep] + nbetapi[irrep];
    for(p=0; p < orbspi[irrep] - nbetapi[irrep]; p++) {
      order_beta[this_offset + p] = cnt_beta++;
    }
  }

  /* do the frozen uocc */
  for (irrep=0; irrep<nirreps; irrep++) {
    this_offset = offset[irrep] + docc[irrep] + socc[irrep] + uocc[irrep];
    for(p=0; p < frozen_uocc[irrep]; p++) {
      order_alpha[this_offset + p] = cnt_alpha++;
      order_beta[this_offset + p] = cnt_beta++;
    }
  }

  /* do a final check */
  for (irrep=0; irrep<nirreps; irrep++) {
    if (cnt_alpha > nmo) {
      outfile->Printf( "(reorder_qt_uhf): on final check, used more orbitals");
      outfile->Printf( "   than were available (%d vs %d) for irrep %d\n",
              cnt_alpha, nmo, irrep);
    }
    if (cnt_beta > nmo) {
      outfile->Printf( "(reorder_qt_uhf): on final check, used more orbitals");
      outfile->Printf( "   than were available (%d vs %d) for irrep %d\n",
              cnt_beta, nmo, irrep);
    }
  }

  free(offset);
  free(uocc);
  delete [] nalphapi;
  delete [] nbetapi;
}

}

