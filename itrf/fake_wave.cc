#include <boost/shared_ptr.hpp>
#include <liboptions/liboptions.h>
#include "libmints/basisset.h"
#include "libmints/molecule.h"
#include "libmints/wavefunction.h"
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include "fake_wave.h"

/*
 * The fake classes doesn't work because the base class constructor
 * assumes many initial values we cannot provide.
 */

namespace psi {

FakeRefWavefunction::FakeRefWavefunction(Options& options)
    : Wavefunction(options)
{
    chkpt_init(PSIO_OPEN_OLD);
    energy_ = chkpt_rd_escf();
    chkpt_close();
}


FakeCCWavefunction::FakeCCWavefunction(Options& options)
    : Wavefunction(options)
{
    chkpt_init(PSIO_OPEN_OLD);
    nirreps_ = 1//chkpt_rd_nirreps();
    nmo_ = chkpt_rd_nmo();
    nso_ = chkpt_rd_nso();
    //labels_ = chkpt_rd_irr_labs();
    enuc_ = chkpt_rd_enuc();
    escf_ = chkpt_rd_escf();

    nsopi_  = Dimension(nirrep_, "SOs per irrep");
    nmopi_  = Dimension(nirrep_, "MOs per irrep");
    doccpi_ = Dimension(nirrep_, "Doubly occupied orbitals per irrep");
    soccpi_ = Dimension(nirrep_, "Singly occupied orbitals per irrep");
    nalphapi_ = Dimension(nirrep_, "Alpha electrons per irrep");
    nbetapi_  = Dimension(nirrep_, "Beta electrons per irrep");
    frzcpi_ = Dimension(nirrep_, "Frozen core orbitals per irrep");
    frzvpi_ = Dimension(nirrep_, "Frozen virtual orbitals per irrep");

    int *rec_sopi = chkpt_rd_sopi();
    int *rec_orbspi = chkpt_rd_orbspi();
    int *rec_clsdpi = chkpt_rd_clsdpi();
    int *rec_openpi = chkpt_rd_openpi();
    int *rec_frzcpi = chkpt_rd_frzcpi();
    int *rec_frzvpi = chkpt_rd_frzvpi();

    for (int k = 0; k < nirrep_; k++) {
        nsopi_[k] = rec_sopi[k];
        nmopi_[k] = rec_sopi[k];
        doccpi_[k] = rec_clsdpi[k];
        soccpi_[k] = rec_openpi[k];
        nalphapi_[k] = rec_clsdpi[k] + rec_openpi[k];
        nbetapi_[k] = rec_clsdpi[k];
        frzcpi_[k] = rec_frzcpi[k];
        frzvpi_[k] = rec_frzvpi[k];
    } 

    basisset_ = FakeBasisSet(chkpt_rd_nao());
    molecule_ = FakeMolecule(chkpt_rd_enuc());

    chkpt_close();
}


FakeBasisSet::FakeBasisSet(int nao) // fixme
    : nao_(nao)
{
}



FakeMolecule::FakeMolecule(double e_nuc) // fixme
    : nuclear_repulsion_(e_nuc)
{
}

char **FakeMolecule::irrep_labels()
{
    int nirreps = 1;
    char **irreplabel = (char **) malloc(sizeof(char *)*nirreps);
    for (int i=0; i<nirreps; i++) {
        irreplabel[i] = (char *) malloc(sizeof(char)*5);
        ::memset(irreplabel[i], 0, sizeof(char)*5);
    }
    return irreplabel;
}

double FakeMolecule::nuclear_repulsion_energy()
{
    return nuclear_repulsion_;
}

}
