#include <boost/shared_ptr.hpp>
#include <liboptions/liboptions.h>
#include "libmints/basisset.h"
#include "libmints/molecule.h"
#include "libmints/wavefunction.h"

/*
 * Fake reference WF and CC WF to mimic the PSI4 default one because
 * Wavefunction constructor calls Wavefunction::common_init which
 * requires more parameters we don't initialized
 */
namespace psi {
class FakeRefWavefunction : public Wavefunction {
    public:
        FakeRefWavefunction(Options &options);
        ~FakeRefWavefunction();
        double compute_energy() { return 0; }
};


class FakeCCWavefunction : public Wavefunction {
    public:
        FakeCCWavefunction(Options &options);
        ~FakeCCWavefunction();
        double compute_energy() { return 0; }
};



class FakeBasisSet : public BasisSet {
    public:
        FakeBasisSet(int nao);
        ~FakeBasisSet();
};

class FakeMolecule : public Molecule {
    public:
        FakeMolecule(double e_nuc);
        ~FakeMolecule();

    private:
        double nuclear_repulsion_;
};

}
