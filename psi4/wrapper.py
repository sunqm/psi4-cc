import os, sys
import shutil
import tempfile
import numpy
import _psi4

class capture_stdout:
    '''redirect stdout to a string

    Examples
    --------
    with capture_stdout() as stdout:
        psi4.function
    print stdout.read()
    '''
    def __enter__(self):
        sys.stdout.flush()
        self._contents = None
        self.old_stdout_fileno = sys.stdout.fileno()
        self.bak_stdout_fd = os.dup(self.old_stdout_fileno)
        self.fd, self.ftmp = tempfile.mkstemp(prefix='tmpsi4')
        os.dup2(self.fd, self.old_stdout_fileno)
        return self
    def __exit__(self, type, value, traceback):
        sys.stdout.flush()
        self._contents = open(self.ftmp, 'r').read()
        os.dup2(self.bak_stdout_fd, self.old_stdout_fileno)
        os.close(self.fd)
        os.remove(self.ftmp)

    def read(self):
        if self._contents:
            return self._contents
        else:
            sys.stdout.flush()
            return open(self.ftmp, 'r').read()

class quite_run:
    '''output nothing

    Examples
    --------
    with quite_run():
        psi4.function
    '''
    def __enter__(self):
        sys.stdout.flush()
        self.dirnow = os.getcwd()
        self.tmpdir = tempfile.mkdtemp(prefix='tmpsi4')
        os.chdir(self.tmpdir)
        self.old_stdout_fileno = sys.stdout.fileno()
        self.bak_stdout_fd = os.dup(self.old_stdout_fileno)
        self.fnull = open(os.devnull, 'wb')
        os.dup2(self.fnull.fileno(), self.old_stdout_fileno)
    def __exit__(self, type, value, traceback):
        sys.stdout.flush()
        os.dup2(self.bak_stdout_fd, self.old_stdout_fileno)
        self.fnull.close()
        shutil.rmtree(self.tmpdir)
        os.chdir(self.dirnow)


class Solver(object):
    __all__ = ['set_psi4_env', 'del_psi4_env', 'prepare_chkpt', 'read_rdm',
               'energy', 'density', 'capture_stdout', 'quite_run',]
    '''
    Interface for Psi4 CC solver

    Examples
    --------
    ps = psi4.Solver
    with psi4.quite_run():
        ps.prepare(mo_coeff, fock_on_mo, hcore_on_mo, eri_on_mo, nelec)
        ecc = ps.energy('CCSD')
        rdm1, rdm2 = ps.density()
    '''

    instances = 0
    def __init__(self, max_memory=(1<<30)):
        self._tmpdir = tempfile.mkdtemp(prefix='tmpsi4')
        if Solver.instances == 0:
            _psi4.set_psi4_env(self._tmpdir, max_memory=max_memory)
        Solver.instances += 1

        self.nmo = 0
        self.ref = 'RHF'

    def __del__(self):
        Solver.instances -= 1
        if Solver.instances == 0:
            _psi4.del_psi4_env()
        shutil.rmtree(self._tmpdir)

    def prepare(self, ref, mo, hcore_mo, eri_mo, nelec, ms=0):
        '''initialize psi4 chkpt and psio files

        Parameters
        ----------
        ref : {'RHF', 'UHF'}
            reference wave functions for CC solver
        mo : array
            MO coefficients
        fock_mo : array
            Fock matrix
        hcore_mo : array or (array_alpha, array_beta)
            if ref is UHF, it is a tuple of arrays
        eri_mo :  array or (array_aa, array_bb, array_ab)
            if ref is UHF, it is a tuple of arrays
        nelec : int
            number of electrons

        Returns
        -------
        None
        '''
        _psi4.psioclean_before_preparing()

        assert(ref.upper() in ('RHF', 'ROHF', 'UHF'))
        self.ref = ref.upper()

        if self.ref == 'RHF':
            self.nmo = mo.shape[1]
            mo_energy = numpy.zeros(self.nmo)
            for i in range(self.nmo):
                mo_energy[i] = hcore_mo[i,i]
                ii = i*(i+1)/2+i
                for j in range(nelec/2):
                    jj = j*(j+1)/2+j
                    if i < j:
                        ij = j*(j+1)/2+i
                    else:
                        ij = i*(i+1)/2+j
                    mo_energy[i] += 2*eri_mo[ii,jj] - eri_mo[ij,ij]
        else: #UHF
            self.nmo = mo[0].shape[1]
            mo_energy = numpy.zeros((2,self.nmo))
            for i in range(self.nmo):
                mo_energy[0,i] = hcore_mo[0][i,i]
                mo_energy[1,i] = hcore_mo[1][i,i]
                ii = i*(i+1)/2+i
                for j in range((nelec+ms)/2): # alpha
                    jj = j*(j+1)/2+j
                    if i < j:
                        ij = j*(j+1)/2+i
                    else:
                        ij = i*(i+1)/2+j
                    mo_energy[0][i] += eri_mo[0][ii,jj] - eri_mo[0][ij,ij]
                    mo_energy[1][i] += eri_mo[2][jj,ii]
                for j in range((nelec-ms)/2): # beta
                    jj = j*(j+1)/2+j
                    if i < j:
                        ij = j*(j+1)/2+i
                    else:
                        ij = i*(i+1)/2+j
                    mo_energy[0][i] += eri_mo[2][ii,jj]
                    mo_energy[1][i] += eri_mo[1][ii,jj] - eri_mo[1][ij,ij]
        assert(self.nmo < 210)
        ref_id = {'RHF': 0, 'ROHF': 1, 'UHF': 2}[self.ref]
        _psi4.prepare_chkpt(ref_id, mo, mo_energy, nelec, ms)
        _psi4.prepare_dpd_ints(self.ref, self.nmo, hcore_mo, eri_mo)

    def energy(self, key):
        '''Calculate correlation energy

        Parameters
        ----------
        key : {'CCSD'}
            CC solver
        '''
        assert ( key.upper() == 'CCSD' )
        return _psi4.energy(key.upper(), self.ref)

    def density(self, notation='D'):
        '''Calculate 1-particle and anti-symmetrized 2-particle density matrices.

        The density matrices correspond to the energy expression
            E = sum_pq Dpq hpq + 1/2 sum_pqrs Gpqrs <pq|rs>
        where hpq is NOT fock matrix
      
        Parameters
        ----------
        notation : {'D', 'M'}, optional
            Dirac notation or Mulliken notation.
            In Dirac notation, Gpqrs means p,r for electron 1 and q,s for
            electron 2.
            In Mulliken notation, Gpqrs means p,q for electron 1 and q,s for
            electron 2.

        Returns
        -------
        rdm1 : 2-D array or 2 2D arrays (alpha, beta)
            One-particle reduced density matrix.
        rdm2 : 4-D array or 4 4D arrays (aa, bb, ab, ba)
            two-particle reduced density matrix.
        '''
        _psi4.density()
        return self.read_dm(notation)

    def read_dm(self, notation='D'):
        '''Read 1-particle and anti-symmetrized 2-particle density matrices if
        density matrices have been calculated.
        * Hartree-Fock 1pdm and 2pdm are NOT included.
        * J,K of Fock operator has been included in 2pdm
        * Only get alpha spin rdm1 for RHF, should * 2 for spin free DM
      
        Parameters
        ----------
        notation : {'D', 'M'}, optional
            Dirac notation or Mulliken notation.
            In Dirac notation, Gpqrs means p,r for electron 1 and q,s for
            electron 2.
            In Mulliken notation, Gpqrs means p,q for electron 1 and q,s for
            electron 2.

        Returns
        -------
        rdm1 : 2-D array or 2 2D arrays (alpha, beta)
            One-particle reduced density matrix.
        rdm2 : 4-D array or 4 4D arrays (aa, bb, ab, ba)
            two-particle reduced density matrix.
        '''
        rdm1, rdm2 = _psi4.read_rdm(self.ref, self.nmo)
        if notation == 'M':
            if self.ref == 'RHF':
                return rdm1, _transform_D2M(rdm2)
            else:
                return rdm1, map(_transform_D2M, rdm2)
        else:
            return rdm1, rdm2

def _transform_D2M(rdm2_D):
    n = rdm2_D.shape[0]
    return rdm2_D.transpose(0,2,1,3)


if __name__ == '__main__':
    pass
