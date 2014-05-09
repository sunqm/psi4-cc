'''
ps = psi4.Solver
with psi4.quite_run():
    ps.prepare_chkpt(mo_coeff, fock_on_mo, nelec, e_scf, nuclear_repulsion)
    ecc = ps.energy('CCSD', c.shape[1], hcore_on_mo, eri_on_mo)
    rdm1, rdm2 = ps.density(mo_coeff.shape[1])

    eccsdt = ps.energy('CCSD(T)', c.shape[1], hcore_on_mo, eri_on_mo)
    rdm1, rdm2 = ps.density(mo_coeff.shape[1])
'''

from wrapper import *

__all__ = filter(lambda s: not s.startswith('_'), dir())
