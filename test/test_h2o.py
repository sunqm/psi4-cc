#!/usr/bin/env python

import numpy

import psi4

import sys
sys.path.append('/home/seba') #Folder in which PySCF is installed
from pyscf import gto, scf, symm, ao2mo

mol = gto.Mole()
mol.verbose = 5
mol.output = 'out_h2o'
mol.atom = [
    [8 , (0. , 0.     , 0.)],
    [1 , (0. , -0.757 , 0.587)],
    [1 , (0. , 0.757  , 0.587)]]

mol.basis = {'H': 'cc-pvdz',
             'O': 'cc-pvdz',}
mol.build()
rhf = scf.RHF(mol)
print rhf.scf()

L = rhf.mo_coeff.shape[1]

hcore_mo = reduce(numpy.dot, (rhf.mo_coeff.T, rhf.get_hcore(), rhf.mo_coeff))
eri = ao2mo.outcore.full_iofree(mol, rhf.mo_coeff, compact=False).reshape(L,L,L,L)
ps = psi4.Solver()
with psi4.capture_stdout():
    ps.prepare('RHF', rhf.mo_coeff, hcore_mo, ao2mo.restore(4,eri,L), mol.nelectron)
    ecc = ps.energy('CCSD')
    rdm1, rdm2 = ps.density() #RDM2 in physics notation!

print rdm1.shape
print rdm2.shape

e1 = numpy.einsum('ij,ij->', rdm1, hcore_mo) * 2
e2 = numpy.einsum('ijkl,ikjl->', rdm2, eri) * 0.5

# ecc should be -0.213343234276
print ecc, e1, e2, e1+e2
