#!/usr/bin/env python

import numpy

import psi4
import gto
import scf
import ao2mo

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

hcore_mo = reduce(numpy.dot, (rhf.mo_coeff.T, rhf.get_hcore(), rhf.mo_coeff))
eri = ao2mo.gen_int2e_ao2mo(mol, rhf.mo_coeff)
ps = psi4.Solver()
with psi4.capture_stdout():
    ps.prepare('RHF', rhf.mo_coeff, hcore_mo, eri, mol.nelectron)
    ecc = ps.energy('CCSD')
    rdm1, rdm2 = ps.density()

n = rhf.mo_coeff.shape[1]
eri_full = numpy.empty((n,n,n,n))
for i in range(n):
    for j in range(i+1):
        ij = i*(i+1)/2 + j
        for k in range(n):
            for l in range(k+1):
                kl = k*(k+1)/2 + l
                eri_full[i,k,j,l] = eri[ij,kl]
                eri_full[j,k,i,l] = eri[ij,kl]
                eri_full[i,l,j,k] = eri[ij,kl]
                eri_full[j,l,i,k] = eri[ij,kl]
e1 = numpy.dot(rdm1.flatten(), hcore_mo.flatten()) * 2
e2 = numpy.dot(rdm2.flatten(), eri_full.flatten()) * .5

# ecc should be -0.213343234276
print ecc, e1, e2, e1+e2
