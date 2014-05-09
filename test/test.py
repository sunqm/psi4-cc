#!/usr/bin/env python
import numpy
import psi4

n = 4; m = 2
c = numpy.linalg.svd(numpy.sin(numpy.arange(n*n).reshape(n,n)))[0]
c = numpy.array(c, order='C')
e = numpy.arange(n)*.1+.1
#fock = numpy.dot(c*e, c.T)
fock = numpy.diag(e)

nn = n*(n+1)/2
eri = numpy.cos(numpy.arange(nn**2)).reshape(nn,nn) * .01
eri = eri + eri.T

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
J = numpy.zeros_like(fock)
K = numpy.zeros_like(fock)
for i in range(m):
    J += eri_full[:,i,:,i]
    K += eri_full[:,i,i,:]
h1e = fock - (2*J-K)

ps = psi4.Solver()
with psi4.quite_run():
    ps.prepare('RHF', c, h1e, eri, m*2, 0)
    ecc = ps.energy('CCSD')
    rdm1, rdm2 = ps.density()

e1 = numpy.dot(rdm1.flatten(), h1e.flatten()) * 2
e2 = numpy.dot(rdm2.flatten(), eri_full.flatten()) * .5
assert(abs(-0.003480884520372-ecc) < 1e-10)
assert(abs( 0.003122970157749-e1) < 1e-10)
assert(abs(-0.006603854678121-e2) < 1e-10)
eijkl = (eri_full[:m,:m,:m,:m]*rdm2[:m,:m,:m,:m]).sum() * .5
eijka =((eri_full[:m,:m,:m,m:]*rdm2[:m,:m,:m,m:]).sum()
      + (eri_full[:m,:m,m:,:m]*rdm2[:m,:m,m:,:m]).sum()
      + (eri_full[:m,m:,:m,:m]*rdm2[:m,m:,:m,:m]).sum()
      + (eri_full[m:,:m,:m,:m]*rdm2[m:,:m,:m,:m]).sum()) * .5
eijab =((eri_full[:m,:m,m:,m:]*rdm2[:m,:m,m:,m:]).sum()
      + (eri_full[m:,m:,:m,:m]*rdm2[m:,m:,:m,:m]).sum()
      + (eri_full[m:,:m,:m,m:]*rdm2[m:,:m,:m,m:]).sum()
      + (eri_full[:m,m:,m:,:m]*rdm2[:m,m:,m:,:m]).sum()) * .5
eibja =((eri_full[:m,m:,:m,m:]*rdm2[:m,m:,:m,m:]).sum()
      + (eri_full[m:,:m,m:,:m]*rdm2[m:,:m,m:,:m]).sum()) * .5
eciab =((eri_full[m:,:m,m:,m:]*rdm2[m:,:m,m:,m:]).sum()
      + (eri_full[m:,m:,:m,m:]*rdm2[m:,m:,:m,m:]).sum()
      + (eri_full[m:,m:,m:,:m]*rdm2[m:,m:,m:,:m]).sum()
      + (eri_full[:m,m:,m:,m:]*rdm2[:m,m:,m:,m:]).sum()) * .5
eabcd = (eri_full[m:,m:,m:,m:]*rdm2[m:,m:,m:,m:]).sum() * .5
assert(abs( 0.000199604546228-eijkl) < 1e-10)
assert(abs(-0.000507942173007-eijka) < 1e-10)
assert(abs(-0.006798942593477-eijab) < 1e-10)
assert(abs( 0.000499123140399-eibja) < 1e-10)
assert(abs(-0.000074579201724-eciab) < 1e-10)
assert(abs( 0.000078881603460-eabcd) < 1e-10)
print 'restrict CCSD PASS'

with psi4.capture_stdout():
    ps.prepare('RHF', c, h1e, eri, m*2, 0)
    ecc = ps.energy('CCSD(T)')
    rdm1, rdm2 = ps.density()

e1 = numpy.dot(rdm1.flatten(), h1e.flatten()) * 2
e2 = numpy.dot(rdm2.flatten(), eri_full.flatten()) * .5
assert(abs(-0.0034808845203719995-ecc) < 1e-10)
assert(abs( 0.003122970157749-e1) < 1e-10)
assert(abs(-0.006603854678121-e2) < 1e-10)
eijkl = (eri_full[:m,:m,:m,:m]*rdm2[:m,:m,:m,:m]).sum() * .5
eijka =((eri_full[:m,:m,:m,m:]*rdm2[:m,:m,:m,m:]).sum()
      + (eri_full[:m,:m,m:,:m]*rdm2[:m,:m,m:,:m]).sum()
      + (eri_full[:m,m:,:m,:m]*rdm2[:m,m:,:m,:m]).sum()
      + (eri_full[m:,:m,:m,:m]*rdm2[m:,:m,:m,:m]).sum()) * .5
eijab =((eri_full[:m,:m,m:,m:]*rdm2[:m,:m,m:,m:]).sum()
      + (eri_full[m:,m:,:m,:m]*rdm2[m:,m:,:m,:m]).sum()
      + (eri_full[m:,:m,:m,m:]*rdm2[m:,:m,:m,m:]).sum()
      + (eri_full[:m,m:,m:,:m]*rdm2[:m,m:,m:,:m]).sum()) * .5
eibja =((eri_full[:m,m:,:m,m:]*rdm2[:m,m:,:m,m:]).sum()
      + (eri_full[m:,:m,m:,:m]*rdm2[m:,:m,m:,:m]).sum()) * .5
eciab =((eri_full[m:,:m,m:,m:]*rdm2[m:,:m,m:,m:]).sum()
      + (eri_full[m:,m:,:m,m:]*rdm2[m:,m:,:m,m:]).sum()
      + (eri_full[m:,m:,m:,:m]*rdm2[m:,m:,m:,:m]).sum()
      + (eri_full[:m,m:,m:,m:]*rdm2[:m,m:,m:,m:]).sum()) * .5
eabcd = (eri_full[m:,m:,m:,m:]*rdm2[m:,m:,m:,m:]).sum() * .5
assert(abs( 0.000199604546228-eijkl) < 1e-10)
assert(abs(-0.000507942173007-eijka) < 1e-10)
assert(abs(-0.006798942593477-eijab) < 1e-10)
assert(abs( 0.000499123140399-eibja) < 1e-10)
assert(abs(-0.000074579201724-eciab) < 1e-10)
assert(abs( 0.000078881603460-eabcd) < 1e-10)
print 'restrict CCSD(T) PASS'

########################################
e = fock.diagonal()
fock_a = fock
fock_b = numpy.diag(e*1.5)

eri_a = eri
eri_b = eri*.8
eri_ab = eri.copy()
eri_ab -= numpy.eye(eri_ab.shape[0])*.02
eri_fulla = numpy.empty((n,n,n,n))
eri_fullb = numpy.empty((n,n,n,n))
eri_fullab = numpy.empty((n,n,n,n))
eri_fullba = numpy.empty((n,n,n,n))
for i in range(n):
    for j in range(i+1):
        ij = i*(i+1)/2 + j
        for k in range(n):
            for l in range(k+1):
                kl = k*(k+1)/2 + l
                eri_fulla[i,k,j,l] = eri_fulla[j,k,i,l] = eri_fulla[i,l,j,k] \
                                   = eri_fulla[j,l,i,k] = eri_a[ij,kl]
                eri_fullb[i,k,j,l] = eri_fullb[j,k,i,l] = eri_fullb[i,l,j,k] \
                                   = eri_fullb[j,l,i,k] = eri_b[ij,kl]
                eri_fullab[i,k,j,l]= eri_fullab[j,k,i,l]= eri_fullab[i,l,j,k] \
                                   = eri_fullab[j,l,i,k]= eri_ab[ij,kl]
                eri_fullba[k,i,l,j]= eri_fullba[k,j,l,i]= eri_fullba[l,i,k,j] \
                                   = eri_fullba[l,j,k,i]= eri_ab[ij,kl]
Ja = numpy.zeros_like(fock)
Jb = numpy.zeros_like(fock)
Ja1= numpy.zeros_like(fock)
Jb1= numpy.zeros_like(fock)
Ka = numpy.zeros_like(fock)
Kb = numpy.zeros_like(fock)
for i in range(m):
    Ja += eri_fulla[:,i,:,i]
    Jb1+= eri_fullab[i,:,i,:]
    Ka += eri_fulla[:,i,i,:]
for i in range(m-1):
    Jb += eri_fullb[:,i,:,i]
    Kb += eri_fullb[:,i,i,:]
    Ja1+= eri_fullab[:,i,:,i]
h1e_a = fock_a - (Ja+Ja1-Ka)
h1e_b = fock_b - (Jb+Jb1-Kb)
with psi4.quite_run():
    ps.prepare('UHF', (c,c), (h1e_a,h1e_b), (eri_a,eri_b,eri_ab), m*2-1, 1)
    ecc = ps.energy('CCSD')
    rdm1, rdm2 = ps.density()

w = m - 1
e1 = numpy.dot(rdm1[0].flatten(), h1e_a.flatten()) \
   + numpy.dot(rdm1[1].flatten(), h1e_b.flatten())
e2 = numpy.dot(rdm2[0].flatten(), eri_fulla.flatten()) * .5 \
   + numpy.dot(rdm2[1].flatten(), eri_fullb.flatten()) * .5 \
   + numpy.dot(rdm2[2].flatten(), eri_fullab.flatten()) * .5 \
   + numpy.dot(rdm2[3].flatten(), eri_fullba.flatten()) * .5
assert(abs(-0.00566250260924-ecc) < 1e-10)
assert(abs( 0.00928927225526-e1) < 1e-10)
assert(abs(-0.01495177480547-e2) < 1e-10)
eIJKL = (eri_fulla [:m,:m,:m,:m]*rdm2[0][:m,:m,:m,:m]).sum() * .5
eijkl = (eri_fullb [:w,:w,:w,:w]*rdm2[1][:w,:w,:w,:w]).sum() * .5
eIjKl = (eri_fullab[:m,:w,:m,:w]*rdm2[2][:m,:w,:m,:w]).sum() * .5
eiJkL = (eri_fullba[:w,:m,:w,:m]*rdm2[3][:w,:m,:w,:m]).sum() * .5
eIJKA =((eri_fulla [:m,:m,:m,m:]*rdm2[0][:m,:m,:m,m:]).sum()
      + (eri_fulla [:m,:m,m:,:m]*rdm2[0][:m,:m,m:,:m]).sum()
      + (eri_fulla [:m,m:,:m,:m]*rdm2[0][:m,m:,:m,:m]).sum()
      + (eri_fulla [m:,:m,:m,:m]*rdm2[0][m:,:m,:m,:m]).sum()) * .5
eijka =((eri_fullb [:w,:w,:w,w:]*rdm2[1][:w,:w,:w,w:]).sum()
      + (eri_fullb [:w,:w,w:,:w]*rdm2[1][:w,:w,w:,:w]).sum()
      + (eri_fullb [:w,w:,:w,:w]*rdm2[1][:w,w:,:w,:w]).sum()
      + (eri_fullb [w:,:w,:w,:w]*rdm2[1][w:,:w,:w,:w]).sum()) * .5
eIjKa =((eri_fullab[:m,:w,:m,w:]*rdm2[2][:m,:w,:m,w:]).sum()
      + (eri_fullab[:m,:w,m:,:w]*rdm2[2][:m,:w,m:,:w]).sum()
      + (eri_fullab[:m,w:,:m,:w]*rdm2[2][:m,w:,:m,:w]).sum()
      + (eri_fullab[m:,:w,:m,:w]*rdm2[2][m:,:w,:m,:w]).sum()) * .5
eiJkA =((eri_fullba[:w,:m,:w,m:]*rdm2[3][:w,:m,:w,m:]).sum()
      + (eri_fullba[:w,:m,w:,:m]*rdm2[3][:w,:m,w:,:m]).sum()
      + (eri_fullba[:w,m:,:w,:m]*rdm2[3][:w,m:,:w,:m]).sum()
      + (eri_fullba[w:,:m,:w,:m]*rdm2[3][w:,:m,:w,:m]).sum()) * .5
assert(abs(-0.001961829412970-eIJKA) < 1e-10)
assert(abs(                   eijka) < 1e-10)
assert(abs(-0.000543696858059-eIjKa) < 1e-10)
assert(abs(-0.000543696858059-eiJkA) < 1e-10)
assert(abs( 0.000000561762988-eIjKl) < 1e-10)
assert(abs( 0.000000561762988-eiJkL) < 1e-10)
assert(abs( 0.000002349945814-eIJKL) < 1e-10)
assert(abs(                   eijkl) < 1e-10)
eIJAB =((eri_fulla [:m,:m,m:,m:]*rdm2[0][:m,:m,m:,m:]).sum()
      + (eri_fulla [m:,m:,:m,:m]*rdm2[0][m:,m:,:m,:m]).sum()
      + (eri_fulla [m:,:m,:m,m:]*rdm2[0][m:,:m,:m,m:]).sum()
      + (eri_fulla [:m,m:,m:,:m]*rdm2[0][:m,m:,m:,:m]).sum()) * .5
eijab =((eri_fullb [:w,:w,w:,w:]*rdm2[1][:w,:w,w:,w:]).sum()
      + (eri_fullb [w:,w:,:w,:w]*rdm2[1][w:,w:,:w,:w]).sum()
      + (eri_fullb [w:,:w,:w,w:]*rdm2[1][w:,:w,:w,w:]).sum()
      + (eri_fullb [:w,w:,w:,:w]*rdm2[1][:w,w:,w:,:w]).sum()) * .5
eIjAb =((eri_fullab[:m,:w,m:,w:]*rdm2[2][:m,:w,m:,w:]).sum()
      + (eri_fullab[m:,w:,:m,:w]*rdm2[2][m:,w:,:m,:w]).sum()
      + (eri_fullab[m:,:w,:m,w:]*rdm2[2][m:,:w,:m,w:]).sum()
      + (eri_fullab[:m,w:,m:,:w]*rdm2[2][:m,w:,m:,:w]).sum()) * .5
eiJaB =((eri_fullba[:w,:m,w:,m:]*rdm2[3][:w,:m,w:,m:]).sum()
      + (eri_fullba[w:,m:,:w,:m]*rdm2[3][w:,m:,:w,:m]).sum()
      + (eri_fullba[w:,:m,:w,m:]*rdm2[3][w:,:m,:w,m:]).sum()
      + (eri_fullba[:w,m:,w:,:m]*rdm2[3][:w,m:,w:,:m]).sum()) * .5
assert(abs(-0.00020898141735-eIJAB) < 1e-10)
assert(abs(                  eijab) < 1e-10)
assert(abs(-0.00558706749653-eIjAb) < 1e-10)
assert(abs(-0.00558706749653-eiJaB) < 1e-10)
eIBJA =((eri_fulla [:m,m:,:m,m:]*rdm2[0][:m,m:,:m,m:]).sum()
      + (eri_fulla [m:,:m,m:,:m]*rdm2[0][m:,:m,m:,:m]).sum()) * .5
eibja =((eri_fullb [:w,w:,:w,w:]*rdm2[1][:w,w:,:w,w:]).sum()
      + (eri_fullb [w:,:w,w:,:w]*rdm2[1][w:,:w,w:,:w]).sum()) * .5
eIbJa =((eri_fullab[:m,w:,:m,w:]*rdm2[2][:m,w:,:m,w:]).sum()
      + (eri_fullab[m:,:w,m:,:w]*rdm2[2][m:,:w,m:,:w]).sum()) * .5
eiBjA =((eri_fullba[:w,m:,:w,m:]*rdm2[3][:w,m:,:w,m:]).sum()
      + (eri_fullba[w:,:m,w:,:m]*rdm2[3][w:,:m,w:,:m]).sum()) * .5
assert(abs(-0.000017868968495-eIBJA) < 1e-10)
assert(abs(                   eibja) < 1e-10)
assert(abs(-0.000158582278839-eIbJa) < 1e-10)
assert(abs(-0.000158582278839-eiBjA) < 1e-10)
eCIAB =((eri_fulla [m:,:m,m:,m:]*rdm2[0][m:,:m,m:,m:]).sum()
      + (eri_fulla [m:,m:,:m,m:]*rdm2[0][m:,m:,:m,m:]).sum()
      + (eri_fulla [m:,m:,m:,:m]*rdm2[0][m:,m:,m:,:m]).sum()
      + (eri_fulla [:m,m:,m:,m:]*rdm2[0][:m,m:,m:,m:]).sum()) * .5
eciab =((eri_fullb [w:,:w,w:,w:]*rdm2[1][w:,:w,w:,w:]).sum()
      + (eri_fullb [w:,w:,:w,w:]*rdm2[1][w:,w:,:w,w:]).sum()
      + (eri_fullb [w:,w:,w:,:w]*rdm2[1][w:,w:,w:,:w]).sum()
      + (eri_fullb [:w,w:,w:,w:]*rdm2[1][:w,w:,w:,w:]).sum()) * .5
eCiAb =((eri_fullab[m:,:w,m:,w:]*rdm2[2][m:,:w,m:,w:]).sum()
      + (eri_fullab[m:,w:,:m,w:]*rdm2[2][m:,w:,:m,w:]).sum()
      + (eri_fullab[m:,w:,m:,:w]*rdm2[2][m:,w:,m:,:w]).sum()
      + (eri_fullab[:m,w:,m:,w:]*rdm2[2][:m,w:,m:,w:]).sum()) * .5
ecIaB =((eri_fullba[w:,:m,w:,m:]*rdm2[3][w:,:m,w:,m:]).sum()
      + (eri_fullba[w:,m:,:w,m:]*rdm2[3][w:,m:,:w,m:]).sum()
      + (eri_fullba[w:,m:,w:,:m]*rdm2[3][w:,m:,w:,:m]).sum()
      + (eri_fullba[:w,m:,w:,m:]*rdm2[3][:w,m:,w:,m:]).sum()) * .5
assert(abs( 0.00000014066422328-eCIAB) < 1e-10)
assert(abs(                     eciab) < 1e-10)
assert(abs(-0.00000844074284043-eCiAb) < 1e-10)
assert(abs(-0.00000844074284043-ecIaB) < 1e-10)
eABCD = (eri_fulla [m:,m:,m:,m:]*rdm2[0][m:,m:,m:,m:]).sum() * .5
eabcd = (eri_fullb [w:,w:,w:,w:]*rdm2[1][w:,w:,w:,w:]).sum() * .5
eAbCd = (eri_fullab[m:,w:,m:,w:]*rdm2[2][m:,w:,m:,w:]).sum() * .5
eaBcD = (eri_fullba[w:,m:,w:,m:]*rdm2[3][w:,m:,w:,m:]).sum() * .5
assert(abs(-0.000000926223721-eABCD) < 1e-10)
assert(abs(                  -eabcd) < 1e-10)
assert(abs(-0.000085104083204-eAbCd) < 1e-10)
assert(abs(-0.000085104083204-eaBcD) < 1e-10)

print 'unrestrict CCSD PASS'
