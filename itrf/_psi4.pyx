import os
from libc.stdlib cimport malloc, free
import numpy
cimport numpy
cimport cython
#cython: boundscheck=False
#cython: wraparound=False

cdef extern from "Python.h":
    char *PyString_AsString(object)

cdef extern from "libchkpt/chkpt.h":
    # moinfo in ccenergy
    int chkpt_init(int status)
    int chkpt_close()
    void chkpt_wt_nirreps(int)
    void chkpt_wt_irr_labs(char **)
    void chkpt_wt_nmo(int)
    void chkpt_wt_nso(int, char *)
    void chkpt_wt_nao(int, char *)
    void chkpt_wt_iopen(int)
    void chkpt_wt_ref(int)
    void chkpt_wt_etot(double)
    void chkpt_wt_escf(double)
    void chkpt_wt_eref(double)
    void chkpt_wt_enuc(double)
    void chkpt_wt_orbspi(int *)
    void chkpt_wt_clsdpi(int *)
    void chkpt_wt_openpi(int *)
    void chkpt_wt_phase_check(int)
    void chkpt_wt_sopi(int *, char *)
    void chkpt_wt_nfzc(int)
    void chkpt_wt_nfzv(int)
    void chkpt_wt_frzcpi(int *)
    void chkpt_wt_frzvpi(int *)
    void chkpt_wt_evals(double *)
    void chkpt_wt_scf(double **)
    void chkpt_wt_alpha_evals(double *)
    void chkpt_wt_beta_evals(double *)
    void chkpt_wt_alpha_scf(double **)
    void chkpt_wt_beta_scf(double **)
    void chkpt_wt_enuc(double)

    # moinfo in transqt2
    void chkpt_wt_efzc(double) # energy of frozen core

cdef extern from "psi4itrf.h":
    void psi4itrf_init_env(char *outfile, char *tmpdir, 
                           unsigned long max_memory,
                           int argc4MPI, char **argv4MPI)
    void psi4itrf_del_env()
    void psio_clean()
    double psi4ccsd_energy(char *ref_wfn)
    double psi4ccsd_t(char *ref_wfn)
    void psi4cc_density()
    void dpd_mo_ints_rhf(double *h1e, double *eri, int nmo)
    void dpd_mo_ints_uhf(double *hcore_a, double *hcore_b,
                         double *eri_aa, double *eri_bb, double *eri_ab,
                         int nmo)
    void read_cc_dm_rhf(double *rdm1, double *rdm2, int nmo)
    void read_cc_dm_uhf(double *rdm1_a, double *rdm1_b,
                        double *rdm2_aa, double *rdm2_bb,
                        double *rdm2_ab, double *rdm2_ba, int nmo)
 

def set_psi4_env(tmpdir, outfile='stdout', argv=[], max_memory=(1<<30)):
    '''
    set_psi4_env()
    
    setup psi4 running environment
    '''
    cdef int argc = len(argv)
    cdef char **argv_buf = <char **>malloc(argc * sizeof(char *))
    for i, a in enumerate(argv):
        argv_buf[i] = PyString_AsString(a)
    psi4itrf_init_env(outfile, tmpdir, max_memory, argc, argv_buf)
    free(argv_buf)

def del_psi4_env():
    '''
    del_psi4_env()

    clean up psi4 running trashes
    '''
    psi4itrf_del_env()
    try:
        os.remove('psi.timer.dat')
        os.remove('psi.ijk.dat')
    except:
        pass

# NOTE: It is really importatnt to call me before prepare anything.  When the
# number of basis changes, psio will raise errors because of the confliction
# of the temporary file size.
def psioclean_before_preparing():
    psio_clean()

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def prepare_chkpt(ref_wfn_id, mo, mo_energy, nelec, ms):
    chkpt_init(0) # PSIO_OPEN_NEW

    if ref_wfn_id == 0: # RHF
        nao, nmo = mo.shape
    else:
        nao, nmo = mo[0].shape
    chkpt_wt_nmo(nmo)
    chkpt_wt_nso(nao, '')
    chkpt_wt_nao(nao, '')

    cdef numpy.ndarray mo_a, mo_b
    cdef numpy.ndarray mo_e_a, mo_e_b
    cdef double **mo_buf = <double **>malloc(sizeof(double *))
    if ref_wfn_id == 0: # RHF
        mo_a = mo
        mo_e_a = mo_energy
        mo_buf[0] = <double *>mo_a.data
        chkpt_wt_scf(mo_buf) # mo coeffs
        chkpt_wt_evals(<double *>mo_e_a.data)
    else:
        mo_a, mo_b = mo
        mo_e_a, mo_e_b = mo_energy
        mo_buf[0] = <double *>mo_a.data
        chkpt_wt_alpha_scf(mo_buf)
        mo_buf[0] = <double *>mo_b.data
        chkpt_wt_beta_scf(mo_buf)
        chkpt_wt_alpha_evals(<double *>mo_e_a.data)
        chkpt_wt_beta_evals(<double *>mo_e_b.data)
    free(mo_buf)

    chkpt_wt_ref(ref_wfn_id)
    chkpt_wt_nirreps(1)

    cdef char **labels = <char **>malloc(1 * sizeof(char *))
    labels[0] = PyString_AsString('A')
    chkpt_wt_irr_labs(labels)
    free(labels)

    chkpt_wt_etot(0)
    chkpt_wt_escf(0)
    chkpt_wt_eref(0)
    chkpt_wt_enuc(0)
    #chkpt_wt_efzc(e_fzc)

    chkpt_wt_orbspi([nmo]) # num. mo per irrep
    chkpt_wt_clsdpi([(nelec-ms)/2]) # num. doubly occupied mo per irrep
    chkpt_wt_openpi([ms]) # num. singly occupied mo per irrep
    chkpt_wt_iopen(ms) # libchkpt/iopen.cc
    chkpt_wt_phase_check(0) # =0 forces psi4 param.restart = 0
    chkpt_wt_sopi([nao], '') # num. symmetrized orbitals per irrep
    chkpt_wt_frzcpi([0]) # num. frozen cores per irrep
    chkpt_wt_frzvpi([0]) # num. frozen virtuals per irrep
    chkpt_close()
    return 0

def prepare_dpd_ints(ref_wfn, nmo, hcore_mo, eri_mo):
    cdef numpy.ndarray hcore_a, hcore_b
    cdef numpy.ndarray eri_aa, eri_bb, eri_ab
    if ref_wfn == 'UHF':
        hcore_a, hcore_b = hcore_mo
        eri_aa, eri_bb, eri_ab = eri_mo
        dpd_mo_ints_uhf(<double *>hcore_a.data, <double *>hcore_b.data,
                        <double *>eri_aa.data, <double *>eri_bb.data,
                        <double *>eri_ab.data, nmo)
    else:
        hcore_a = hcore_mo
        eri_aa = eri_mo
        dpd_mo_ints_rhf(<double *>hcore_a.data, <double *>eri_aa.data, nmo)

# return 2-rdm in Dirac notation
def read_rdm(ref_wfn, nmo):
    cdef numpy.ndarray rdm1a, rdm1b
    cdef numpy.ndarray rdm2aa, rdm2bb, rdm2ab, rdm2ba
    if ref_wfn == 'UHF':
        rdm1a = numpy.zeros((nmo,nmo), dtype=numpy.double)
        rdm1b = numpy.zeros((nmo,nmo), dtype=numpy.double)
        rdm2aa = numpy.zeros((nmo,nmo,nmo,nmo), dtype=numpy.double)
        rdm2bb = numpy.zeros((nmo,nmo,nmo,nmo), dtype=numpy.double)
        rdm2ab = numpy.zeros((nmo,nmo,nmo,nmo), dtype=numpy.double)
        rdm2ba = numpy.zeros((nmo,nmo,nmo,nmo), dtype=numpy.double)
        read_cc_dm_uhf(<double *>rdm1a.data, <double *>rdm1b.data,
                       <double *>rdm2aa.data, <double *>rdm2bb.data,
                       <double *>rdm2ab.data, <double *>rdm2ba.data, nmo)
        return (rdm1a,rdm1b), (rdm2aa,rdm2bb,rdm2ab,rdm2ba)
    else:
        rdm1a = numpy.zeros((nmo,nmo), dtype=numpy.double)
        rdm2aa = numpy.zeros((nmo,nmo,nmo,nmo), dtype=numpy.double)
        read_cc_dm_rhf(<double *>rdm1a.data, <double *>rdm2aa.data, nmo)
        return rdm1a, rdm2aa

def energy(key, ref_wfn):
    cdef char *refname = PyString_AsString(ref_wfn.upper())
    if key == 'CCSD(T)':
        return psi4ccsd_t(refname);
    else: # CCSD
        return psi4ccsd_energy(refname);

def density():
    psi4cc_density()
