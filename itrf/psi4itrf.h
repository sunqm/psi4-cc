#define CACHELEVEL 4



#if defined __cplusplus
extern "C" {
#endif
void psi4itrf_init_env(const char *outfile, const char *tmpdir,
                       int argc4MPI, char **argv4MPI);
void psi4itrf_del_env();
void psio_clean();
double psi4ccsd_energy(char *ref_wfn);
double psi4ccsd_t(char *ref_wfn);
void psi4cc_density();

void dpd_mo_ints_rhf(const double *h1e, const double *eri, int nmo);
void dpd_mo_ints_uhf(double *hcore_a, double *hcore_b,
                     double *eri_aa, double *eri_bb, double *eri_ab,
                     int nmo);
void read_cc_dm_rhf(double *rdm1, double *rdm2, int nmo);
void read_cc_dm_uhf(double *rdm1_a, double *rdm1_b,
                    double *rdm2_aa, double *rdm2_bb,
                    double *rdm2_ab, double *rdm2_ba, int nmo);
#if defined __cplusplus
}
#endif
