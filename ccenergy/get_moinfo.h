
namespace psi { namespace ccenergy {

template <class T>
void read_moinfo_from_chkpt(T& moinfo)
{
    moinfo.nirreps = chkpt_rd_nirreps();
    moinfo.nmo = chkpt_rd_nmo();
    moinfo.nso = chkpt_rd_nso();
    moinfo.nao = chkpt_rd_nao();
    for (int i=0; i<moinfo.nirreps; i++) {
        moinfo.labels[i] = (char *) malloc(sizeof(char)*5);
        ::memset(moinfo.labels[i], 0, sizeof(char)*5);
    }
    moinfo.enuc = chkpt_rd_enuc();
    moinfo.escf = chkpt_rd_escf();

    moinfo.sopi = chkpt_rd_sopi();
    moinfo.orbspi = chkpt_rd_orbspi();
    moinfo.openpi = chkpt_rd_openpi();
    moinfo.clsdpi = chkpt_rd_clsdpi();
    moinfo.frdocc = chkpt_rd_frzcpi();
    moinfo.fruocc = chkpt_rd_frzvpi();
}

}}
