import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

WORKDIR = os.path.realpath('.')
PSIBUILDIR = '/home/sunqm/Software/src/QC/psi4/obj'

# lack Boost
PSILIB = '''mints_wrapper dfmp2 dfocc scf psimrcc ccenergy ccsort ccdensity
  transqt2 cctriples scf_solver fock dcft lmp2 mcscf sapt dftsapt
  sapt_solver cchbar cclambda transqt ccresponse detci occ
  mrcc fnocc cceom adc thermo functional disp thce 3index deriv_wrapper
  optking findif mints trans dpd chkpt iwl psio qt ciomr options
  moinfo util stable scfgrad util diis plugin parallel'''.split()
ITRFLIB = '''psi4itrf psi4itrfccdensity psi4itrfccenergy'''.split()
BLASLAPCK = ['blas','lapack']
print PSILIB+ITRFLIB

ext_modules = [Extension('psi4',
                         ['psi4.pyx'],
                         language='c++',
                         #libraries=['psi4itrf'],
                         libraries=ITRFLIB+PSILIB+ITRFLIB+BLASLAPCK,
                         library_dirs = [WORKDIR+'/../lib',
                                         PSIBUILDIR+'/lib'],
                         #define_macros = [('DEBUG_ON',1)],
                         include_dirs=[numpy.get_include(),
                                       '/home/sunqm/Software/opt/psi4/include',
                                       '/home/sunqm/Software/opt/psi4/src/bin',
                                       '/home/sunqm/Software/opt/psi4/src/lib'],
                        )]

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    script_args = ['build_ext', '--inplace']
)

# build .so with    python setup.py build_ext -i
