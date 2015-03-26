psi4cc: a python wrapper for psi4's coupled cluster module
==========================================================
Written by Qiming Sun and Sebastian Wouters

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

Information
-----------

This version of psi4cc is adapted from Qiming Sun's version:
[https://github.com/sunqm/psi4-cc](https://github.com/sunqm/psi4-cc).
Qiming's version works for [psi4public](https://github.com/psi4/psi4public)
commit a6a30638bb19448528ff016f6ff37960560d8402 (April 6, 2014).

This repository is intended to keep psi4cc working for later versions as well.
I have tested it on commit 2c73a1f3751103fe2b7abc2b756d11f2d6556207
(March 24, 2015).

Installation
------------

### 1. Get psi4 (commit 2c73a1f3751103fe2b7abc2b756d11f2d6556207 or later)

    > git clone 'https://github.com/psi4/psi4public.git'
    > cd psi4public
    > mkdir objects
    > cd objects
    > BLA_VENDOR=Intel10_64lp CXX=icpc CC=icc CXXFLAGS="-fPIC" F77FLAGS="-fPIC" cmake .. -DENABLE_PLUGINS=ON -DENABLE_DUMMY_PLUGIN=ON
    > BLA_VENDOR=Intel10_64lp CXX=icpc CC=icc make
    > make install

### 2. Compile psi4cc

In order to compile psi4cc, please change CMakeLists.txt at the beginning:

    set(PSI4DIR /path/to/psi4/source)
    set(PSI4BUILDIR ${PSI4DIR}/obj)

Then configure and compile:

    > BLA_VENDOR=Intel10_64lp CXX=icpc CC=icc cmake .
    > BLA_VENDOR=Intel10_64lp CXX=icpc CC=icc make

### 3. Test psi4cc

[test.py](test/test.py) relies only on psi4 and psi4cc.
[test_h2o.py](test/test_h2o.py) requires
[pyscf](https://github.com/sunqm/pyscf).

    > cd tests
    > python test.py
    > python test_h2o.py

