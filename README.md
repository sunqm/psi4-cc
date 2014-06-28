psi4-cc
=======

python wrapper for Psi4 Coupled Cluster module

Installation
------------

* Compile Psi4
    - Download Psi4 from

        https://github.com/psi4/psi4public

    - Change Psi4 CMakeLists.txt. Before the line "find_package(Boost ...",
      add BOOST path as follows

        set(Boost_INCLUDE_DIR /path/to/boost/include)
        set(Boost_LIBRARY_DIR /path/to/boost/lib)

    - configure and make

        mkdir obj && cd obj
        ../configure.cmake \
        --with-lapack="-L$MKLROOT/lib/intel64 -I$MKLROOT/include
        -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lmkl_avx" \
        --with-ldflags="-lrt -lm" --with-cxxflags="-fPIC" \
        --with-f77flags="-fPIC" \
        --with-cxx=/usr/bin/g++ --with-plugin
        make

* Compile psi4-cc
    - Change CMakeLists.txt at the beginning

        set(PSI4DIR /path/to/psi4/source)
        set(PSI4BUILDIR ${PSI4DIR}/obj)

    - configure and make
        
        cmake28 . && make

    - tests

        cd test && python test.py

* Install psi4-cc to the place where python can find

    cd /path/to/psi4-cc/source; cp -a psi4 /path/that/is/included/in/$PYTHONPATH

