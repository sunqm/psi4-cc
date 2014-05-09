set(SRC _psi4.c
  ${PSI4DIR}/src/bin/psi4/read_options.cc
  dpd_mo_ints.cc fake_function.cc psi_start.cc psi_stop.cc
  read_den_mat.cc psi4itrf.cc)

#add_library(psi4itrf SHARED ${SRC})
#add_library(psi4itrf ${SRC})
add_library(_psi4 SHARED ${SRC})
set_target_properties(_psi4 PROPERTIES PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/psi4)

target_link_libraries(_psi4
  ${ITRFLIB} ${PSILIB} ${ITRFLIB} ${Boost_LIBRARIES} ${LAPACK_LIBRARIES})

