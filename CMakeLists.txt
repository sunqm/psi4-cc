cmake_minimum_required (VERSION 2.8)
project (psi4cc)

set(PSI4DIR $ENV{HOME}/psi4public)
set(PSI4BUILDIR ${PSI4DIR}/objects)
message("psi4 src dir is  ${PSI4DIR}\n")

#get_filename_component(ROOT_SRC_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
#set(CMAKE_MODULE_PATH "${ROOT_SRC_DIR}/cmake/")

set(LAPACK_LIBRARIES "-L${MKL_LIB_PATH} -I${MKLINCLUDE} -lmkl")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LIBRARIES}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -fPIC")

#set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(
  ${PROJECT_BINARY_DIR}/include 
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}
  ${PSI4DIR}/include
  ${PSI4DIR}/src/lib
  ${PSI4DIR}/src/bin)

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(PythonLibs REQUIRED)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.get_include()"
  OUTPUT_VARIABLE NUMPY_INCLUDE)
include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE})

set(Boost_USE_STATIC_LIBS off)
set(Boost_USE_MULTITHREADED off)
set(Boost_COMPONENTS filesystem python regex serialization system thread)
find_package(Boost 1.47.0 COMPONENTS ${Boost_COMPONENTS})
#target_link_libraries(${Boost_LIBRARY_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

configure_file(
  "${PROJECT_SOURCE_DIR}/include/psiconfig.h.in"
  "${PROJECT_BINARY_DIR}/include/psiconfig.h")
#
# The location of compiled libraries and executables
#
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(PSILIB mints_wrapper dfmp2 dfocc scf psimrcc ccenergy ccsort ccdensity
  transqt2 cctriples scf_solver fock dcft lmp2 mcscf sapt dftsapt
  sapt_solver cchbar cclambda transqt transqt2 ccresponse detci occ
  mrcc fnocc cceom adc thermo functional disp thce 3index deriv_wrapper
  optking findif mints trans dpd chkpt iwl psio qt ciomr options
  moinfo util stable scfgrad util diis plugin parallel parallel2
  psi4util PsiUtil efp efp_solver interface_libefp)
set(ITRFLIB psi4itrfccsort psi4itrfccdensity psi4itrfccenergy psi4itrfcclambda psi4itrfcchbar psi4cctriples psi4qt)
link_directories(
  ${PSI4BUILDIR}/lib
  ${PROJECT_BINARY_DIR}/lib)

add_subdirectory(ccenergy)
add_subdirectory(cchbar)
add_subdirectory(cclambda)
add_subdirectory(ccsort)
add_subdirectory(ccdensity)
add_subdirectory(cctriples)
add_subdirectory(itrf)
add_subdirectory(qt)

