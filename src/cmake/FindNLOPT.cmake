#stuff to find and download/check out nlopt or find one installed on the system.


#include(ExternalProject)
#ExternalProject_Add(nlopt_project
#  GIT_REPOSITORY https://github.com/stevengj/nlopt.git
#  STEP_TARGETS   build
#  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/nlopt
#  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#)
#ExternalProject_Get_Property(nlopt_project install_dir)
#include_directories(${install_dir}/include)

#message("NLOPT LIBS ")

#include(ExternalProject)
#ExternalProject_Add(foobar
#  GIT_REPOSITORY https://github.com/stevengj/nlopt.git
#  GIT_TAG        origin/release/1.2.3
#  STEP_TARGETS   build
#)

#TODO: detect from system, automatically download, etc. etc.
#TODO: Make sure to statically link

option(NLOPT_CXX ON)
option (NLOPT_FORTRAN OFF)
option (BUILD_SHARED_LIBS OFF) #We want to force static linkage
option (NLOPT_PYTHON OFF)
option (NLOPT_OCTAVE OFF)
option (NLOPT_MATLAB OFF)
option (NLOPT_GUILE OFF)
option (NLOPT_SWIG OFF)
option (NLOPT_TESTS OFF)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/nlopt ${CMAKE_BINARY_DIR}/nlopt)
#what about the .hpp?
