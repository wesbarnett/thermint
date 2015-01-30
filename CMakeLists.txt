cmake_minimum_required (VERSION 2.8)
project (ThermodynamicIntegration)

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=gnu++0x)
endif()

find_package(Boost)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
else(Boost_FOUND)
    message(FATAL_ERROR "Could not find Boost.")
endif()

set(GROMACS_ROOT $ENV{GMXPREFIX})
if (NOT GROMACS_ROOT)
    message(FATAL_ERROR "Could not find GROMACS directory. Ensure that GMXPREFIX environmental variable is set")
endif()
message("GROMACS installation directory = ${GROMACS_ROOT}")
find_library(GROMACS_LIBRARY NAMES libgromacs PATH ${GROMACS_ROOT}/lib64)
include_directories(${GROMACS_ROOT}/include)

add_executable(thermint ti.cpp ${GROMACS_ROOT}/include/gromacs/fileio/enxio.h)
target_link_libraries(thermint gromacs)
target_link_libraries(thermint boost_program_options)

install (TARGETS thermint DESTINATION bin)
install (FILES ti.h DESTINATION include)
