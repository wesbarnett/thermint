cmake_minimum_required (VERSION 2.8)
project (ThermodynamicIntegration)

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=gnu++0x)
endif()

find_package(Boost)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
endif()

find_path(GROMACS_GROMACS_INCLUDE gromacs/fileio/enxio.h)
set(GROMACS_INCLUDE ${GROMACS_GROMACS_INCLUDE}/include)
find_library(GROMACS_LIBRARY NAMES libgromacs PATH ${GROMACS_ROOT}/lib64)
include_directories(${GROMACS_INCLUDE})
message("GROMACS include dir = ${GROMACS_INCLUDE}")

add_executable(thermint ti.cpp ${GROMACS_INCLUDE}/gromacs/fileio/enxio.h)
target_link_libraries(thermint gromacs)
target_link_libraries(thermint boost_program_options)

install (TARGETS thermint DESTINATION bin)
install (FILES ti.h DESTINATION include)
