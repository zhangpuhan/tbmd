cmake_minimum_required(VERSION 3.13)
project(project)

set(CMAKE_CXX_STANDARD 11)
set(Torch_DIR /opt/libtorch/share/cmake/Torch)

add_executable(project main.cpp analysis.cpp analysis.hpp potential.cpp
        potential.hpp tbmd.cpp tbmd.hpp test.cpp test.hpp units.hpp util.cpp util.hpp vec3.hpp)

find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Torch REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(project ${LAPACK_LIBRARIES} ${BLAS_LIBARIES} ${ARMADILLO_LIBRARIES} ${TORCH_LIBRARIES})