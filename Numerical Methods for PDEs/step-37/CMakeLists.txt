cmake_minimum_required(VERSION 3.16)

project(step-37 CXX)

include(../common/cmake-common.cmake)

# In Apptainer, deal.II is already installed in the container environment.
# CMake will find it automatically. If not, you may set DEAL_II_DIR manually.
find_package(deal.II REQUIRED HINTS ${DEAL_II_DIR})

deal_ii_initialize_cached_variables()

add_executable(step-37 step-37.cc)

deal_ii_setup_target(step-37)
