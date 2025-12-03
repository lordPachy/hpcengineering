cmake_minimum_required(VERSION 3.12.0)
project(lab06 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(exercise-01 src/exercise-01.cpp src/LinearElasticity.cpp)
deal_ii_setup_target(exercise-01)
