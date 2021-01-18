cmake_minimum_required(VERSION 3.17)
project(pd_pn_generators)

set(CMAKE_CXX_STANDARD 14)

include_directories(ThirdParties/eigen-3.3.9/)

add_executable(
        pd_pn_generators
        Src/main.cpp
        Src/Simulators/PNSimulator.cpp
        Src/Simulators/PNSimulator.h Src/macro.h)

