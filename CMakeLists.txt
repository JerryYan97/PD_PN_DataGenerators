cmake_minimum_required(VERSION 3.17)

project(pd_pn_generators)

set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(ThirdParties/eigen-3.3.9/)
include_directories(/opt/intel/oneapi/mkl/2021.1.1/include)

# Find Intel OneAPI libraries
find_library(MKL_LIB NAMES libmkl_intel_lp64 PATHS /opt/intel/oneapi/mkl/2021.1.1/lib/intel64)
if(NOT MKL_LIB)
    message(FATAL_ERROR "mkl library not found")
endif()

find_library(MKL_TBB NAMES libmkl_tbb_thread PATHS /opt/intel/oneapi/mkl/2021.1.1/lib/intel64)
if(NOT MKL_TBB)
    message(FATAL_ERROR "MKL_TBB library not found")
endif()
message(${MKL_TBB})

find_library(MKL_CORE NAMES libmkl_core PATHS /opt/intel/oneapi/mkl/2021.1.1/lib/intel64)
if(NOT MKL_CORE)
    message(FATAL_ERROR "mkl CORE library not found")
endif()

find_library(TBB_ONEAPI NAMES libtbb.so.12 PATHS /opt/intel/oneapi/tbb/2021.1.1/lib/intel64/gcc4.8)
if(NOT TBB_ONEAPI)
    message(FATAL_ERROR "TBB library not found")
endif()
message(${TBB_ONEAPI})

# Add tbb library to rpath
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags,-rpath,/opt/intel/oneapi/tbb/2021.1.1/lib/intel64/gcc4.8")

FIND_PACKAGE(Threads REQUIRED)

add_executable(
        pd_pn_generators
        Src/main.cpp
        Src/Simulators/PNSimulator.cpp
        Src/Simulators/PNSimulator.h Src/macro.h Src/app.cpp Src/app.h)

target_link_libraries(pd_pn_generators
                      PUBLIC ${TBB_ONEAPI}
                      PUBLIC ${MKL_LIB}
                      PUBLIC ${MKL_TBB}
                      PUBLIC ${MKL_CORE}
                      PUBLIC ${CMAKE_THREAD_LIBS_INIT}
                      )
