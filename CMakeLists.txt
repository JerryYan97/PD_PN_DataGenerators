cmake_minimum_required(VERSION 3.17)

project(pd_pn_generators)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# set(CMAKE_VERBOSE_MAKEFILE ON)

# Check environment variables
if(NOT DEFINED ENV{LIBIGL_DIR} OR NOT DEFINED ENV{MKL_ROOT} OR NOT DEFINED ENV{TBB_ROOT})
    message(FATAL_ERROR "Your environment variables are not set.")
endif()

include_directories($ENV{MKL_ROOT}/include)
include_directories($ENV{TBB_ROOT}/include)

# Find Intel OneAPI libraries
find_library(MKL_LIB NAMES mkl_intel_lp64 HINTS $ENV{MKL_ROOT}/lib/intel64 REQUIRED)
find_library(MKL_TBB NAMES mkl_tbb_thread HINTS $ENV{MKL_ROOT}/lib/intel64 REQUIRED)
find_library(MKL_CORE NAMES mkl_core HINTS $ENV{MKL_ROOT}/lib/intel64 REQUIRED)
find_library(TBB_LIB NAMES tbb HINTS $ENV{TBB_ROOT}/lib/intel64/gcc4.8 REQUIRED)

# Add tbb library to rpath
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags,-rpath,$ENV{TBB_ROOT}/lib/intel64/gcc4.8")

# libigl settings
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)


FIND_PACKAGE(Threads REQUIRED)
FIND_PACKAGE(LIBIGL REQUIRED)

add_executable(
        pd_pn_generators
        Src/main.cpp
        Src/macro.h
        Src/app.cpp
        Src/app.h
        Src/Simulators/DataGenerator.cpp
        Src/Simulators/DataGenerator.h
        Src/Simulators/PNSimulator.cpp
        Src/Simulators/PNSimulator.h
        Src/Simulators/PNTrainingDataGenerator.cpp
        Src/Simulators/PNTrainingDataGenerator.h
        Src/common/utils/AppReader.cpp
        Src/common/utils/AppReader.h
        Src/common/utils/AppWriter.cpp
        Src/common/utils/AppWriter.h
        Src/common/utils/InfoStructs.h Src/common/utils/Triplet.h Src/Simulators/ForceFields/ForceField.h Src/Simulators/ForceFields/DirectForceField.h)

target_link_libraries(pd_pn_generators
                      PUBLIC ${TBB_LIB}
                      PUBLIC ${MKL_LIB}
                      PUBLIC ${MKL_TBB}
                      PUBLIC ${MKL_CORE}
                      PUBLIC ${CMAKE_THREAD_LIBS_INIT}
                      PUBLIC igl::core
                      PUBLIC igl::opengl_glfw
                      )
