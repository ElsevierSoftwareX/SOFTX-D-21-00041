###
# @file   CMakeLists.txt
#
# @author David S. Kammer <dkammer@ethz.ch>
# @author Gabriele Albertini <ga288@cornell.edu>
# @author Chun-Yu Ke <ck659@cornell.edu>
#
# @date creation: Fri Feb 5 2021
# @date last modification: Fri Feb 5 2021
#
# @brief  TODO
#
#
# Copyright (C) 2021 ETH Zurich (David S. Kammer)
#
# This file is part of uguca.
#
# uguca is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# uguca is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with uguca.  If not, see <https://www.gnu.org/licenses/>.
####
cmake_minimum_required(VERSION 3.1)
project(uguca)
enable_language(CXX)

set(UGUCA_MAJOR_VERSION 0)
set(UGUCA_MINOR_VERSION 9)
set(UGUCA_PATCH_VERSION 0)
set(UGUCA_COPYRING "2021 ETH Zurich (David S. Kammer)")
set(UGUCA_AUTHOR "David S. Kammer")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -std=c++11 " CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG " CACHE STRING "" FORCE )
set(CMAKE_MACOSX_RPATH TRUE)

SET(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
  "Build options: None Debug Release RelWithDebInfo MinSizeRel."
  FORCE )

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")
mark_as_advanced(
  BUILD_SHARED_LIBS)

find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)
find_package(Git)

# ------------------------------------------------------
# for KernelCollection to find kernels
set(KERNELS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/kernels/")

# config file
configure_file(uca_config.hh.in
  "${CMAKE_CURRENT_BINARY_DIR}/src/uca_config.hh" @ONLY)

#===============================================================================
# PERFORMANCE FEATURES
#===============================================================================

#-------------------------------------------------------
# BLAS
option(UCA_USE_BLAS "Use BLAS dot product for kernel convolution" OFF)
if (UCA_USE_BLAS)
  find_package(BLAS REQUIRED)
endif (UCA_USE_BLAS)

#-------------------------------------------------------
# OPENMP
option(UCA_USE_OPENMP "Use OPENMP: kernel convolution on 4 threads: 1 kernel per thread (shared memory)" OFF)
if (UCA_USE_OPENMP)
   set(CMAKE_CXX_FLAGS "-std=c++11 -fopenmp" CACHE STRING "" FORCE)
endif (UCA_USE_OPENMP)

#-------------------------------------------------------
# MPI :
option(UCA_USE_MPI "Use MPI for kernel convolutions on N processes (Distributed Memory => Message Passing)" ON)
if (UCA_USE_MPI)
   find_package(MPI REQUIRED)
endif (UCA_USE_MPI)

#-------------------------------------------------------
# VERBOSE:
option(UCA_VERBOSE "write extra information to console" OFF)

# ------------------------------------------------------
set(UGUCA_INCLUDE_DIRS
  "${CMAKE_CURRENT_BINARY_DIR}/src"

  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/interface_laws"
  )
include_directories(${UGUCA_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_PATH})
include_directories(${FFTW_INCLUDE_PATH})
include_directories(${BLAS_INCLUDE_PATH})
include_directories(${MPI_INCLUDE_PATH})

add_subdirectory(src)

configure_file(cmake/ugucaConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/ugucaConfig.cmake" @ONLY)
configure_file(cmake/ugucaBuildTreeSettings.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/ugucaBuildTreeSettings.cmake" @ONLY)


#===============================================================================
# Release file target
#===============================================================================
set(RELEASE_VFILE ${CMAKE_BINARY_DIR}/release.version)
set(RELEASE_PFILE ${CMAKE_BINARY_DIR}/release.patch)
set(__stamp_file ${CMAKE_CURRENT_BINARY_DIR}/release.stamp)

# get full path to all files
unset(_sources)
foreach(_file ${UGUCA_DEPEND_FILES})
  get_filename_component(_abs_file "${_file}" ABSOLUTE
    BASE_DIR ${PROJECT_SOURCE_DIR}/src)
  if(EXISTS ${_abs_file})
    list(APPEND _sources "${_abs_file}")
  endif()
endforeach()

if(GIT_FOUND)
  add_custom_command(
    OUTPUT ${__stamp_file}
    COMMAND ${CMAKE_COMMAND} -E touch ${__stamp_file}
    COMMAND ${GIT_EXECUTABLE} diff HEAD > ${RELEASE_PFILE}
    COMMAND ${GIT_EXECUTABLE} remote -v > ${RELEASE_VFILE}
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD >> ${RELEASE_VFILE}
    COMMAND ${GIT_EXECUTABLE} status >> ${RELEASE_VFILE}
    COMMENT "Building Release Info"
    DEPENDS ${_sources}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
else()
  add_custom_command(
    OUTPUT ${__stamp_file}
    COMMAND ${CMAKE_COMMAND} -E touch ${__stamp_file}
    COMMAND echo "git is not installed" > ${RELEASE_VFILE}
    )
endif()

add_custom_target(RELEASE_INFO ALL DEPENDS ${__stamp_file})
unset(__stamp_file)

#===============================================================================
# DOCUMENTATION
#===============================================================================
option(UCA_DEVELOPER_DOC "developer documentation" OFF)
if(UCA_DEVELOPER_DOC)
  add_subdirectory(doc)
endif()

#===============================================================================
# TESTS / EXAMPLES / BENCHMARKS / SIMULATIONS
#===============================================================================

# define add_simulation function
function(add_simulation SIM_EXE SIM_SRC)
  add_executable(${SIM_EXE} ${SIM_SRC})
  target_link_libraries(${SIM_EXE} uguca)
endfunction()

# ------------------------------------------------------
# test
option(UCA_TEST "test" ON)
if(${UCA_TEST})
  enable_testing()
  add_subdirectory(test)
endif()

# ------------------------------------------------------
# benchmarks
option(UCA_BENCHMARKS "benchmarks" OFF)
if(UCA_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

# ------------------------------------------------------
# examples
option(UCA_EXAMPLES "examples" ON)
if(UCA_EXAMPLES)
  add_subdirectory(examples)
endif()

# ------------------------------------------------------
# add simulation folder (if you like)
# this folder does not exist
option(UCA_SIMULATIONS "simulations" OFF)
if (UCA_SIMULATIONS)
  add_subdirectory(simulations)
endif()
