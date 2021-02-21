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
include_directories(${uguca_SOURCE_DIR}/src)

function(uguca_add_test TEST_NAME)
  add_simulation(${TEST_NAME} "${TEST_NAME}.cc")
  # if there is a .sh file use it for test (not for mpi though)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.sh")
    configure_file("${TEST_NAME}.sh" "${TEST_NAME}.sh" COPYONLY)
    add_test(${TEST_NAME} "${TEST_NAME}.sh")
  else()
    add_test(${TEST_NAME} ${TEST_NAME})
  endif()
  # if mpi is on and there is a .mpi.sh file use it for test
  if(UCA_USE_MPI AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.mpi.sh")
    configure_file("${TEST_NAME}.mpi.sh" "${TEST_NAME}.mpi.sh" COPYONLY)
    add_test("${TEST_NAME}-mpi" "${TEST_NAME}.mpi.sh")
  endif()
  # if file needed to verify solution
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.txt")
    configure_file("${TEST_NAME}.txt" "${TEST_NAME}.txt" COPYONLY)
  endif()
endfunction()

# basics
uguca_add_test(test_parameter_reader)
uguca_add_test(test_nodal_field)
uguca_add_test(test_fftable_nodal_field)
uguca_add_test(test_limited_history)
uguca_add_test(test_mesh)
if(UCA_USE_MPI)
  uguca_add_test(test_static_communicator_mpi)
endif()

# kernel
uguca_add_test(test_precomputed_kernel)
uguca_add_test(test_kernel_collection)
uguca_add_test(test_preint_kernel)

# core
uguca_add_test(test_material)
uguca_add_test(test_dumper)
uguca_add_test(test_half_space)
uguca_add_test(test_infinite_boundary)

# test interfaces
uguca_add_test(test_interface)
uguca_add_test(test_unimat_shear_interface)
uguca_add_test(test_bimat_interface)
uguca_add_test(test_defrig_interface)

# test interface laws
uguca_add_test(test_linear_shear_cohesive_law)
uguca_add_test(test_linear_coulomb_friction_law)
uguca_add_test(test_barras_law)
uguca_add_test(test_rate_and_state_law)
