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
set(UGUCA_HEADER
  uca_common.hh
  uca_parameter_reader.hh
  material.hh
  nodal_field.hh
  fftable_nodal_field.hh
  limited_history.hh
  kernel.hh
  precomputed_kernel.hh
  kernel_collection.hh
  preint_kernel.hh
  half_space.hh
  uca_dumper.hh
  uca_mesh.hh

  # MPI static communicator
  static_communicator_mpi.hh

  # use method as infinite boundary condition for FEM
  infinite_boundary.hh

  # interfaces
  interface.hh
  bimat_interface.hh
  defrig_interface.hh
  unimat_shear_interface.hh

  # interface laws
  interface_law.hh
  interface_laws/linear_coulomb_friction_law.hh
  interface_laws/linear_shear_cohesive_law.hh
  interface_laws/barras_law.hh
  interface_laws/rate_and_state_law.hh
  )

set(UGUCA_SRC
  uca_parameter_reader.cc
  material.cc
  nodal_field.cc
  fftable_nodal_field.cc
  limited_history.cc
  kernel.cc
  precomputed_kernel.cc
  kernel_collection.cc
  preint_kernel.cc
  half_space.cc
  uca_dumper.cc
  uca_mesh.cc

  # MPI static communicator
  static_communicator_mpi.cc

  # use method as infinite boundary condition for FEM
  infinite_boundary.cc

  # interfaces
  interface.cc
  bimat_interface.cc
  defrig_interface.cc
  unimat_shear_interface.cc

  # interface laws
  interface_law.cc
  interface_laws/linear_coulomb_friction_law.cc
  interface_laws/linear_shear_cohesive_law.cc
  interface_laws/barras_law.cc
  interface_laws/rate_and_state_law.cc
  )
set(UGUCA_DEPEND_FILES ${UGUCA_SRC} ${UGUCA_HEADER} PARENT_SCOPE)

add_library(uguca ${UGUCA_SRC})

target_link_libraries(uguca ${MPI_LIBRARIES} ${FFTW_LIBRARIES}
  ${GSL_LIBRARIES} ${BLAS_LIBRARIES})


list(APPEND UGUCA_EXPORT_LIST uguca)

if("${UGUCA_TARGETS_EXPORT}" STREQUAL "ugucaLibraryDepends")
  install(EXPORT ugucaLibraryDepends DESTINATION  lib/uguca
    COMPONENT dev)
endif()

#Export for build tree
export(TARGETS ${UGUCA_EXPORT_LIST}
  FILE "${CMAKE_BINARY_DIR}/ugucaLibraryDepends.cmake")
export(PACKAGE uguca)
