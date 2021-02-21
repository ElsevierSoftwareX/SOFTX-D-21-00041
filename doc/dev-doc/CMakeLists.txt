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
# @brief  CMakeLists for online documentation
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

set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")
set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")
set(SPHINX_OUTPUT "${SPHINX_HTML_DIR}/index.html")


find_package(Sphinx REQUIRED)

set(SPHINX_VERBOSE_FLAG "-q")
if(CMAKE_VERBOSE_MAKEFILE)
  set(SPHINX_VERBOSE_FLAG)
endif(CMAKE_VERBOSE_MAKEFILE)

set(UGUCA_IN_READTHEDOC FALSE)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
  @ONLY)


set(_sphinx_command ${SPHINX_BUILD_EXECUTABLE}
  -j auto
  ${SPHINX_VERBOSE_FLAG} -b html
  -c "${CMAKE_CURRENT_BINARY_DIR}"
  -d "${SPHINX_CACHE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${SPHINX_HTML_DIR}"
  )

file(GLOB_RECURSE _SPHINX_SRCS
  "*.rst")

add_custom_command(
  OUTPUT ${SPHINX_OUTPUT}
  COMMAND ${_sphinx_command}
  DEPENDS ${SPHINX_INPUT} ${_SPHINX_SRCS} #doxygen-doc
  COMMENT "Building HTML documentation with Sphinx in ${SPHINX_HTML_DIR}"
  )

add_custom_target(sphinx-doc ALL
  DEPENDS ${SPHINX_OUTPUT})

add_custom_target(sphinx-doc-forced
  COMMAND ${_sphinx_command}
#  DEPENDS doxygen-doc
  COMMENT "Building HTML documentation with Sphinx (forced) in ${SPHINX_HTML_DIR}"
  )
