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
function(uguca_add_benchmark NAME)
  add_subdirectory(${NAME})
  configure_file("plotting_utils.py" "${NAME}/plotting_utils.py" COPYONLY)
  configure_file("${NAME}/inspect_results.py" "${NAME}/inspect_results.py" COPYONLY)
  if(EXISTS ${CMAKE_SOURCE_DIR}/benchmarks/${NAME}/ref)
    file(COPY "${NAME}/ref" DESTINATION "${NAME}")
  endif()
endfunction()

uguca_add_benchmark(TPV3)
uguca_add_benchmark(TPV101)
uguca_add_benchmark(TPV102)
uguca_add_benchmark(TPV103)
uguca_add_benchmark(TPV104)
uguca_add_benchmark(TPV205)
