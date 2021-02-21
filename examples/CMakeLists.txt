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

# basic example
add_simulation(basic_example basic_example.cc)
configure_file(basic_example_plot.py basic_example_plot.py COPYONLY)

# 2d fracture example
add_simulation(fracture_2d_example fracture_2d_example.cc)
configure_file(fracture_2d_example.in fracture_2d_example.in COPYONLY)
configure_file(fracture_2d_example.sh fracture_2d_example.sh COPYONLY)
configure_file(fracture_2d_example_plot.py fracture_2d_example_plot.py COPYONLY)
configure_file(fracture_2d_example_fig.py fracture_2d_example_fig.py COPYONLY)

# 3d fracture example
add_simulation(fracture_3d_example fracture_3d_example.cc)
configure_file(fracture_3d_example.in fracture_3d_example.in COPYONLY)
configure_file(fracture_3d_example.sh fracture_3d_example.sh COPYONLY)
configure_file(fracture_3d_example_plot.py fracture_3d_example_plot.py COPYONLY)

# 2d friction example
add_simulation(friction_2d_example friction_2d_example.cc)
configure_file(friction_2d_example.in friction_2d_example.in COPYONLY)
configure_file(friction_2d_example.sh friction_2d_example.sh COPYONLY)
configure_file(friction_2d_example_plot.py friction_2d_example_plot.py COPYONLY)
