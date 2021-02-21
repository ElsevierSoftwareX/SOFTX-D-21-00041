/**
 * @file   test_fftable_nodal_field.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Fri Feb 5 2021
 * @date last modification: Fri Feb 5 2021
 *
 * @brief  TODO
 *
 *
 * Copyright (C) 2021 ETH Zurich (David S. Kammer)
 *
 * This file is part of uguca.
 *
 * uguca is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * uguca is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with uguca.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <math.h>

#include "fftable_nodal_field.hh"
#include "static_communicator_mpi.hh"

using namespace uguca;

/* This is serial
 * however it can be executed in parallel and only the 0 rank would perform the operations
 */

int main() {

  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  if (world_rank==0)
    std::cout << "start test: fftable_nodal_field" << std::endl;

  int nb_points_x = 64*2;
  int nb_points_z = 64;


  double length_x = 3.0;
  double length_z = 1.0;


  Mesh mesh(length_x,nb_points_x,
	    length_z,nb_points_z);

  FFTableNodalField ftf(mesh);
  if (world_rank==0)
    std::cout<<"fftable field contructed"<< std::endl;

  NodalField field(nb_points_x*nb_points_z);

  std::vector<NodalField *>  coords = mesh.getCoords();

  // initiate field and ftf
  for (int i=0; i<mesh.getGlobalNbNodesX(); ++i) {
    for (int j=0; j<mesh.getGlobalNbNodesZ(); ++j) {
      int ij = i*mesh.getGlobalNbNodesZ()+j;
      double x = (*coords[0])(ij);
      double z = (*coords[2])(ij);
      if (world_rank!=0) continue;
      field(ij) = 1.0+cos(x*4*M_PI)*sin(z*8*M_PI)+cos(x*2*M_PI);
      ftf(ij) = field(ij);
    }
  }

  if (world_rank==0)
    std::cout << "fftable field initialized" << std::endl
	      << "check forward and backward fft" << std::endl;

  ftf.forwardFFT();
  ftf.setAllValuesTo(0.0);
  ftf.backwardFFT();

  // check
  for (int i=0; i<mesh.getGlobalNbNodesX(); ++i) {
    for (int j=0; j<mesh.getGlobalNbNodesZ(); ++j) {
      int ij = i*mesh.getGlobalNbNodesZ()+j;
      if (world_rank!=0) continue;
      if (fabs(field(ij) - ftf(ij))>1e-12) {
	std::cout <<"("<<i<<","<<j<<") "<< fabs(field(ij) - ftf(ij)) << "\n";
	std::cout <<"2D FFT failed"<< std::endl;
	return 1;
      }
    }
  }

  if (world_rank==0)
  std::cout << "forward and backward FFT correct"<< std::endl
	    << "all checks passed -> overall success" << std::endl;

  StaticCommunicatorMPI::getInstance()->finalize();
  return 0; // success
}
