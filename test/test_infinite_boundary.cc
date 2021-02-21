/**
 * @file   test_inf_bc.cc
 *
 * @author David S. Kammer <dkammer@ethz.ch>
 * @author Gabriele Albertini <ga288@cornell.edu>
 * @author Chun-Yu Ke <ck659@cornell.edu>
 *
 * @date creation: Tue Jun 30 2020
 * @date last modification: Tue Jun 30 2020
 *
 * @brief  TODO
 *
 *
 * Copyright (C) 2020 ETH Zurich (David S. Kammer)
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
#include <iomanip>      // std::setprecision

#include <cmath>
#include <sstream>
#include <vector>

/* ------------------------------------------------------------------------ */
// HYBRID METHOD

#include "static_communicator_mpi.hh"
#include "material.hh"
#include "uca_mesh.hh"
#include "infinite_boundary.hh"

using namespace uguca;

//define private = public;
//define protected = public;


int main() {

  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();

  //----------------------------------------------------------------//
  // initialize infinite boundary

  //material
  double inf_E   = 5e6;
  double inf_nu  = 0.25;
  double inf_rho = 1e3;

  Material mat = Material(inf_E,inf_nu,inf_rho);
  mat.readPrecomputedKernels();

  //dimension
  int dim = 3;
  unsigned int inf_nb_nodes_x = 3;
  double inf_length_x = 3.0;

  double time_step=0.001;

  int side_factor=1;

  unsigned int inf_nb_nodes_y = 4;
  double inf_length_y = 4.0;

  std::vector<std::vector<double>> coord_tmp;

  if (world_size==1) {
    coord_tmp=
      {{2.0, 0.0, 1.0},  // 9
       {0.0, 0.0, 1.0},  // 1
       {0.0, 0.0, 0.0},  // 0
       {0.0, 0.0, 2.0},  // 2
       {0.0, 0.0, 3.0},  // 3
       {1.0, 0.0, 0.0},  // 4
       {1.0, 0.0, 1.0},  // 5
       {1.0, 0.0, 2.0},  // 6
       {1.0, 0.0, 3.0},  // 7
       {2.0, 0.0, 0.0},  // 8
       {2.0, 0.0, 2.0},  // 10
       {2.0, 0.0, 3.0}}; // 11
  }
  else {
    if (world_rank==0) {
      coord_tmp=
	{{2.0, 0.0, 1.0},   // 9
	 {0.0, 0.0, 1.0},   // 1
	 {0.0, 0.0, 0.0},   // 0
	 {0.0, 0.0, 2.0},   // 2
	 {0.0, 0.0, 3.0}};  // 3
    }
    if (world_rank==1) {
      coord_tmp=
	{{1.0, 0.0, 0.0},   // 4
	 {1.0, 0.0, 1.0},   // 5
	 {1.0, 0.0, 2.0},   // 6
	 {1.0, 0.0, 3.0},   // 7
	 {2.0, 0.0, 0.0},   // 8
	 {2.0, 0.0, 2.0},   // 10
	 {2.0, 0.0, 3.0}};  // 11
    }
  }
  int nb_nodes_local = coord_tmp.size();

  std::vector<NodalField *> coords_local;
  coords_local.resize(dim);
  for (int d=0; d<dim; ++d) {
    coords_local[d] = new NodalField(nb_nodes_local);
    for (int n=0; n<nb_nodes_local; ++n) {
      (*coords_local[d])(n)=coord_tmp[n][d];
    }
  }
  std::cout<<"init mesh"<<std::endl;
  Mesh mesh(inf_length_x,inf_nb_nodes_x,
	    inf_length_y,inf_nb_nodes_y,coords_local);

  std::cout<<"init infinite boundary"<<std::endl;
  InfiniteBoundary infinite_boundary(mesh,
				     side_factor,
				     mat);

  std::cout<<"set time"<<std::endl;
  // time step
  infinite_boundary.setTimeStep(time_step);
  std::cout<<"init predictor corrector"<<std::endl;
  infinite_boundary.initPredictorCorrector();
  std::cout<<"init convolutions"<<std::endl;
  infinite_boundary.initConvolutions();
  

  // populate fields
  for (int i=0; i<dim; ++i) {
    NodalField * inf_bc_ext = infinite_boundary.getExternal(i);
    FFTableNodalField * inf_bc_dsp = infinite_boundary.getDisp(i);
    NodalField * inf_bc_vel = infinite_boundary.getVelo(i);
      
    inf_bc_ext->zeros();
    inf_bc_dsp->zeros();
    for (int n=0; n<inf_bc_vel->getNbNodes(); ++n){
      (*inf_bc_vel)(n)=coord_tmp[n][0]*coord_tmp[n][2]+coord_tmp[n][0]+coord_tmp[n][2];
    }
  }

  std::cout<<"test advance time step Neumann"<<std::endl;
  
  infinite_boundary.advanceTimeStepNeumann();
  double mu = mat.getShearModulus();
  double Cs = mat.getCs();
  double Cp = mat.getCp();
  std::vector<double> eta = {1.0, Cp / Cs, 1.0};
  for (int i=0; i<dim; ++i) {
    NodalField * inf_bc_ext = infinite_boundary.getExternal(i);
    NodalField * inf_bc_vel = infinite_boundary.getVelo(i);
    for (int n=0; n<inf_bc_ext->getNbNodes(); ++n){
      if(std::abs((*inf_bc_ext)(n)-
		  (- side_factor * mu/Cs*eta[i]*(*inf_bc_vel)(n)))>1e-12){
	std::cout<<"error "<<std::endl;
	return 1;
      }
    }
  }
  
  std::cout<<"advance time step Neumann success!"<<std::endl;


  std::cout<<"test advance time step Dirichlet"<<std::endl;
  infinite_boundary.advanceTimeStepDirichlet();
  {
    NodalField * u0 = infinite_boundary.getDisp(0);

    if (false) {
      std::cout<<"solution"<<std::endl
	       << std::setprecision(12)
	       << (*u0)(4) << ", "
	       << (*u0)(9) << std::endl;
    }
    else {
      if (world_size==1) {
	if (std::abs((*u0)(4)- 0.003)>1e-6 ||
	    std::abs((*u0)(9)- 0.002)>1e-6) {
	  std::cout << "failed" << std::endl
		    << (*u0)(4) << ", "
		    << (*u0)(9) << std::endl;
	  return 1; // failure
	}
      }
      else {
	if ((world_rank==0 && std::abs((*u0)(4)- 0.003)>1e-6) ||
	    (world_rank==1 && std::abs((*u0)(4)- 0.002)>1e-6)) {
	  std::cout << "failed" << std::endl;
	  return 1; // failure
	}
      }
    }
  }
  std::cout<<"test predict time step Dirichlet"<<std::endl;
  
  infinite_boundary.predictTimeStepDirichlet();
  
  {
    NodalField * u0 = infinite_boundary.getDisp(0,1);

    if (false) {
      std::cout<<"solution"<<std::endl
	       << std::setprecision(12)
	       << (*u0)(4) << ", "
	       << (*u0)(9) << std::endl;
    }
    else {
      if (world_size==1) {
	if (std::abs((*u0)(4)- (-2.00839640955e-05))>1e-12 ||
	    std::abs((*u0)(9)- (-2.86716308883e-07))>1e-12) {
	  std::cout << "failed" << std::endl
		    << (*u0)(4) << ", "
		    << (*u0)(9) << std::endl;
	  return 1; // failure
	}
      }
      else {
	if ((world_rank==0 && std::abs((*u0)(4)- (-2.00839640955e-05))>1e-12) ||
	    (world_rank==1 && std::abs((*u0)(4)- (-2.86716308883e-07))>1e-12)) {
	  std::cout << "failed" << std::endl;
	  return 1; // failure
	}
      }
    }
  }

  // free memory
  for (int d=0; d<dim; ++d) {
    delete coords_local[d];
  }

  std::cout<<"went til the end rank "<<world_rank<<"\n";

  std::cout << "all checks passed -> overall success" << std::endl;
  
  StaticCommunicatorMPI::getInstance()->finalize();
  return 0;
}
