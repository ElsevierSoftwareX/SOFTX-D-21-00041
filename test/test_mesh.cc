/**
 * @file   test_mesh.cc
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
#include "uca_mesh.hh"
#include "static_communicator_mpi.hh"
#include "fftable_nodal_field.hh"

#include <random>
#include <iostream>
#include <unistd.h>

using namespace uguca;

void print(std::vector<NodalField *> & fld){
  int dim = fld.size();
  int size = fld[0]->getNbNodes();
  for (int n=0; n<size;n++) {
    for (int d=0;d<dim;d++) {
      std::cout<<(*fld[d])(n)<<", ";
    }
    std::cout<<std::endl;
  }
  std::cout<<"----"<<std::endl;
}

int main(){//int argc, char *argv[]) {

  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  std::cout<<"world_rank "<<world_rank<<std::endl;

  unsigned int nb_nodes_x = 3;
  double length_x = 3;

  unsigned int nb_nodes_z = 4;
  double length_z = 4.0;
  int dim=3;

  // construct 2D mesh with default coords
  printf("==============\n2D mesh\n==============\n");
  Mesh mesh2d(length_x,10);//nb_nodes_x);

  std::vector<NodalField *>  coords2d = mesh2d.getCoords();
  std::vector<NodalField *>  wavenb2d = mesh2d.getWaveNumbers();
  sleep(world_rank);
  print(coords2d);
  print(wavenb2d);
  if (world_size==0) {
    std::cout<<"test 2d coords\n";

    int nx=mesh2d.getGlobalNbNodesX();
    for (int i=0; i<nx;++i){
      if (std::abs((*coords2d[0])(i)-mesh2d.getLengthX()*i/nx)>1e-6){
	std::cout<<"error0\n";
	return 1;
      }
      if (std::abs((*coords2d[1])(i)-0.0)>1e-6){
	std::cout<<"error1\n";
	return 1;
      }
    }

    std::cout<<"2d coords success!\n";

    std::cout<<std::flush;

    std::cout<<"test 2d wave nb\n";

    nx=mesh2d.getGlobalNbFFTX();
    for (int i=0; i<nx;++i){
      if (std::abs((*wavenb2d[0])(i)- 2*M_PI / mesh2d.getLengthX()*i)>1e-6){
	std::cout<<"error\n";
	return 1;
      }
    }

    std::cout<<"2d wave nb success!\n";
  }

  // construct 3D mesh with default coords
  std::cout<<std::flush;

  printf("==============\n3D mesh\n==============\n");
  Mesh mesh3d(length_x,nb_nodes_x,
	      length_z,nb_nodes_z);


  std::vector<NodalField *>  coords3d = mesh3d.getCoords();
  std::vector<NodalField *>  wavenb3d = mesh3d.getWaveNumbers();
  sleep(world_rank);

  print(coords3d);
  print(wavenb3d);

  if (world_size==0) {
    std::cout<<"test 3d coords\n";

    int nx=mesh3d.getGlobalNbNodesX();
    int nz=mesh3d.getGlobalNbNodesZ();
    for (int i=0; i<nx; ++i){
      for (int j=0; j<nz; ++j){
	int ij=i*nz+j;
	if (std::abs((*coords3d[0])(ij)-mesh3d.getLengthX()*i/nx)>1e-6){
	  std::cout<<"error0\n";
	  return 1;
	}
	if (std::abs((*coords3d[1])(ij)-0.0)>1e-6){
	  std::cout<<"error1\n";
	  return 1;
	}
	if (std::abs((*coords3d[2])(ij)-mesh3d.getLengthZ()*j/nz)>1e-6){
	  std::cout<<"error2\n";
	  return 1;
	}
      }
    }
    std::cout<<"3d coords success!\n";

    std::cout<<std::flush;

    std::cout<<"test 3d wave nb\n";

    nx=mesh3d.getGlobalNbFFTX();
    nz=mesh3d.getGlobalNbFFTZ();
    int f_ny_x;
    if (mesh3d.getDim()==2)
      f_ny_x = nx;
    else
      f_ny_x = nx/2+1;
    for (int i=0; i<nx;++i){
      for (int j=0; j<nz; ++j){
	int ij=i*nz+j;
	if (std::abs((*wavenb3d[0])(ij)- 2*M_PI / mesh3d.getLengthX()*(i-(i/f_ny_x)*nx))>1e-6){
	  std::cout<<"error0\n";
	  return 1;
	}
	if (std::abs((*wavenb3d[2])(ij)- 2*M_PI / mesh3d.getLengthZ()*j)>1e-6){
	  std::cout<<"error2\n";
	  return 1;
	}
      }
    }
    std::cout<<"3d wave nb success!\n";
  }

  // construct 3D mesh with given coords
  std::cout<<std::flush;

  printf("==============\n3D mesh costum\n==============\n");
  if (world_size==1){
    std::cout<<"test check costum mesh "<<std::endl;
    std::vector<std::vector<double>> coord_tmp;

    coord_tmp=
      {{2.0, 0.0, 1.1},  // 9
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

    int nb_nodes_local = coord_tmp.size();

    std::vector<NodalField *> coords_local;
    coords_local.resize(dim);
    for (int d=0; d<dim;d++) {
      coords_local[d] = new NodalField(nb_nodes_local);
      for (int n=0; n<nb_nodes_local; n++) {
	(*coords_local[d])(n)=coord_tmp[n][d];
	//std::cout<<(*coords_local[d])(n)<<", ";
      }
      //std::cout<<std::endl;
    }

    bool caught_exception = true;
    try{
      // construct mesh
      Mesh mesh3dC(length_x,nb_nodes_x,
		   length_z,nb_nodes_z,coords_local);
      std::cout<<"mesh done"<<world_rank<<std::endl;
      caught_exception = false;
    }
    catch (std::runtime_error &e) {
      std::cout << "caught exception -> success" << std::endl;
    }
    if (!caught_exception) {
      std::cout << "failed" << std::endl;
      return 1; // failure
    }

    //for (int d=0; d<dim;d++)
    // delete coords_local[d];
  }// end test check costum mesh

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
	 {1.0, 0.0, 3.0},   // 7
	 {0.0, 0.0, 3.0}};  // 3
    }
    if (world_rank==1) {
      coord_tmp=
	{{1.0, 0.0, 0.0},   // 4
	 {1.0, 0.0, 1.0},   // 5
	 {1.0, 0.0, 2.0},   // 6
	 {2.0, 0.0, 0.0},   // 8
	 {2.0, 0.0, 2.0},   // 10
	 {2.0, 0.0, 3.0}};  // 11
    }
  }
  int nb_nodes_local = coord_tmp.size();
  
  std::vector<NodalField *> coords_local;
  coords_local.resize(dim);
  for (int d=0; d<dim;d++) {
    coords_local[d] = new NodalField(nb_nodes_local);
    for (int n=0; n<nb_nodes_local; n++) {
      (*coords_local[d])(n)=coord_tmp[n][d];
      //std::cout<<(*coords_local[d])(n)<<", ";
    }
    //std::cout<<std::endl;
  }
  // construct mesh
  Mesh mesh3dC(length_x,nb_nodes_x,
	       length_z,nb_nodes_z,coords_local);

  if (world_rank==0) {
    std::vector<NodalField *>  coords3dC = mesh3dC.getCoords();
    std::vector<NodalField *>  wavenb3dC = mesh3dC.getWaveNumbers();
    print(coords3dC);
    print(wavenb3dC);
  }

  std::vector<NodalField *> coords_global;
  coords_global.resize(dim);
  for (int d=0; d<dim;d++) {
    coords_global[d] = new NodalField(nb_nodes_x*nb_nodes_z);
  }

  if (world_size>1) {
    mesh3dC.gatherAndSortCostumNodes(coords_global[0]->storage(),
				     coords_local[0]->storage());
    mesh3dC.gatherAndSortCostumNodes(coords_global[2]->storage(),
				     coords_local[2]->storage());
  }

  std::cout<<std::flush;
  if (world_rank==0 && world_size>1){
    std::cout<<"coords global sorted"<<std::endl;
    print(coords_global);
  }
  // test
  if (world_rank==0 && world_size>1){
    bool test_pass=true;
    int size = coords_global[0]->getNbNodes();
    for (int n=0; n<size;n++) {
      int idx=n/nb_nodes_z;
      int idz=n%nb_nodes_z;
      if ((* coords_global[0])(n) != length_x*idx/nb_nodes_x ||
	  (* coords_global[2])(n) != length_z*idz/nb_nodes_z )
	test_pass=false;
    }
    if (test_pass)
      std::cout<<"Success 3D mesh costum gatherAndSortNodes \n";
    else
      return 1;
  }

  std::vector<NodalField *> coords_local_tmp;
  coords_local_tmp.resize(dim);
  for (int d=0; d<dim;d++) {
    coords_local_tmp[d] = new NodalField(nb_nodes_x*nb_nodes_z);
  }


  if (world_size>1) {
    mesh3dC.sortAndScatterCostumNodes(coords_global[0]->storage(),
				      coords_local_tmp[0]->storage());
    mesh3dC.sortAndScatterCostumNodes(coords_global[2]->storage(),
				      coords_local_tmp[2]->storage());
  }

  sleep(world_rank);
  if (world_size>1) {
    printf("local coords scatter\n");
    print(coords_local_tmp);
  }


  if (world_size>1){
    bool test_pass=true;
    int size = mesh3dC.getNbNodes();
    for (int n=0; n<size;n++) {
      if ((* coords_local[0])(n) != (* coords_local_tmp[0])(n) ||
	  (* coords_local[2])(n) != (* coords_local_tmp[2])(n) ) {
	printf("%g != %g or %g != %g \n",
	       (*coords_local[0])(n), (*coords_local_tmp[0])(n),
	       (*coords_local[2])(n), (*coords_local_tmp[2])(n));
	test_pass=false;
      }
    }
    if (test_pass)
      std::cout<<"Success 3D mesh costum sortAndScatterNodes \n";
    else {
      std::cout<<"error"<<std::endl;
      return 1;
    }
  }

  for (int d=0; d<dim;d++) {
    delete coords_local[d];
    delete coords_global[d];
    delete coords_local_tmp[d];
  }


  std::cout<<"went til the end rank "<<world_rank<<"\n";
  std::cout << "all checks passed -> overall success" << std::endl;
  StaticCommunicatorMPI::getInstance()->finalize();

  return 0;
}

