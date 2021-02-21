/**
 * @file   test_half_space.cc
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
#include <iomanip>      // std::setprecision

#include "half_space.hh"
#include "static_communicator_mpi.hh"

using namespace uguca;

int main(){

  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  std::cout << "start test: test_half_space" << std::endl;

  /* no need to test:
     initConvolutions
     forwardFFT
     backwardFFT
     gatherCostumMeshForwardFFT
     backwardFFTscatterCostumMesh
  */

  // --------------------------------------------------------------
  // init 2D half_space

  {
    double length = 2.0;
    int nb_elements = 16;

    Mesh msh(length,nb_elements);
    Material mat = Material(71e9, 0.33, 2777);
    mat.readPrecomputedKernels();

    HalfSpace hs2(msh,1);

    hs2.setMaterial(&mat);



    // --------------------------------------------------------------

    std::cout << "check   getStableTimeStep " << std::endl;

    double delta_x = msh.getDeltaX();
    double delta_z = msh.getDeltaZ();

    if (msh.getDim()==2)
      delta_z = 1e100;

    if (std::abs(hs2.getStableTimeStep() -
		 std::min(delta_x,delta_z) / mat.getCs())>1e-12) {
      std::cout << "wrong stable time step" << std::endl;
      return 1; // failure
    }


    double dt=hs2.getStableTimeStep();
    bool caught_exception = true;

    try {
      hs2.setTimeStep(2*dt);
      caught_exception = false;
    }
    catch (std::runtime_error &e) {
      std::cout << "caught exception -> success" << std::endl;
    }
    if (!caught_exception) {
      std::cout << "failed" << std::endl;
      return 1; // failure
    }

    std::cout << "check set/get TimeStep" << std::endl;

    dt*=0.25;
    hs2.setTimeStep(dt);

    std::cout << "check computeDisplacement" << std::endl;

    NodalField * v0 = hs2.getVelo(0);
    NodalField * v1 = hs2.getVelo(1);

    for (int i=0; i<nb_elements; ++i){
      (*v0)(i)=0.3*cos(i*3)+sin(i*2);
      (*v1)(i)=0.5*cos(i*6)+0.4*(sin(i));
      // std::cout << (*v0)(i)<<std::endl;
    }

    hs2.computeDisplacement();
    NodalField * u0 = hs2.getDisp(0);
    NodalField * u1 = hs2.getDisp(1);

    for (int i=0; i<nb_elements; ++i){
      if (std::abs((*u0)(i) - dt*(*v0)(i))>1e-12) {
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      if (std::abs((*u1)(i) - dt*(*v1)(i))>1e-12) {
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }

    std::cout << "check initPredictorCorrector" << std::endl;

    hs2.initPredictorCorrector();
    NodalField * up0 = hs2.getDisp(0,true);
    NodalField * up1 = hs2.getDisp(1,true);
    NodalField * vp0 = hs2.getVelo(0,true);
    NodalField * vp1 = hs2.getVelo(1,true);

    std::cout << "initPredictorCorrector success" << std::endl;

    std::cout << "check updateVelocity" << std::endl;

    hs2.updateVelocity();
    for (int i=0; i<nb_elements; ++i){
      if (std::abs((*vp0)(i) - (*v0)(i))>1e-12) {
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      if (std::abs((*vp1)(i) - (*v1)(i))>1e-12) {
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }


    std::cout << "check computeDisplacement pred" << std::endl;

    hs2.computeDisplacement(true);

    for (int i=0; i<nb_elements; ++i){
      if (std::abs((*up0)(i) - dt*(*v0)(i)-(*u0)(i))>1e-12) {
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
      if (std::abs((*up1)(i) - dt*(*v1)(i)-(*u1)(i))>1e-12) {
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }

    std::cout << "check correct velocity" << std::endl;

    for (int i=0; i<nb_elements; ++i){
      (*v0)(i)*=3;
      (*v1)(i)*=3;
    }
    hs2.correctVelocity(false);
    for (int i=0; i<nb_elements; ++i){
      if (std::abs((*vp0)(i)/2-(*v0)(i)/3)>1e-12) {
	  std::cout << "failed" << std::endl;
	  return 1; // failure
      }
      if (std::abs((*vp1)(i)/2-(*v1)(i)/3)>1e-12) {
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }

    std::cout << "check correct velocity last" << std::endl;

    hs2.correctVelocity(true);
    for (int i=0; i<nb_elements; ++i){
      if (std::abs((*vp0)(i)/2-(*v0)(i)/2.5)>1e-12) {
	  std::cout << "failed" << std::endl;
	  return 1; // failure
      }
      if (std::abs((*vp1)(i)/2-(*v1)(i)/2.5)>1e-12) {
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }


    std::cout << "check initConvolutions" << std::endl;

    hs2.initConvolutions();

    std::cout << "initConvolutions success" << std::endl;

    std::cout << "check computeStressFourierCoeff 2D" << std::endl;

    hs2.forwardFFT();
    hs2.computeStressFourierCoeff();

    FFTableNodalField * s0 = hs2.getInternal(0);
    FFTableNodalField * s1 = hs2.getInternal(1);

    if (world_rank==0) { // real space computations are on 0 rank process
      if (false) {
	std::cout<<"solution"<<std::endl
		 << std::setprecision(12) << s0->fd(4)[0] << ", " << s0->fd(4)[1] << std::endl
		 << s1->fd(2)[0] << ", " << s1->fd(2)[1] << std::endl;
      }
      else {
	if (std::abs(s0->fd(4)[0]- (-138070.216931))>1e-6 ||
	    std::abs(s0->fd(4)[1]- (731156.525273))>1e-6) {
	  std::cout << "failed" << std::endl
		    << s0->fd(4)[0] << ", " << s0->fd(4)[1] << std::endl;

	  return 1; // failure
	}
	if (std::abs(s1->fd(2)[0]- (-44634.1170066))>1e-6 ||
	    std::abs(s1->fd(2)[1]- (7529.73118974))>1e-6) {
	  std::cout << "failed" << std::endl
		    << s1->fd(2)[0] << ", " << s1->fd(2)[1] << std::endl;

	  return 1; // failure
	}
      }
    }

    std::cout << "check computeResidual" << std::endl;

    std::vector<NodalField *> u= {u0,u1};
    hs2.computeResidual(u);

    NodalField * r0 = hs2.getResidual(0);
    NodalField * r1 = hs2.getResidual(1);

    for (int i=0;i<nb_elements;++i){
      if (std::abs((*r1)(i)-((*s1)(i)+(*u1)(i)))>1e-12) {
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }

    std::cout << "check computeVelocity" << std::endl;
    hs2.computeVelocity();
    if (world_rank==0) // real space computations are on 0 rank process
    for (int i=0;i<nb_elements;++i){
      if (std::abs((*v0)(i) - ( mat.getCs()/mat.getShearModulus()*(*r0)(i)) )>1e-12) {
	std::cout << "failed 0" << std::endl
		  << (*v0)(i) << " != " << ( mat.getCs()/mat.getShearModulus()*(*r0)(i)) << std::endl;
	return 1; // failure
      }
      if (std::abs((*v1)(i) - ( mat.getCs()/mat.getShearModulus()/(mat.getCp()/mat.getCs())*(*r1)(i)) )>1e-12) {
	std::cout << "failed 1" << std::endl
		  << (*v1)(i) << " != " << ( mat.getCs()/mat.getShearModulus()/ (mat.getCp()/mat.getCs()) *(*r1)(i) ) << std::endl;
	return 1; // failure
      }
    }



  // --------------------------------------------------------------

  }
  {
    std::cout << "check computeStressFourierCoeff 3D" << std::endl;
    // init 3D half_space

    double lx = 2.0;
    double lz = 1.5;
    int nb_x = 16;
    int nb_z = 8;
    std::cout<<"init mesh 3D" << std::endl;
    Mesh msh3(lx,nb_x,
	      lz,nb_z);
    Material mat = Material(71e9, 0.33, 2777);
    mat.readPrecomputedKernels();

    HalfSpace hs3(msh3,1);

    hs3.setMaterial(&mat);
    double dt=hs3.getStableTimeStep()*0.1;

    hs3.setTimeStep(dt);

    hs3.initPredictorCorrector();
    hs3.initConvolutions();

    std::cout<<"init disp"<<std::endl;

    NodalField * u0 = hs3.getDisp(0);
    NodalField * u1 = hs3.getDisp(1);
    NodalField * u2 = hs3.getDisp(1);

    for (int i=0; i<nb_x; ++i){
      for (int j=0; j<nb_z; ++j){
	int ij =i*nb_z+j;
	(*u0)(ij)=0.3*cos(i*3+6)+sin(i*2+1) +2*cos(j*2)+sin(j*6-2);
	(*u1)(ij)=0.5*cos(i*6)+0.4*(sin(i+2))-1*cos(j*7)+sin(j*9-5);
	(*u2)(ij)=0.5*cos(i*6+5)+0.4*(sin(i-2))+0.5*cos(5-j*2)+sin(j);;

	// std::cout << (*u0)(ij)<<std::endl;
      }
    }

    hs3.forwardFFT();
    hs3.computeStressFourierCoeff();

    std::cout<<"test computeStressFourierCoeff 3D"<<std::endl;
    if (world_rank==0) // complete data is gathered to process 0
    {
      FFTableNodalField * s0 = hs3.getInternal(0);
      FFTableNodalField * s1 = hs3.getInternal(1);
      FFTableNodalField * s2 = hs3.getInternal(2);

      if (false) {
	std::cout<<"solution"<<std::endl
		 << std::setprecision(12)
		 << s0->fd(4)[0]  << ", " << s0->fd(4)[1]  << std::endl
		 << s1->fd(62)[0] << ", " << s1->fd(62)[1] << std::endl
		 << s2->fd(47)[0] << ", " << s2->fd(47)[1] << std::endl;
      }
      else {
	if (std::abs(s0->fd(4)[0]- (-1.12212163586e+12))>1e3 ||
	    std::abs(s0->fd(4)[1]- (0))>1e-3) {
	  std::cout << "failed" << std::endl
		    << s0->fd(4)[0] << ", " << s0->fd(4)[1] << std::endl;

	  return 1; // failure
	}
	if (std::abs(s1->fd(62)[0]- (-1.41414746235e-05))>1e-12 ||
	    std::abs(s1->fd(62)[1]- (1.17845621863e-05))>1e-12) {
	  std::cout << "failed" << std::endl
		    << s1->fd(62)[0] << ", " << s1->fd(62)[1] << std::endl;

	  return 1; // failure
	}
	if (std::abs(s2->fd(47)[0]- (-0.000161704136674))>1e-12 ||
	    std::abs(s2->fd(47)[1]- (4.27856693215e-05))>1e-12) {
	  std::cout << "failed" << std::endl
		    << s2->fd(47)[0] << ", " << s1->fd(47)[1] << std::endl;

	  return 1; // failure
	}
      }
    }
  }

  std::cout << "all checks passed -> overall success" << std::endl;
  StaticCommunicatorMPI::getInstance()->finalize();
  return 0; // success
}
