/**
 * @file   infinite_boundary.cc
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
#include "infinite_boundary.hh"
#include "static_communicator_mpi.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */

InfiniteBoundary::InfiniteBoundary(Mesh & mesh,
				   int side_factor,
				   Material & material) :
  HalfSpace(mesh, side_factor) {

  this->setMaterial(&material);

  this->external.    resize(this->mesh.getDim());
  this->scratch.     resize(this->mesh.getDim());

  for (int d=0; d<this->mesh.getDim(); ++d){
    this->external[d]     = new NodalField(mesh.getNbNodes());
    this->scratch[d]      = new NodalField(mesh.getNbNodes());
  }
}

/* -------------------------------------------------------------------------- */
InfiniteBoundary::~InfiniteBoundary() {
  for (int d=0; d<this->mesh.getDim(); ++d) {
    delete this->external[d];
    delete this->scratch[d];
  }
}


/* -------------------------------------------------------------------------- */
void InfiniteBoundary::gatherCostumMeshForwardFFT(bool predicting) {
  HalfSpace::gatherCostumMeshForwardFFT(this->scratch,
					predicting);
}

/* -------------------------------------------------------------------------- */
void InfiniteBoundary::predictTimeStepDirichlet() {

  // copy v to scratch memory
  this->updateVelocity();
  
  // Predict
  // u* = u + v * dt
  this->computeDisplacement(true);

  // tau_ext* -> compute residual -> tau_res*
  this->computeResidual(this->external);
  // tau_res* -> compute velocity -> v*
  this->computeVelocity(true);

  // Correct
  // v** = (v + v*) / 2 ---> overwrite storage if reached last prediction step
  this->correctVelocity(true);
}
void InfiniteBoundary::advanceTimeStepDirichlet() {

  // copy FEM displacement and internal

  this->computeDisplacement();

  this->gatherCostumMeshForwardFFT();

  this->computeStressFourierCoeff();

  this->backwardFFTscatterCostumMesh();

  this->computeResidual(this->external);

  this->computeVelocity();

}

// --------------------------------------------------------------------------
// NEUMANN BC on FEM
void InfiniteBoundary::computeExternal() {
  double mu = this->material->getShearModulus();
  double Cs = this->material->getCs();
  double Cp = this->material->getCp();
  std::vector<double> eta = {1.0, Cp / Cs, 1.0};

  for (int d = 0; d < this->mesh.getDim(); ++d) {
    double *int_p =  this->internal[d]->storage();
    double *ext_p =  this->external[d]->storage();
    double *velo_p = this->velo[d]->storage();
    double eta_d = eta[d];

    for (int n = 0; n < this->mesh.getNbNodes(); ++n) {
      ext_p[n] = - this->side_factor * mu / Cs * eta_d * velo_p[n] + int_p[n];
    }
  }
}


/* --------------------------------------------------------------------------
 *
 * external = traction to keep gap closed
 *
 * displacements and velocities are use in FEM as Dirichlet boundary condition
 * for the next time step
 *
 */

void InfiniteBoundary::advanceTimeStepNeumann() {

  // set displacement and velocity from FEM

  this->gatherCostumMeshForwardFFT(false);

  this->computeStressFourierCoeff();

  this->backwardFFTscatterCostumMesh();

  this->computeExternal();
}


__END_UGUCA__
