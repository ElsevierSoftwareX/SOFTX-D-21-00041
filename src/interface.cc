/**
 * @file   interface.cc
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
#include "interface.hh"
#include "static_communicator_mpi.hh"
#include "cmath"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */

Interface::Interface(Mesh & mesh,
		     InterfaceLaw & law) :
  Dumper(mesh),
  mesh(mesh),
  time_step(0.),
  law(law) {

  this->load.resize(this->mesh.getDim());
  this->cohesion.resize(this->mesh.getDim());
  this->scratch_field.resize(this->mesh.getDim());

  for (int d = 0; d < this->mesh.getDim(); ++d) {
    this->load[d]          = new NodalField(this->mesh.getNbNodes());
    this->cohesion[d]      = new NodalField(this->mesh.getNbNodes());
    this->scratch_field[d] = new NodalField(this->mesh.getNbNodes());
  }

  this->law.setInterface(this);
}

/* -------------------------------------------------------------------------- */
Interface::~Interface() {
  for (int d = 0; d < this->mesh.getDim(); ++d) {
    delete this->load[d];
    delete this->cohesion[d];
    delete this->scratch_field[d];
  }
}

/* -------------------------------------------------------------------------- */
void Interface::init(bool velocity_initial_conditions) {
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  for (unsigned int i=0;i<this->half_space.size();++i)
    this->half_space[i]->initConvolutions();

  // like a typical time step (advanceTimeStep)
  // but displacement is potentially already imposed as initial condition
  bool exec_spatial_domain=true;
  if (this->mesh.isSpatialDomainParallel()==false)
    if (world_rank!=0)      // skip this part if spatial domain serial and not root
      exec_spatial_domain=false;

  if (exec_spatial_domain) {
    for (int i = 0; i < this->nb_pc; ++i) {
      if (i == 0)
        this->updateVelocity();
      this->computeDisplacement(true);
      this->computeCohesion(true);
      this->computeResidual();
      if (!velocity_initial_conditions) {
        this->computeVelocity(true);
        this->correctVelocity(i == this->nb_pc - 1);
      }
    }
  }
  this->forwardFFT();
  this->computeStressFourierCoeff();

  if (exec_spatial_domain) {
    this->backwardFFT();
    this->computeCohesion();
    this->computeResidual();
    if (!velocity_initial_conditions)
      this->computeVelocity();
  }
}

/* -------------------------------------------------------------------------- */
void Interface::initDump(const std::string &bname,
			 const std::string &path,
                         const Dumper::Format format) {

  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  if (world_rank!=0) return;

  Dumper::initDump(bname, path, format);
}

/* -------------------------------------------------------------------------- */
void Interface::initPredictorCorrector(int iterations) {
  this->nb_pc = iterations;
  if (iterations>0)
    for (unsigned int i=0;i<this->half_space.size();++i)
      this->half_space[i]->initPredictorCorrector();
}

/* -------------------------------------------------------------------------- */
void Interface::setTimeStep(double time_step) {
  this->time_step = time_step;
  for (unsigned int i=0;i<this->half_space.size();++i)
    this->half_space[i]->setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
double Interface::getStableTimeStep() {

  if (this->half_space.size()==2)
    return std::min(this->half_space[0]->getStableTimeStep(),
		    this->half_space[1]->getStableTimeStep());
  else
    return this->half_space[0]->getStableTimeStep();
}

/* -------------------------------------------------------------------------- */
void Interface::setDynamic(bool fully_dynamic) {
  for (unsigned int i = 0; i < this->half_space.size(); ++i)
    this->half_space[i]->setDynamic(fully_dynamic);
}

/* -------------------------------------------------------------------------- */
void Interface::computeDisplacement(bool predicting) {
  for (unsigned int i=0;i<this->half_space.size();++i)
    this->half_space[i]->computeDisplacement(predicting);
}
/* -------------------------------------------------------------------------- */
void Interface::forwardFFT(bool predicting) {
  for (unsigned int i=0;i<this->half_space.size();++i)
    this->half_space[i]->forwardFFT(predicting);
}

/* -------------------------------------------------------------------------- */
void Interface::computeStressFourierCoeff(bool predicting, bool correcting) {
  for (unsigned int i=0;i<this->half_space.size();++i)
    this->half_space[i]->computeStressFourierCoeff(predicting, correcting);
}

/* -------------------------------------------------------------------------- */
void Interface::backwardFFT() {
  for (unsigned int i=0;i<this->half_space.size();++i)
    this->half_space[i]->backwardFFT();
}

/* -------------------------------------------------------------------------- */
void Interface::computeResidual() {
  // doesn't matter if predicting or not
  std::vector<NodalField *> load_and_cohesion = this->scratch_field;
  this->combineLoadAndCohesion(load_and_cohesion);
  for (unsigned int i=0;i<this->half_space.size();++i)
    this->half_space[i]->computeResidual(load_and_cohesion);
}
/* -------------------------------------------------------------------------- */
void Interface::computeVelocity(bool predicting) {
  for (unsigned int i=0;i<this->half_space.size();++i)
    this->half_space[i]->computeVelocity(predicting);
}

/* -------------------------------------------------------------------------- */
void Interface::gatherCostumMeshForwardFFT() {
  std::vector<NodalField *> disp_copy = this->scratch_field;
  for (unsigned int i=0;i<this->half_space.size();++i)
    this->half_space[i]->gatherCostumMeshForwardFFT(disp_copy);
}

/* -------------------------------------------------------------------------- */
void Interface::backwardFFTscatterCostumMesh() {
  for (unsigned int i=0;i<this->half_space.size();++i)
    this->half_space[i]->backwardFFTscatterCostumMesh();
}

/* -------------------------------------------------------------------------- */
void Interface::updateVelocity() {
  for (unsigned int i=0;i<this->half_space.size();++i)
    this->half_space[i]->updateVelocity();
}

/* -------------------------------------------------------------------------- */
void Interface::correctVelocity(bool last_step) {
  for (unsigned int i=0;i<this->half_space.size();++i)
    this->half_space[i]->correctVelocity(last_step);
}


/* -------------------------------------------------------------------------- */
void Interface::advanceTimeStep() {
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  bool exec_spatial_domain=true;
  if (this->mesh.isSpatialDomainParallel()==false)
    if (world_rank !=0 ) // skip this part if spatial domain serial and not root
      exec_spatial_domain=false;


  if (exec_spatial_domain) {
    // predictor-corrector
    for (int i = 0; i < this->nb_pc; ++i) {
      if (i == 0) {
        // copy v to scratch memory
        this->updateVelocity();
      }
      // Predict
      // u* = u + v * dt
      this->computeDisplacement(true);

      // Evaluate

      // this->forwardFFT(true);
      // this->computeStressFourierCoeff(true, i > 0);
      // this->backwardFFT();

      // f* -> compute cohesion -> tau_coh*
      this->computeCohesion(true);
      // tau_coh* -> compute residual -> tau_res*
      this->computeResidual();
      // tau_res* -> compute velocity -> v*
      this->computeVelocity(true);

      // Correct
      // v** = (v + v*) / 2 ---> overwrite storage if reached last step
      this->correctVelocity(i == this->nb_pc - 1);
    }

    // compute displacement
    this->computeDisplacement();
  }

  // to fourier space
  // if fftw_mpi executed by all processes else automatically by root
  if (this->mesh.isCostum())
    this->gatherCostumMeshForwardFFT();
  else
    this->forwardFFT();

  // compute convolutions
  //this->computeStressFourierCoeff(false, this->nb_pc > 0); // use this if evaluated when predicting
  this->computeStressFourierCoeff(false,false);

  // back to normal space
  // if fftw_mpi executed by all processes else automatically by root
  if (this->mesh.isCostum())
    this->backwardFFTscatterCostumMesh();
  else
    this->backwardFFT();

  if (exec_spatial_domain) {
    // compute forces due to interface law
    this->computeCohesion();

    // compute residual force
    this->computeResidual();

    // compute velocity
    this->computeVelocity();
  }
}
/* -------------------------------------------------------------------------- */
void Interface::computeCohesion(bool predicting) {
  this->law.computeCohesiveForces(this->cohesion, predicting);
}

/* -------------------------------------------------------------------------- */
void Interface::computeNorm(const std::vector<NodalField *> & field,
			    NodalField & norm, bool shear_only) {

  norm.zeros();
  double * norm_p = norm.storage();

  for (int d=0; d<this->mesh.getDim(); ++d) {
    if (shear_only && d == 1) continue;

    double * field_d_p = field[d]->storage();

    for (int n=0;n<this->mesh.getNbNodes(); ++n) {
      norm_p[n] += field_d_p[n] * field_d_p[n];
    }
  }

  for (int n=0;n<this->mesh.getNbNodes(); ++n) {
    norm_p[n] = std::sqrt(norm_p[n]);
  }
}

/* -------------------------------------------------------------------------- */
void Interface::multiplyFieldByScalar(std::vector<NodalField *> & field,
				      const NodalField & scalar,
				      bool shear_only) {

  const double * scalar_p = scalar.storage();

  for (int d=0; d<this->mesh.getDim(); ++d) {
    if (shear_only && d == 1) continue;

    double * field_d_p = field[d]->storage();

    for (int n=0;n<this->mesh.getNbNodes(); ++n) {
      field_d_p[n] *= scalar_p[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void Interface::combineLoadAndCohesion(std::vector<NodalField *> &load_and_cohesion) {
  for (int d=0; d<this->mesh.getDim(); ++d) {
    // tau_0 - tau_coh
    double * load_and_cohesion_p = load_and_cohesion[d]->storage();
    double * coh_p = this->cohesion[d]->storage();
    double * load_p = this->load[d]->storage();

    for (int n=0; n<this->mesh.getNbNodes(); ++n) {
      load_and_cohesion_p[n] = load_p[n] - coh_p[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void Interface::registerDumpField(const std::string & field_name) {

  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

   if (world_rank!=0) return;

   int d = std::atoi(&field_name[field_name.length()-1]);

  if (d >= this->mesh.getDim())
    throw std::runtime_error("Field "+field_name
			     +" cannot be dumped, too high dimension");

  if (field_name == "load_"+std::to_string(d))
    this->registerForDump(field_name,
			  this->load[d]);

  else if (field_name == "cohesion_"+std::to_string(d))
    this->registerForDump(field_name,
			  this->cohesion[d]);

  else
    this->law.registerDumpField(field_name);

}

__END_UGUCA__
