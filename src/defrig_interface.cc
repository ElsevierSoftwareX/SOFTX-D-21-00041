/**
 * @file   defrig_interface.cc
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
#include "defrig_interface.hh"
#include "static_communicator_mpi.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
DefRigInterface::DefRigInterface(Mesh & mesh,
				 Material & top_material,
				 InterfaceLaw & law) :
  Interface(mesh, law),
  top(mesh, 1)
{
  this->half_space.resize(1);
  this->half_space[0] = &this->top;
  this->top.setMaterial(&top_material);

  // top material information
  const Material & mat_t = this->top.getMaterial();
  double cp_t = mat_t.getCp();
  double cs_t = mat_t.getCs();
  double mu_t = mat_t.getShearModulus();
  double eta_t = cp_t / cs_t;
  this->fact_t_2 = cs_t / mu_t / eta_t;
}

/* -------------------------------------------------------------------------- */
DefRigInterface::~DefRigInterface() {}

/* -------------------------------------------------------------------------- */
void DefRigInterface::closingNormalGapForce(NodalField * close_force, bool predicting) {
  // C factor of notes
  double fact_t = this->time_step * this->fact_t_2;

  // accessors
  double * f_1_t = this->top.getInternal(1)->storage();
  double * t0_1 = this->load[1]->storage();
  double * cf = close_force->storage();

  std::vector<NodalField *> gap = this->scratch_field;
  this->computeGap(gap, predicting);
  double * gap_1_p = gap[1]->storage();

  for (int n=0; n<this->mesh.getNbNodes(); ++n) {
    double u_1_gap = gap_1_p[n] / fact_t;
    double du_1_t = t0_1[n] + f_1_t[n];
    cf[n] = u_1_gap + du_1_t;
  }
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::maintainShearGapForce(std::vector<NodalField *> &maintain_force) {
  // doesn't matter if predicting or not
  for (int d=0; d<this->mesh.getDim();d+=2) {
    // accessors
    double * f_t = this->top.getInternal(d)->storage();
    double * t0 = this->load[d]->storage();
    double * mf = maintain_force[d]->storage();

    for (int n=0; n<this->mesh.getNbNodes(); ++n) {
      mf[n] = t0[n] + f_t[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::computeGap(std::vector<NodalField *> &gap,
                                 bool predicting) {
  for (int d=0;d<this->mesh.getDim(); ++d) {
    double * top_disp = this->top.getDisp(d, predicting)->storage();
    double * gap_p = gap[d]->storage();

    for (int n=0; n<this->mesh.getNbNodes(); ++n) {
      gap_p[n] = top_disp[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::computeGapVelocity(std::vector<NodalField *> &gap_velo,
                                         bool predicting) {
  for (int d = 0; d < this->mesh.getDim(); ++d) {
    double *top_velo = this->top.getVelo(d, predicting)->storage();
    double *gap_velo_p = gap_velo[d]->storage();

    for (int n = 0; n < this->mesh.getNbNodes(); ++n) {
      gap_velo_p[n] = top_velo[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::registerDumpField(const std::string & field_name) {
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  if (world_rank!=0) return;

  int d = std::atoi(&field_name[field_name.length()-1]);

  if (d >= this->mesh.getDim())
    throw std::runtime_error("Field "+field_name
			     +" cannot be dumped, too high dimension");

  bool registered = false;
  // field_name starts with "top"
  if (field_name.rfind("top", 0) == 0) {
    // cut away "top_" from field_name and give interface as dumper
    registered = this->top.registerDumpFieldToDumper(field_name.substr(4),
						     field_name,
						     this);
  }

  if (!registered)
    Interface::registerDumpField(field_name);
}

__END_UGUCA__
