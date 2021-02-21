/**
 * @file   bimat_interface.cc
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
#include "bimat_interface.hh"
#include "static_communicator_mpi.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */

BimatInterface::BimatInterface(Mesh & mesh,
			       Material & top_material,
			       Material & bot_material,
			       InterfaceLaw & law) :
  Interface(mesh, law),
  top(mesh,  1),
  bot(mesh, -1)
{
  this->half_space.resize(2);
  this->half_space[0] = &this->top;
  this->half_space[1] = &this->bot;

  this->top.setMaterial(&top_material);
  this->bot.setMaterial(&bot_material);
}

/* -------------------------------------------------------------------------- */
BimatInterface::~BimatInterface() {}

/* -------------------------------------------------------------------------- */
void BimatInterface::closingNormalGapForce(NodalField *close_force,
                                           bool predicting) {
  // top material information
  const Material & mat_t = this->top.getMaterial();
  double cp_t = mat_t.getCp();
  double cs_t = mat_t.getCs();
  double mu_t = mat_t.getShearModulus();
  double eta_t = cp_t / cs_t;
  double fact_t = this->time_step * cs_t / mu_t / eta_t;

  // bot material information
  const Material & mat_b = this->bot.getMaterial();
  double cp_b = mat_b.getCp();
  double cs_b = mat_b.getCs();
  double mu_b = mat_b.getShearModulus();
  double eta_b = cp_b / cs_b;
  double fact_b = this->time_step * cs_b / mu_b / eta_b;

  // accessors
  double * u_1_t = this->top.getDisp(1, predicting)->storage();
  double * u_1_b = this->bot.getDisp(1, predicting)->storage();
  double * f_1_t = this->top.getInternal(1)->storage();
  double * f_1_b = this->bot.getInternal(1)->storage();
  double * t0_1 = this->load[1]->storage();
  double * cf = close_force->storage();

  double fact_cf = fact_t + fact_b;

  fact_t /= fact_cf;
  fact_b /= fact_cf;

  for (int n=0; n<this->mesh.getNbNodes(); ++n) {
    double u_1_gap = u_1_t[n] - u_1_b[n];
    double du_1_t = fact_t * (t0_1[n] + f_1_t[n]);
    double du_1_b = fact_b * (t0_1[n] + f_1_b[n]);
    cf[n] = (u_1_gap/fact_cf + du_1_t + du_1_b);
  }
}

/* -------------------------------------------------------------------------- */
void BimatInterface::maintainShearGapForce(std::vector<NodalField *> &maintain_force) {
  // top material information
  const Material & mat_t = this->top.getMaterial();
  double cs_t = mat_t.getCs();
  double mu_t = mat_t.getShearModulus();
  double fact_t = cs_t / mu_t;

  // bot material information
  const Material & mat_b = this->bot.getMaterial();
  double cs_b = mat_b.getCs();
  double mu_b = mat_b.getShearModulus();
  double fact_b = cs_b / mu_b;

  double fact_mf = fact_t + fact_b;

  fact_t/=fact_mf;
  fact_b/=fact_mf;

  for (int d=0; d<this->mesh.getDim();d+=2) {
    // accessors
    double * f_t = this->top.getInternal(d)->storage();
    double * f_b = this->bot.getInternal(d)->storage();
    double * t0 = this->load[d]->storage();
    double * mf = maintain_force[d]->storage();

    for (int n=0; n<this->mesh.getNbNodes(); ++n) {
      //better for numerical error
      mf[n] = t0[n] + (fact_t * f_t[n] + fact_b * f_b[n]);
    }
  }
}

/* -------------------------------------------------------------------------- */
void BimatInterface::computeGap(std::vector<NodalField *> & gap,
                                bool predicting) {
  for (int d=0;d<this->mesh.getDim();++d) {
    double * top_disp = this->top.getDisp(d, predicting)->storage();
    double * bot_disp = this->bot.getDisp(d, predicting)->storage();
    double * gap_p = gap[d]->storage();

    for (int n=0; n<this->mesh.getNbNodes(); ++n) {
      gap_p[n] = top_disp[n] - bot_disp[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void BimatInterface::computeGapVelocity(std::vector<NodalField *> & gap_velo,
                                        bool predicting) {
  for (int d=0;d<this->mesh.getDim();++d) {
    double * top_velo = this->top.getVelo(d, predicting)->storage();
    double * bot_velo = this->bot.getVelo(d, predicting)->storage();
    double * gap_velo_p = gap_velo[d]->storage();

    for (int n = 0; n < this->mesh.getNbNodes(); ++n) {
      gap_velo_p[n] = top_velo[n] - bot_velo[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void BimatInterface::registerDumpField(const std::string &field_name) {

  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  if (world_rank != 0) return;

  int d = std::atoi(&field_name[field_name.length() - 1]);

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
  // field_name starts with "bot"
  else if (field_name.rfind("bot", 0) == 0) {
    // cut away "bot_" from field_name and give interface as dumper
    registered = this->bot.registerDumpFieldToDumper(field_name.substr(4),
						     field_name,
						     this);
  }

  if (!registered)
    Interface::registerDumpField(field_name);
}

__END_UGUCA__
