/**
 * @file   linear_coulomb_friction_law.cc
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
#include "linear_coulomb_friction_law.hh"
#include "interface.hh"

#include <cmath>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
LinearCoulombFrictionLaw::LinearCoulombFrictionLaw(Mesh & mesh,
						   double mu_s_default,
						   double mu_k_default,
						   double d_c_default,
						   double char_reg_time) :
  InterfaceLaw(mesh),
  reg_contact_pressure(mesh.getNbNodes()),
  mu_s(mesh.getNbNodes()),
  mu_k(mesh.getNbNodes()),
  d_c(mesh.getNbNodes()),
  char_time(mesh.getNbNodes()),
  reg_cont_pres_tmp(mesh.getNbNodes())
{
  if (d_c_default < 1e-12) {
    std::cerr << "d_c cannot be zero, and it is currently: " << d_c_default << std::endl;
    throw;
  }

  this->initialized = false;
  this->reg_contact_pressure.setAllValuesTo(0. );
  this->mu_s.setAllValuesTo(mu_s_default);
  this->mu_k.setAllValuesTo(mu_k_default);
  this->d_c.setAllValuesTo(d_c_default);
  this->char_time.setAllValuesTo(char_reg_time);
}

/* -------------------------------------------------------------------------- */
void LinearCoulombFrictionLaw::computeCohesiveForces(std::vector<NodalField *> & cohesion,
						     bool predicting) {

  unsigned int nb = this->mesh.getNbNodes();

  // find forces needed to close normal gap
  this->interface->closingNormalGapForce(cohesion[1], predicting);
  double * coh1 = cohesion[1]->storage();

  // find force needed to maintain shear gap
  this->interface->maintainShearGapForce(cohesion);

  // get norm of shear cohesion
  NodalField shear_trac_norm(nb);
  this->interface->computeNorm(cohesion, shear_trac_norm, true);
  double * tau_shear = shear_trac_norm.storage();

  // find current gap
  std::vector<NodalField *> gap = this->interface->getBufferField();
  this->interface->computeGap(gap, predicting);

  // compute norm of shear gap
  NodalField shear_gap_norm(nb);
  this->interface->computeNorm(gap, shear_gap_norm, true);
  double * shear_gap = shear_gap_norm.storage();

  // interface properties
  double * mus = this->mu_s.storage();
  double * muk = this->mu_k.storage();
  double * dc  = this->d_c.storage();

  // initialize regularized contact pressure
  if (!this->initialized) {
    for (unsigned int n = 0; n<nb; ++n) {
      this->reg_contact_pressure(n) = coh1[n];
    }
    this->initialized = true;
  }

  // coh1 > 0 is a adhesive force
  // coh1 < 0 is a contact pressure
  for (unsigned int n = 0; n<nb; ++n) {
    // avoid penetration "at any cost"
    // apply no normal cohesive force
    coh1[n] = std::min(coh1[n], 0.);
  }

  // regularized contact pressure
  double * reg_sig = NULL;
  if (predicting) {
    for (unsigned int n = 0; n<nb; ++n)
      this->reg_cont_pres_tmp(n) = this->reg_contact_pressure(n);
    this->computeRegContactPressure(cohesion[1],
				    &(this->reg_cont_pres_tmp));
    reg_sig = this->reg_cont_pres_tmp.storage();
  }
  else {
    this->computeRegContactPressure(cohesion[1],
				    &(this->reg_contact_pressure));
    reg_sig = this->reg_contact_pressure.storage();
  }

  // to be filled
  NodalField alpha_field(nb);
  double * alpha = alpha_field.storage();

  for (unsigned int n = 0; n<nb; ++n) {

    // compute friction coefficient
    double dmu = (1 - shear_gap[n] / dc[n]) * (mus[n] - muk[n]);
    double mu = std::max(muk[n], muk[n] + dmu);

    // maximal shear cohesive force given by strength.
    // keep orientation of shear force
    double strength = std::abs(mu * reg_sig[n]);

    // maximal shear cohesive force given by strength.
    // keep orientation of shear force
    alpha[n] = std::min(1.,std::abs(strength / tau_shear[n]));
  }

  // only in shear direction
  this->interface->multiplyFieldByScalar(cohesion, alpha_field, true);
}

/* -------------------------------------------------------------------------- */
void LinearCoulombFrictionLaw::computeRegContactPressure(NodalField * cohesion_1,
							 NodalField * reg_cont_pres) {

  unsigned int nb = this->mesh.getNbNodes();
  double dt = this->interface->getTimeStep();

  double * reg_sig = reg_cont_pres->storage();
  double * coh1 = cohesion_1->storage();
  double * tc  = this->char_time.storage();

  for (unsigned int n = 0; n<nb; ++n) {

    // regularized
    if (tc[n] > 0) {
      // interface opened -> no history to preserve
      if (std::abs(coh1[n]) < 1e-12)
	reg_sig[n] = 0.;
      else
	reg_sig[n] = (reg_sig[n] + dt / tc[n] * coh1[n]) / (1 + dt / tc[n]);
    }
    // not regularized
    else {
      reg_sig[n] = coh1[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void LinearCoulombFrictionLaw::registerDumpField(const std::string & field_name) {

  // mu_s
  if (field_name == "mu_s") {
    this->interface->registerForDump(field_name,
				     &(this->mu_s));
  }

  // mu_k
  else if (field_name == "mu_k") {
    this->interface->registerForDump(field_name,
				     &(this->mu_k));
  }

  // d_c
  else if (field_name == "d_c") {
    this->interface->registerForDump(field_name,
				     &(this->d_c));
  }

  // reg_cont_pres
  else if (field_name == "reg_cont_pres") {
    this->interface->registerForDump(field_name,
				     &(this->reg_contact_pressure));
  }

  // do not know this field
  else {
    InterfaceLaw::registerDumpField(field_name);
  }

}

__END_UGUCA__

