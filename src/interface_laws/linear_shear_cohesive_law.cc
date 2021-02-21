/**
 * @file   linear_shear_cohesive_law.cc
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
#include "linear_shear_cohesive_law.hh"
#include "interface.hh"

#include <cmath>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */

LinearShearCohesiveLaw::LinearShearCohesiveLaw(Mesh & mesh,
					       double Gc_default,
					       double tau_c_default,
					       double tau_r_default) :
  InterfaceLaw(mesh),
  G_c(mesh.getNbNodes()),
  tau_c(mesh.getNbNodes()),
  tau_r(mesh.getNbNodes())
{
  this->G_c.setAllValuesTo(Gc_default);
  this->tau_c.setAllValuesTo(tau_c_default);
  this->tau_r.setAllValuesTo(tau_r_default);
}

/* -------------------------------------------------------------------------- */
void LinearShearCohesiveLaw::computeCohesiveForces(std::vector<NodalField *> & cohesion,
                                                   bool predicting) {

  int nb = this->mesh.getNbNodes();

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
  double * Gc = this->G_c.storage();
  double * tauc = this->tau_c.storage();
  double * taur = this->tau_r.storage();

  // to be filled
  NodalField alpha_field(nb);
  double * alpha = alpha_field.storage();

  // coh1 > 0 is a adhesive force
  // coh1 < 0 is a contact pressure
  for (int n = 0; n<nb; ++n) {

    // avoid penetration "at any cost"
    // apply no normal cohesive force
    coh1[n] = std::min(coh1[n], 0.);

    double slope = pow(tauc[n] - taur[n],2) / 2. / Gc[n];
    double strength = std::max(tauc[n] - shear_gap[n] * slope,
			       taur[n]);

    // maximal shear cohesive force given by strength.
    // keep orientation of shear force
    alpha[n] = std::min(1.,std::abs(strength / tau_shear[n]));
  }

  // only in shear direction
  this->interface->multiplyFieldByScalar(cohesion, alpha_field, true);
}

/* -------------------------------------------------------------------------- */
void LinearShearCohesiveLaw::registerDumpField(const std::string & field_name) {

  // G_c
  if (field_name == "G_c") {
    this->interface->registerForDump(field_name,
				     &(this->G_c));
  }

  // tau_c
  else if (field_name == "tau_c") {
    this->interface->registerForDump(field_name,
				     &(this->tau_c));
  }

  // tau_r
  else if (field_name == "tau_r") {
    this->interface->registerForDump(field_name,
				     &(this->tau_r));
  }

  // do not know this field
  else {
    InterfaceLaw::registerDumpField(field_name);
  }

}

__END_UGUCA__
