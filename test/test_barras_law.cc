/**
 * @file   test_barras_law.cc
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

#include "material.hh"
#include "barras_law.hh"
#include "defrig_interface.hh"

using namespace uguca;

int main(){

  std::cout << "start test: test_barras_law" << std::endl;

  // information for checks
  Mesh mesh(1., 2);
  double tauc = 4e6;
  double dc = 0.01;
  Material material(1e9, 0.25, 1000, false);

  std::cout << "check constructor" << std::endl;
  BarrasLaw law(mesh, tauc, dc);
  DefRigInterface interface(mesh, material, law);
  interface.setTimeStep(1e-8);
  std::cout << "constructor correct -> success" << std::endl;

  std::cout << "check data" << std::endl;
  NodalField * tmp = law.getTauMax();
  if (std::abs((*tmp)(0) - tauc) / tauc > 1e-5){
    std::cout << "wrong tau_max: " << (*tmp)(0) << std::endl;
    return 1; // failure
  }
  tmp = law.getDc();
  if (std::abs((*tmp)(0) - dc) / dc > 1e-5){
    std::cout << "wrong d_c: " << (*tmp)(0) << std::endl;
    return 1; // failure
  }
  std::cout << "data correct -> success" << std::endl;

  std::cout << "check computeCohesiveForces" << std::endl;

  // fill empty cohesion vector for testing
  std::vector<NodalField *> cohesion;
  NodalField coh0(mesh.getNbNodes());
  cohesion.push_back(&coh0);
  NodalField coh1(mesh.getNbNodes());
  cohesion.push_back(&coh1);

  // access to various properties needed to apply values
  NodalField * tau0 = interface.getShearLoad();
  NodalField * sig0 = interface.getNormalLoad();
  NodalField * u0 = interface.getTop().getDisp(0);
  NodalField * u1 = interface.getTop().getDisp(1);

  // check shear: tau0 < tauc & u=0
  double tau0v = 0.9*tauc;
  double sig0v = 0;
  tau0->setAllValuesTo(tau0v);
  sig0->setAllValuesTo(sig0v);
  double u0v = 0.;
  double u1v = 0.;
  u0->setAllValuesTo(u0v);
  u1->setAllValuesTo(u1v);
  double val = tau0v;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0(0) - val) / val > 1e-5) || (coh0(0) * tau0v < 0)) {
    std::cout << "shear failed (" << val << "): " << coh0(0) << std::endl;
    return 1; // failure
  }

  // check normal: 0 < sig0 < tauc & u=0
  tau0v = 0.;
  sig0v = 0.9*tauc;
  tau0->setAllValuesTo(tau0v);
  sig0->setAllValuesTo(sig0v);
  u0v = 0.;
  u1v = 0.;
  u0->setAllValuesTo(u0v);
  u1->setAllValuesTo(u1v);
  val = sig0v;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh1(0) - val) / val > 1e-5) || (coh1(0) * sig0v < 0)) {
    std::cout << "normal failed (" << val << "): " << coh1(0) << std::endl;
    return 1; // failure
  }

  // check contact: sig0 < 0 & u=0
  tau0v = 0.;
  sig0v = -1.5*tauc;
  tau0->setAllValuesTo(tau0v);
  sig0->setAllValuesTo(sig0v);
  u0v = 0.;
  u1v = 0.;
  u0->setAllValuesTo(u0v);
  u1->setAllValuesTo(u1v);
  val = sig0v;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh1(0) - val) / val > 1e-5) || (coh1(0) * sig0v < 0)) {
    std::cout << "contact failed (" << val << "): " << coh1(0) << std::endl;
    return 1; // failure
  }

  // check decohesion: tau0 > tauc & sig0 > tauc & u=0
  tau0v = 1.2*tauc;
  sig0v = 1.5*tauc;
  tau0->setAllValuesTo(tau0v);
  sig0->setAllValuesTo(sig0v);
  u0v = 0.;
  u1v = 0.;
  u0->setAllValuesTo(u0v);
  u1->setAllValuesTo(u1v);
  val = tauc;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0(0) - val) / val > 1e-5) || (coh0(0) * tau0v < 0) ||
      (std::abs(coh1(0) - val) / val > 1e-5) || (coh1(0) * sig0v < 0)) {
    std::cout << "decohesion failed (" << val << "): "
	      << coh0(0) << "," << coh1(0) << std::endl;
    return 1; // failure
  }

  // check weakening: tau0 > tauc & sig0 > tauc & u < dc
  tau0v = 1.2*tauc;
  sig0v = 1.5*tauc;
  tau0->setAllValuesTo(tau0v);
  sig0->setAllValuesTo(sig0v);
  u0v = 0.1*dc;
  u1v = 0.1*dc;
  u0->setAllValuesTo(u0v);
  u1->setAllValuesTo(u1v);
  val = tauc * (1 - std::sqrt(u0v*u0v+u1v*u1v)/dc);
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0(0) - val) / val > 1e-5) || (coh0(0) * tau0v < 0) ||
      (std::abs(coh1(0) - val) / val > 1e-5) || (coh1(0) * sig0v < 0)) {
    std::cout << "weakening failed (" << val << "): "
	      << coh0(0) << "," << coh1(0) << std::endl;
    return 1; // failure
  }

  // check open: tau0 > tauc & sig0 > tauc & u > dc
  tau0v = 1.2*tauc;
  sig0v = 1.5*tauc;
  tau0->setAllValuesTo(tau0v);
  sig0->setAllValuesTo(sig0v);
  u0v = 0.9*dc;
  u1v = 0.9*dc;
  u0->setAllValuesTo(u0v);
  u1->setAllValuesTo(u1v);
  val = 0.;
  law.computeCohesiveForces(cohesion, false);
  if ((std::abs(coh0(0) - val) / val > 1e-5) || (coh0(0) * tau0v < 0) ||
      (std::abs(coh1(0) - val) / val > 1e-5) || (coh1(0) * sig0v < 0)) {
    std::cout << "open failed (" << val << "): "
	      << coh0(0) << "," << coh1(0) << std::endl;
    return 1; // failure
  }

  std::cout << "computeCohesiveForces -> success" << std::endl;

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
