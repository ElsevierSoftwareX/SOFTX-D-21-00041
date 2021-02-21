/**
 * @file   test_rate_and_state_law.cc
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
#include <math.h>

#include <iostream>

#include "unimat_shear_interface.hh"
#include "material.hh"
#include "rate_and_state_law.hh"

using namespace uguca;

int main(){
  std::cout << "start test: test_rate_and_state_law" << std::endl;

  // information for checks
  Mesh mesh(1., 2);
  double a = 0.008;
  double b = 0.012;
  double Dc = 0.02;
  double V0 = 1.0e-6;
  double f0 = 0.6;
  double sigma = 10e6;
  double tau = 6e6;
  double Vinit = 1.0e-6;
  double theta = Dc / V0 *
    std::exp((a * std::log(2 * std::sinh(tau / a / std::abs(sigma))) - f0 -
              a * std::log(Vinit / V0)) / b);
  Material material(1e9, 0.25, 1000, false);

  std::cout << "check constructor" << std::endl;
  RateAndStateLaw law(mesh, a, b, Dc, V0, f0, theta, sigma,
                      RateAndStateLaw::EvolutionLaw::AgingLaw, true, 0);
  UnimatShearInterface interface(mesh, material, law);
  interface.setTimeStep(1e-8);
  interface.getTop().getVelo(0)->setAllValuesTo(Vinit);
  interface.getLoad(0)->setAllValuesTo(tau);
  std::cout << "constructor correct -> success" << std::endl;

  std::cout << "check data" << std::endl;
  double tol = 1e-5;
  NodalField *tmp = law.getA();
  if (std::abs((*tmp)(0) - a) / a > tol) {
    std::cout << "wrong a: " << (*tmp)(0) << std::endl;
    return 1;  // failure
  }
  tmp = law.getB();
  if (std::abs((*tmp)(0) - b) / b > tol) {
    std::cout << "wrong b: " << (*tmp)(0) << std::endl;
    return 1;  // failure
  }
  tmp = law.getTheta();
  if (std::abs((*tmp)(0) - theta) / theta > tol) {
    std::cout << "wrong theta: " << (*tmp)(0) << std::endl;
    return 1;  // failure
  }
  tmp = law.getSigma();
  if (std::abs((*tmp)(0) - sigma) / sigma > tol) {
    std::cout << "wrong theta: " << (*tmp)(0) << std::endl;
    return 1;  // failure
  }
  std::cout << "data correct -> success" << std::endl;

  std::cout << "check computeCohesiveForces (steady state)" << std::endl;

  // fill empty cohesion vector for testing
  std::vector<NodalField *> cohesion;
  NodalField coh0(mesh.getNbNodes());
  cohesion.push_back(&coh0);
  NodalField coh1(mesh.getNbNodes());
  cohesion.push_back(&coh1);

  law.computeCohesiveForces(cohesion);
  if ((std::abs(coh0(0) - tau) / tau > tol) || (coh0(0) * tau < 0)) {
    std::cout << "stick failed (" << tau << "): " << coh0(0) << std::endl;
    return 1; // failure
  }

  std::cout << "computeCohesiveForces -> success" << std::endl;

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
