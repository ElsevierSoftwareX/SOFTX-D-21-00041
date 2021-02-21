/**
 * @file   test_material.cc
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

using namespace uguca;

int main(){

  std::cout << "start test: test_material" << std::endl;

  // information about material as reference for checks below
  // --------------------------------------------------------------
  double E  = 2.;
  double nu = 0.25;
  double rho = 10.;

  // results
  double mu_pstrain = 0.8; // shear modulus
  double mu_pstress = 0.8; // shear modulus
  double Cp_pstrain = 0.489898; // p-wave speed
  double Cp_pstress = 0.461880; // p-wave speed
  double Cs_pstrain = 0.282843; // s-wave speed
  double Cs_pstress = 0.282843; // s-wave speed
  double Cr_pstrain = 0.260046; // R-wave speed
  double Cr_pstress = 0.257669; // R-wave speed
  // --------------------------------------------------------------

  // PLANE STRAIN
  std::cout << "check plane strain material values" << std::endl;
  Material pstrain(E, nu, rho, false);

  if ((pstrain.getShearModulus() - mu_pstrain) / mu_pstrain > 1e-5) {
    std::cout << "wrong shear modulus: "
	      << pstrain.getShearModulus() << std::endl;
    return 1; // failure
  }

  if ((pstrain.getCp() - Cp_pstrain) / Cp_pstrain > 1e-5) {
    std::cout << "wrong p-wave speed: "
	      << pstrain.getCp() << std::endl;
    return 1; // failure
  }

  if ((pstrain.getCs() - Cs_pstrain) / Cs_pstrain > 1e-5) {
    std::cout << "wrong s-wave speed: "
	      << pstrain.getCs() << std::endl;
    return 1; // failure
  }

  if ((pstrain.getCr() - Cr_pstrain) / Cr_pstrain > 1e-5) {
    std::cout << "wrong R-wave speed: "
	      << pstrain.getCr() << std::endl;
    return 1; // failure
  }
  std::cout << "plane strain material correct -> success" << std::endl;

  // PLANE STRESS
  std::cout << "check plane stress material values" << std::endl;
  Material pstress(E, nu, rho, true);

  if ((pstress.getShearModulus() - mu_pstress) / mu_pstress > 1e-5) {
    std::cout << "wrong shear modulus: "
	      << pstress.getShearModulus() << std::endl;
    return 1; // failure
  }

  if ((pstress.getCp() - Cp_pstress) / Cp_pstress > 1e-5) {
    std::cout << "wrong p-wave speed: "
	      << pstress.getCp() << std::endl;
    return 1; // failure
  }

  if ((pstress.getCs() - Cs_pstress) / Cs_pstress > 1e-5) {
    std::cout << "wrong s-wave speed: "
	      << pstress.getCs() << std::endl;
    return 1; // failure
  }

  if ((pstress.getCr() - Cr_pstress) / Cr_pstress > 1e-5) {
    std::cout << "wrong R-wave speed: "
	      << pstress.getCr() << std::endl;
    return 1; // failure
  }
  std::cout << "plane stress material correct -> success" << std::endl;

  // KERNELS
  std::cout << "check if kernels are read correctly" << std::endl;
  pstrain.readPrecomputedKernels();
  Kernel * kernel = pstrain.getH00();
  double exp = 1.232050807569e+00;
  if ((kernel->at(0) - exp) / exp > 1e-5) {
    std::cerr << "wrong H00 kernel: " << kernel->at(0) << std::endl;
    return 1; // failure
  }

  kernel = pstrain.getH01();
  exp = 4.641016151378e-13;
  if ((kernel->at(0) - exp) / exp > 1e-5) {
    std::cerr << "wrong H01 kernel: " << kernel->at(0) << std::endl;
    return 1; // failure
  }

  kernel = pstrain.getH11();
  exp = 4.019237886467e-01;
  if ((kernel->at(0) - exp) / exp > 1e-5) {
    std::cerr << "wrong H11 kernel: " << kernel->at(0) << std::endl;
    return 1; // failure
  }

  kernel = pstrain.getH22();
  exp = 5.000000000000e-01;
  if ((kernel->at(0) - exp) / exp > 1e-5) {
    std::cerr << "wrong H22 kernel: " << kernel->at(0) << std::endl;
    return 1; // failure
  }
  std::cout << "kernels correct -> success" << std::endl;

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
