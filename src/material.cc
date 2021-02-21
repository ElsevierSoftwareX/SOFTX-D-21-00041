/**
 * @file   material.cc
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
#include <cmath>
#include "material.hh"

__BEGIN_UGUCA__

Material::Material(double E, double nu, double rho, bool pstress) :
  E(E), nu(nu), rho(rho), pstress(pstress) {

  this->computeShearModulus();
  this->computeFirstLame();
  this->computeCp();
  this->computeCs();
  this->computeCr();
}

void Material::computeShearModulus() {
  this->mu = 0.5*this->E / (1+this->nu);
}

void Material::computeFirstLame() {
  // 3D or plane-strain
  if (!this->pstress)
    this->lambda = this->nu*this->E / (1+this->nu) / (1-2*this->nu);
  else
    this->lambda = this->nu*this->E / (1 - this->nu*this->nu);
}

void Material::computeCp() {
  this->cp = std::sqrt((this->lambda + 2*this->mu)/this->rho);
}

void Material::computeCs() {
  this->cs = std::sqrt(this->mu/this->rho);
}

void Material::computeCr() {
  // Freund p .162 eq 4.3.8 : D = 0 for v = cR
  double cr_min = 0.20 * this->cs; // according to estimate: cR = 0.32*cs for nu = 0.5
  double cr_max = 1.00 * this->cs; // according to estimate: cR = 0.87*cs for nu = 0
  this->cr = this->computeCr_bisect(cr_min, cr_max);
}

double Material::computeCr_bisect(double a, double b) {
  if (this->computeCr_fcn(a) * this->computeCr_fcn(b) > 0) {
    throw "Failed to find Rayleigh wave speed. Check bounds.";
  }
  double xtol = 1e-12;
  double atol = 1e-12;
  unsigned maxiter = 100;
  unsigned iter = 0;
  double c = a;
  while ((b - a) >= xtol && iter <= maxiter) {
    c = (a + b) / 2;
    double fc = this->computeCr_fcn(c);
    if (std::abs(fc) <= atol) {
      return c;
    } else if (fc * this->computeCr_fcn(a) <= 0) {
      b = c;
    } else {
      a = c;
    }
    ++iter;
  }
  if (iter == maxiter) {
    throw "Failed to find Rayleigh wave speed. Bisection algorithm did not converge.";
  }
  return c;
}

double Material::computeCr_fcn(double v) {
  double alpha_s = std::sqrt(1.0 - v * v / (this->cs * this->cs));
  double alpha_p = std::sqrt(1.0 - v * v / (this->cp * this->cp));
  double temp = (1.0 + alpha_s * alpha_s);
  return 4.0 * alpha_p * alpha_s - temp * temp;
}


__END_UGUCA__
