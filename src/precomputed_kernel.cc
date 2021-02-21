/**
 * @file   precomputed_kernel.cc
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
#include "precomputed_kernel.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
PrecomputedKernel::PrecomputedKernel(const std::string & fname) {
  this->readKernelFromFile(fname);
}

/* -------------------------------------------------------------------------- */
PrecomputedKernel::~PrecomputedKernel() {
}

/* -------------------------------------------------------------------------- */
double PrecomputedKernel::at(double T) const {

  unsigned int index = (unsigned int)(T/this->delta_t);

  if (index >= this->values.size()-1) {
    std::cerr << "Try to get value of PrecomputedKernel for T=" << T
	      << " index="<<index<<" size="<<this->values.size()
	      << " dt="<<this->delta_t
	      << ", which is beyond the PrecomputedKernel's range" << std::endl;
    throw T;
  }

  double dk = this->values[index+1] - this->values[index];
  double dT = T - index*this->delta_t;

  return this->values[index] +  dk * dT/this->delta_t;
}

/* -------------------------------------------------------------------------- */
void PrecomputedKernel::readKernelFromFile(const std::string & fname) {

  std::ifstream file(fname);

  if (!file.is_open()){
    std::stringstream e;
    e << "Cannot find file named: " << fname << std::endl;
    throw std::runtime_error(e.str());
  }

  std::string line;

  std::getline(file,line);
  this->delta_t = std::stod(line);

  while ( std::getline(file,line,',') ) {
    this->values.push_back(std::stod(line));
  }

  file.close();

  this->trunc = this->delta_t *(this->values.size()-1);

}

__END_UGUCA__
