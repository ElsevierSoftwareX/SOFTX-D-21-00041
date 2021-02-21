/**
 * @file   test_nodal_field.cc
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

#include "nodal_field.hh"

using namespace uguca;

int main(){
  std::cout << "start test: nodal_field" << std::endl;

  std::cout << "check initialization and getSize" << std::endl;
  int sz = 4;
  NodalField nf1(sz);
  if (nf1.getNbNodes() != sz) {
    std::cerr << "wrong size" << std::endl;
    return 1; // failure
  }
  std::cout << "size correct -> success" << std::endl;

  std::cout << "check zeros" << std::endl;
  nf1.zeros();
  double * nf1p = nf1.storage();
  for (int i=0; i<nf1.getNbNodes(); ++i) {
    if (nf1p[i] != 0) {
      std::cerr << "should be zero: " << nf1p[i] << std::endl;
      return 1; // failure
    }
  }
  std::cout << "zeros correct -> success" << std::endl;

  std::cout << "check setAllValuesTo" << std::endl;
  double vl=3.5;
  nf1.setAllValuesTo(vl);
  for (int i=0; i<nf1.getNbNodes(); ++i) {
    if (nf1p[i] != vl) {
      std::cerr << "should be " << vl << ": " << nf1p[i] << std::endl;
      return 1; // failure
    }
  }
  std::cout << "setAllValuesTo correct -> success" << std::endl;

  // set non-uniform values to check
  double fct = 2.;
  for (int i=0; i<nf1.getNbNodes(); ++i)
    nf1p[i] = i*fct;

  std::cout << "check accessing data" << std::endl;
  for (int i=0; i<nf1.getNbNodes(); ++i) {
    if ((nf1(i) != i*fct) && (nf1.at(i) != i*fct)) {
      std::cerr << "data access failed" << std::endl;
      return 1; // failure
    }
  }
  std::cout << "data accessing correct -> success" << std::endl;

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
