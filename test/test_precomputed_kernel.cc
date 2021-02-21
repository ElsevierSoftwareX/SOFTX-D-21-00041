/**
 * @file   test_precomputed_kernel.cc
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
#include <iomanip> // for std::setprecision
#include <limits> // for std::numeric_limits<>
#include <stdexcept> // for std::runtime_error

#include "precomputed_kernel.hh"
#include <cmath>

using namespace uguca;

int main(){

  std::cout << "start test: test_precomputed_kernel" << std::endl;

  // for checking values
  double rel_error = std::numeric_limits<double>::max();
  double expected, found;

  // information about kernel as reference for checks below
  // --------------------------------------------------------------
  std::string kernel_fname = "test_precomputed_kernel.txt";
  unsigned int kernel_length = 2001;
  double kernel_dt = 0.05;
  int nb_kernel_value = 4;
  double kernel_values[] = { 1.23205203845265,
			     1.23085853340430,
			     1.22728287570183,
			     1.22133966063049 };
  // --------------------------------------------------------------

  // test that it fails if files does not exist
  std::cout << "check reading non-existing kernel file" << std::endl;
  bool caught_exception = true;
  try {
    PrecomputedKernel pk1("non-existing-kernel-file.txt");
    caught_exception = false;
  }
  catch (std::runtime_error &e) {
    std::cout << "caught exception -> success" << std::endl;
  }
  if (!caught_exception) {
    std::cout << "failed" << std::endl;
    return 1; // failure
  }

  // test reading file
  std::cout << "check reading kernel file" << std::endl;
  PrecomputedKernel pk2(kernel_fname);
  std::cout << "read without exception -> success" << std:: endl;

  // test content of read information: length
  std::cout << "check length of kernel" << std::endl;
  if (pk2.getSize() != kernel_length) {
    std::cout << "wrong length" << std::endl;
    return 1; // failure
  }
  std::cout << "kernel has correct length -> success" << std::endl;

  // test content of read information: dt
  std::cout << "check delta t of kernel" << std::endl;
  expected = kernel_dt;
  found = pk2.getDt();
  rel_error = std::abs(found - expected) / expected;
  std::cout << std::setprecision(6)
	    << "expected: " << expected
	    << " found: " << found << std::endl;
  if (rel_error > 1e-6) {
    std::cout << "failed" << std::endl;
    return 1; // failure
  }
  std::cout << "kernel has correct delta t -> success" << std::endl;

  // test content of read information: values
  std::cout << "check kernel values" << std::endl;
  std::vector<double> & values = pk2.getValues();
  for (int i=0; i<nb_kernel_value; i++) {
    rel_error = std::abs(values[i] - kernel_values[i]) / kernel_values[i];
    std::cout << std::setprecision(6)
	      << "expect: " << kernel_values[i]
	      << " found: " << values[i] << std::endl;
    if (rel_error > 1e-6) {
      std::cout << "failed" << std::endl;
      return 1; // failure
    }
  }
  std::cout << "kernel has correct values -> success" << std::endl;

  // test access by time
  std::cout << "check kernel access with .at(T)" << std::endl;
  int idx = 1;
  expected = kernel_values[idx];
  found = pk2.at(idx * kernel_dt);
  rel_error = std::abs(found - expected) / expected;
  if (rel_error > 1e-6) {
    std::cout << "failed" << std::endl;
    return 1; // failure
  }
  std::cout << "kernel accesses correctly -> success" << std::endl;

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}

