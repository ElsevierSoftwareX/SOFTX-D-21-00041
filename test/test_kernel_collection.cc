/**
 * @file   test_kernel_collection.cc
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
#include <algorithm> // for std::min
#include <iomanip> // for std::setprecision
#include <cmath> // for std::abs(double)

#include "kernel_collection.hh"
#include "precomputed_kernel.hh"

using namespace uguca;

int main(){

  std::cout << "start test: test_kernel_collection" << std::endl;

  // information about kernel as reference for checks below
  // --------------------------------------------------------------
  double nu = 0.25;
  bool pstress = false;

  std::string path = "../../kernels/";
  const char *fnames[4] = { "nu0.25_h00.txt",
			    "nu0.25_h01.txt",
			    "nu0.25_h11.txt",
			    "h22.txt" };

  // just for cout information
  const char *knames[4] = { "H00",
			    "H01",
			    "H11",
			    "H22" };
  // --------------------------------------------------------------

  // load data
  KernelCollection kc1;
  kc1.readPrecomputedKernels(nu, pstress);

  Kernel * kc1H[4] = { kc1.getH00(),
		       kc1.getH01(),
		       kc1.getH11(),
		       kc1.getH22() };


  // --------------------------------------
  // check if correct kernels was loaded
  for (int k=0; k<4; ++k) {
    std::cout << "check if correct " << knames[k]
	      << " kernel was loaded" << std::endl;

    // access kernel data
    PrecomputedKernel * kc1h = dynamic_cast<PrecomputedKernel *>(kc1H[k]);
    std::vector<double> & hval = kc1h->getValues();

    // load verified data
    PrecomputedKernel expH(path + fnames[k]);
    std::vector<double> & exp_val = expH.getValues();

    // check length of kernel
    std::cout << "check length" << std::endl;
    if (kc1h->getSize() != expH.getSize()) {
      std::cout << "failed" << std::endl;
      return 1; // failure
    }
    std::cout << "length correct -> success" << std::endl;

    // check first 10 values in the kernel
    std::cout << "check values" << std::endl;
    int nb_checks = std::min(10,(int)expH.getSize());
    for (int i=0; i<nb_checks; ++i) {
      double expected = exp_val[i];
      double found = hval[i];
      std::cout << std::setprecision(9)
		<< "expect: " << expected
		<< " found: " << found << std::endl;
      double rel_error = std::abs(found - expected) / expected;
      if (rel_error > 1e-8) {
	std::cout << "failed" << std::endl;
	return 1; // failure
      }
    }
    std::cout << "value are correct -> success" << std::endl;
  }

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
