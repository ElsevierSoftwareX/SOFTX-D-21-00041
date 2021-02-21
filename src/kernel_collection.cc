/**
 * @file   kernel_collection.cc
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

#include "kernel_collection.hh"
#include "precomputed_kernel.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
KernelCollection::KernelCollection() :
  H00(NULL),
  H01(NULL),
  H11(NULL),
  H22(NULL) {}

/* -------------------------------------------------------------------------- */
KernelCollection::~KernelCollection() {
  delete this->H00;
  delete this->H01;
  delete this->H11;
  delete this->H22;
}

/* -------------------------------------------------------------------------- */
void KernelCollection::readPrecomputedKernels(double nu, bool pstress,
					      const std::string & path) {

  std::string nu_str = std::to_string((int)(100*nu));
  std::string ps_str = "";
  if (pstress)
    ps_str = "_pstress";

  this->H00 = new PrecomputedKernel(path+"nu0." + nu_str + ps_str + "_h00.txt");
  this->H01 = new PrecomputedKernel(path+"nu0." + nu_str + ps_str + "_h01.txt");
  this->H11 = new PrecomputedKernel(path+"nu0." + nu_str + ps_str + "_h11.txt");
  this->H22 = new PrecomputedKernel(path+"h22.txt");

}

__END_UGUCA__
