/**
 * @file   kernel_collection.hh
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
#ifndef __KERNEL_COLLECTION_H__
#define __KERNEL_COLLECTION_H__
/* -------------------------------------------------------------------------- */
#include <string>

#include "uca_common.hh"
#include "kernel.hh"

__BEGIN_UGUCA__

class KernelCollection {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  KernelCollection();

  virtual ~KernelCollection();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // read precomputed Kernels from files
  // default path is replaced by cmake
  void readPrecomputedKernels(double nu, bool pstress = false,
			      const std::string & path = global_kernel_path);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  Kernel * getH00() { return this->H00; }
  Kernel * getH01() { return this->H01; }
  Kernel * getH11() { return this->H11; }
  Kernel * getH22() { return this->H22; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  Kernel * H00;
  Kernel * H01;
  Kernel * H11;
  Kernel * H22;

};

__END_UGUCA__

//#include "kernel_collection_impl.cc"

#endif /* __KERNEL_COLLECTION_H__ */
