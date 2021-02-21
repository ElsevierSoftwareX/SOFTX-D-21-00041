/**
 * @file   precomputed_kernel.hh
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
/* -------------------------------------------------------------------------- */
#ifndef __PRECOMPUTED_KERNEL_H__
#define __PRECOMPUTED_KERNEL_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include <fstream>
#include <iostream>
#include <sstream>

#include <vector>

#include "kernel.hh"

__BEGIN_UGUCA__

class PrecomputedKernel : public Kernel {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PrecomputedKernel(const std::string & fname);
  virtual ~PrecomputedKernel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // read kernel from file generated before with fortran script
  void readKernelFromFile(const std::string & fname);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get value of Kernel at given T (by interpolation)
  virtual double at(double T) const;

  // get T for truncation
  virtual double getTruncation() const { return this->trunc; };

  // get number of nodes
  unsigned int getSize() const { return this->values.size(); };

  // get delta t
  double getDt() const { return this->delta_t; };

  // get direct access to values (only used to testing)
  std::vector<double> & getValues() { return this->values; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // nodal field
  std::vector<double> values;

  // truncation of kernel
  double trunc;

  // delta time of kernel entries
  double delta_t;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

__END_UGUCA__

#endif /* __PRECOMPUTED_KERNEL_H__ */
