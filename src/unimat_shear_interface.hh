/**
 * @file   unimat_shear_interface.hh
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
#ifndef __UNIMAT_SHEAR_INTERFACE_H__
#define __UNIMAT_SHEAR_INTERFACE_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "half_space.hh"
#include "interface_law.hh"
#include "interface.hh"
#include <iostream>

/*
  assume antisymmetry - works if interface remains in contact

  u_0_top == -u_0_bot
  f_0_top == -f_0_bot

  u_1_top == u_1_bot -> no opening allowed
  f_1_top == f_1_bot -> convolution response is symmetric

  u_2_top == -u_2_bot
  f_2_top == -f_2_bot

  does not work for mode I crack

  for combined mode I + mode II use bimaterial interface

  verified using TPV3
*/

__BEGIN_UGUCA__

class UnimatShearInterface : public Interface {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  UnimatShearInterface(Mesh & mesh,
		       Material & top_material,
		       InterfaceLaw & law);

  virtual ~UnimatShearInterface();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // compute force needed to close normal gap
  virtual void closingNormalGapForce(NodalField * close_force,
				     bool predicting = false);

  // compute force needed to maintain current shear gap
  virtual void maintainShearGapForce(std::vector<NodalField *> &maintain_force);

  // compute gap in displacement
  virtual void computeGap(std::vector<NodalField *> &gap,
                          bool predicting = false);

  // compute gap relative velocity
  virtual void computeGapVelocity(std::vector<NodalField *> &gap_velo,
                                  bool predicting = false);

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  virtual HalfSpace & getTop() { return this->top; };
  virtual HalfSpace & getBot() {
#ifdef UCA_VERBOSE
    std::cout << "Warning: UnimatShearInterface::getBot() returns the same "
	      << "HalfSpace as UnimatShearInterface::getTop()" << std::endl;
#endif /* UCA_VERBOSE */
    return this->top;
  };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // half spaces
  HalfSpace top;
};

__END_UGUCA__

//#include "unimat_shear_interface_impl.cc"

#endif /* __UNIMAT_SHEAR_INTERFACE_H__ */
