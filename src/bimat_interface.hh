/**
 * @file   bimat_interface.hh
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
#ifndef __BIMAT_INTERFACE_H__
#define __BIMAT_INTERFACE_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "half_space.hh"
#include "interface_law.hh"
#include "interface.hh"

__BEGIN_UGUCA__

class BimatInterface : public Interface {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  BimatInterface(Mesh & mesh,
		 Material & top_material,
		 Material & bot_material,
		 InterfaceLaw & law);

  virtual ~BimatInterface();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // compute force needed to close normal gap
  virtual void closingNormalGapForce(NodalField * close_force,
				     bool predicting = false);

  // compute force needed to maintain current shear gap
  virtual void maintainShearGapForce(std::vector<NodalField *> & maintain_force);

  // compute gap in displacement
  virtual void computeGap(std::vector<NodalField *> & gap,
			  bool predicting = false);

  // compute gap relative velocity
  virtual void computeGapVelocity(std::vector<NodalField *> & gap_velo,
				  bool predicting = false);

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  virtual HalfSpace & getTop() { return this->top; };
  virtual HalfSpace & getBot() { return this->bot; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  // half spaces
  HalfSpace top;
  HalfSpace bot;
};

__END_UGUCA__

//#include "bimat_interface_impl.cc"

#endif /* __BIMAT_INTERFACE_H__ */
