/**
 * @file   barras_law.hh
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
#ifndef __BARRAS_LAW_H__
#define __BARRAS_LAW_H__
/* -------------------------------------------------------------------------- */
#include "interface_law.hh"

/*
   This is not exactly the same law than used in Barras et al. 2014.
   Within the cohesive zone, this law can close in normal direction,
   whereas Barras' law will only maintain the normal gap.
   Also, here no friction behind the cohesive zone is applied.

   strength = tau_max ( 1 - |delta| / delta_c )
 */

__BEGIN_UGUCA__

class BarrasLaw : public InterfaceLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  BarrasLaw(Mesh & mesh,
	    double tau_max_default, double delta_c_default);
  virtual ~BarrasLaw() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
 void computeCohesiveForces(std::vector<NodalField *> &cohesion,
                            bool predicting = false);
 virtual void registerDumpField(const std::string &field_name);

 /* ------------------------------------------------------------------------ */
 /* Accessors                                                                */
 /* ------------------------------------------------------------------------ */
public:
  NodalField * getTauMax() { return &(this->tau_max); };
  NodalField * getDc() { return &(this->delta_c); };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  NodalField tau_max;
  NodalField delta_c;
  NodalField gap_norm;
};

__END_UGUCA__

//#include "barras_law_impl.cc"

#endif /* __BARRAS_LAW_H__ */
