/**
 * @file   linear_shear_cohesive_law.hh
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
#ifndef __LINEAR_SHEAR_COHESIVE_LAW_H__
#define __LINEAR_SHEAR_COHESIVE_LAW_H__
/* -------------------------------------------------------------------------- */
#include "interface_law.hh"
/*
   Linear cohesive law in shear direction only.
   No interpenetration allowed
   but also no opening allowed.
   Thus: should only be used for pure mode II fracture

   Parameter:
   Gc - fracture energy
   tau_c peak shear strength
   tau_r residual shear strength
 */

__BEGIN_UGUCA__

class LinearShearCohesiveLaw : public InterfaceLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:


  LinearShearCohesiveLaw(Mesh & mesh,
			 double Gc_default,
			 double tau_c_default,
			 double tau_r_default = 0.);

  virtual ~LinearShearCohesiveLaw() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
 void computeCohesiveForces(std::vector<NodalField *> &cohesion,
                            bool predicting = false);

 // dumper function
 virtual void registerDumpField(const std::string &field_name);

 /* ------------------------------------------------------------------------ */
 /* Accessors                                                                */
 /* ------------------------------------------------------------------------ */
public:
  NodalField * getGc() { return &(this->G_c); };
  NodalField * getTauc() { return &(this->tau_c); };
  NodalField * getTaur() { return &(this->tau_r); };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  NodalField G_c;
  NodalField tau_c;
  NodalField tau_r;
};

__END_UGUCA__

//#include "linear_shear_cohesive_law_impl.cc"

#endif /* __LINEAR_SHEAR_COHESIVE_LAW_H__ */
