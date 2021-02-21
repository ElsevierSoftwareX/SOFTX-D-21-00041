/**
 * @file   interface_law.hh
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
#ifndef __INTERFACE_LAW_H__
#define __INTERFACE_LAW_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "uca_mesh.hh"
#include <sstream>

#include "nodal_field.hh"

__BEGIN_UGUCA__

class Interface; // <--- don't know if this works --------------------------

class InterfaceLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  InterfaceLaw(Mesh & mesh) : mesh(mesh) {};
  virtual ~InterfaceLaw() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void computeCohesiveForces(std::vector<NodalField *> & cohesion,
                                     bool predicting = false) = 0;

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
 virtual void setInterface(Interface *interface) {
   this->interface = interface;
 };

 /* ------------------------------------------------------------------------ */
 /* Class Members                                                            */
 /* ------------------------------------------------------------------------ */
protected:
  Mesh & mesh;
  Interface * interface;
};

__END_UGUCA__

//#include "interface_law_impl.cc"

#endif /* __INTERFACE_LAW_H__ */
