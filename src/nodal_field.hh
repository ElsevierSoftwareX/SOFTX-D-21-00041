/**
 * @file   nodal_field.hh
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
#ifndef __NODAL_FIELD_H__
#define __NODAL_FIELD_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include <vector>
/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__

class NodalField {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NodalField(){};

  NodalField(int nb_nodes);

  virtual ~NodalField();

private:
  // private copy constructor: NodalField cannot be copied (for now to debug)
  NodalField(NodalField & to_copy);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void zeros();

  void setAllValuesTo(double value);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get number of nodes
  int getNbNodes() const { return this->nb_nodes; };

  // access the value of node n (reading and writing)
  inline double & operator()(int node);

  // access the value of node n (only reading for const nodal field)
  inline double at(int node) const ;

  // access to storage
  inline double * storage() { return this->field; };
  inline const double * storage() const { return this->field; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // number of nodes
  int nb_nodes;

  // nodal field
  double * field;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline double & NodalField::operator()(int node) {
  return this->field[node];
}

inline double NodalField::at(int node) const {
  return this->field[node];
}

__END_UGUCA__

#endif /* __NODAL_FIELD_H__ */
