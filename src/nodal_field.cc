/**
 * @file   nodal_field.cc
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
#include "nodal_field.hh"
#include <cstring>
#include <cassert>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
NodalField::NodalField(int nb_nodes) :
  nb_nodes(nb_nodes),
  field(NULL) {

  this->field = new double[nb_nodes];
  this->zeros();
}

/* -------------------------------------------------------------------------- */
NodalField::~NodalField() {
  delete[] this->field;
}

/* -------------------------------------------------------------------------- */
void NodalField::zeros() {
  this->setAllValuesTo(0.);
}

/* -------------------------------------------------------------------------- */
void NodalField::setAllValuesTo(double value) {
  std::fill_n(this->field,this->nb_nodes,value);

}

__END_UGUCA__
