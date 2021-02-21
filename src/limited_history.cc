/**
 * @file   limited_history.cc
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
#include "limited_history.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
LimitedHistory::LimitedHistory(unsigned int size) :
  nb_history_points(0),
  size(size) {

  this->index_now = 0;
  this->values = new double[size];

  for (unsigned int i=0; i<this->size; ++i)
    this->values[i] = 0.;
}

/* -------------------------------------------------------------------------- */
LimitedHistory::~LimitedHistory() {
  delete[] this->values;
}

__END_UGUCA__
