/**
 * @file   limited_history.hh
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
#ifndef __LIMITED_HISTORY_H__
#define __LIMITED_HISTORY_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include <iostream>
#include <fstream>

__BEGIN_UGUCA__

class LimitedHistory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  LimitedHistory(unsigned int size);
  virtual ~LimitedHistory();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // at the current value of the history
  inline void addCurrentValue(double value);
  inline void changeCurrentValue(double value);

  // get history value at index with index=0 : now
  inline double at(unsigned int index) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  unsigned int getSize() const { return this->size; };
  unsigned int getNbHistoryPoints() const { return std::min(this->nb_history_points,
							    this->size); };
  unsigned int getIndexNow() const {return this->index_now; }
  double * getValues() const {return this->values; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // number of accumulated history points
  unsigned int nb_history_points;

  // number of history entries
  unsigned int size;

  // index pointing to the newest entry
  unsigned int index_now;

  // values
  double * values;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline void LimitedHistory::addCurrentValue(double value) {

  if (this->index_now == 0)
    this->index_now = this->size;

  this->index_now -= 1;

  this->values[this->index_now] = value;

  // increase the counter of history points
  this->nb_history_points = std::min(this->nb_history_points + 1,
				     this->size);
}

/* -------------------------------------------------------------------------- */
inline void LimitedHistory::changeCurrentValue(double value) {
  this->values[this->index_now] = value;
}

/* -------------------------------------------------------------------------- */
inline double LimitedHistory::at(unsigned int index) const {
  if (index >= this->size) {
    std::cerr << "try to access history value beyond existence" << std::endl;
    throw index;
  }

  unsigned int i = (this->index_now + index) % this->size;
  return this->values[i];
}

__END_UGUCA__

#endif /* __LIMITED_HISTORY_H__ */
