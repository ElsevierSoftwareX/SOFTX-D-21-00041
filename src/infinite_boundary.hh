/**
 * @file   infinite_boundary.hh
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
#ifndef __INFINITE_BOUNDARY_H__
#define __INFINITE_BOUNDARY_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "half_space.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class InfiniteBoundary : public HalfSpace {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  // side factor top=1 bot=-1

  InfiniteBoundary(Mesh & mesh,
		   int side_factor,
		   Material & material);

  virtual ~InfiniteBoundary();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void advanceTimeStepDirichlet();
  void predictTimeStepDirichlet();
  void advanceTimeStepNeumann();
 
 

private:

  void gatherCostumMeshForwardFFT(bool predicting = false);
  void computeExternal();
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  NodalField * getExternal(int d) { return this->external[d]; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  std::vector< NodalField *> external;
  std::vector< NodalField *> scratch;
  int nb_pc = 0;
};

__END_UGUCA__

//#include "infinite_boundary_impl.cc"

#endif /* __INFINITE_BOUNDARY_H__ */
