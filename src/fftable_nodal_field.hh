/**
 * @file   fftable_nodal_field.hh
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
#ifndef __FFTABLE_NODAL_FIELD_H__
#define __FFTABLE_NODAL_FIELD_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "uca_mesh.hh"

#include <fftw3.h>

#include "nodal_field.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__
/* **
 * MPI:
 * use default FFTW MPI datastructure N0x(N1/2+1)*2
 * note integer division rounds down
 *
 * Serial:
 * use default FFTW datastructure N0xN1
 */
class FFTableNodalField : public NodalField {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  FFTableNodalField(Mesh & mesh);

  virtual ~FFTableNodalField();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void forwardFFT(unsigned int root=0);
  void backwardFFT(unsigned int root=0);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  // get one value of frequency domain
  inline fftw_complex & fd(int f);

  // get access directly to frequency domain
  // WARNING: convert it to double (assuming that fftw_complex is double[2])
  inline fftw_complex * fd_storage() { return this->freq_dom_field; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  Mesh & mesh;

  // values in frequency domain in complex form
  fftw_complex * freq_dom_field;

  // the forward and backward plan
  fftw_plan forward_plan;
  fftw_plan backward_plan;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline fftw_complex & FFTableNodalField::fd(int f) {
  return this->freq_dom_field[f];
}


__END_UGUCA__

#endif /* __FFTABLE_NODAL_FIELD_H__ */
