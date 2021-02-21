/**
 * @file   fftable_nodal_field.cc
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
#include "fftable_nodal_field.hh"
#include "static_communicator_mpi.hh"
#include <cstring>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
FFTableNodalField::FFTableNodalField(Mesh & mesh) :
  /* use default FFTW datastrucutre N0xN1
   */
  NodalField(mesh.getGlobalNbNodesX()*mesh.getGlobalNbNodesZ()),
  mesh(mesh)
{

  if (this->mesh.getDim()==3) {

    /* use default FFTW datastrucutre N0xN1
     */
    int nb_fft_x = this->mesh.getGlobalNbFFTX();
    int nb_fft_z = this->mesh.getGlobalNbFFTZ();

    this->freq_dom_field = new fftw_complex[nb_fft_x*nb_fft_z];
    // init freq dom field with zeros
    memset(this->freq_dom_field,0.,nb_fft_x*nb_fft_z*sizeof(fftw_complex));

    this->forward_plan = fftw_plan_dft_r2c_2d(this->mesh.getGlobalNbNodesX(),
					      this->mesh.getGlobalNbNodesZ(),
					      this->field,
					      this->freq_dom_field,
					      FFTW_MEASURE);

    this->backward_plan = fftw_plan_dft_c2r_2d(this->mesh.getGlobalNbNodesX(),
					       this->mesh.getGlobalNbNodesZ(),
					       this->freq_dom_field,
					       this->field,
					       FFTW_MEASURE);
  }
  else {//2D is only serial
    int nb_fft_x = this->mesh.getGlobalNbFFTX();

    this->freq_dom_field = new fftw_complex[nb_fft_x];
    memset(this->freq_dom_field, 0., nb_fft_x*sizeof(fftw_complex));

    this->forward_plan = fftw_plan_dft_r2c_1d(this->mesh.getGlobalNbNodesX(),
					      this->field,
					      this->freq_dom_field,
					      FFTW_MEASURE);

    this->backward_plan = fftw_plan_dft_c2r_1d(this->mesh.getGlobalNbNodesX(),
					       this->freq_dom_field,
					       this->field,
					       FFTW_MEASURE);
  }
}

/* -------------------------------------------------------------------------- */
FFTableNodalField::~FFTableNodalField() {

  delete[] this->freq_dom_field;

  fftw_destroy_plan(this->forward_plan);
  fftw_destroy_plan(this->backward_plan);
  //fftw_cleanup(); // not needed because we don't know if all plans are destroyed
}

/* -------------------------------------------------------------------------- */
void FFTableNodalField::forwardFFT(unsigned int root) {

  if (StaticCommunicatorMPI::getInstance()->whoAmI()!=root) return;

  fftw_execute(this->forward_plan);

}
/* -------------------------------------------------------------------------- */
void FFTableNodalField::backwardFFT(unsigned int root) {

  if (StaticCommunicatorMPI::getInstance()->whoAmI()!=root) return;

  fftw_execute(this->backward_plan);

  double nb_nodes_global = this->mesh.getGlobalNbNodesX()*this->mesh.getGlobalNbNodesZ();
  for (int i=0; i<this->nb_nodes; ++i) {//for fftw_mpi it includes padding
    this->field[i] /= nb_nodes_global;
  }
}

__END_UGUCA__
