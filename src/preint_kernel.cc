/**
 * @file   preint_kernel.cc
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
#include "preint_kernel.hh"

__BEGIN_UGUCA__

// BLAS INTERFACE
#ifdef UCA_USE_BLAS
extern "C"{
double ddot_(const int    * __restrict__ N,
             const double * __restrict__ a,
             const int    * __restrict__ inca,
             const double * __restrict__ b,
             const int    * __restrict__ incb);
}

double blas_dot(const int & N, const double * __restrict__ a, const int inca, const double * __restrict__ b, const int incb) {
  return ddot_(&N, a, &inca, b, &incb);
}
#endif /* UCA_USE_BLAS */

/* -------------------------------------------------------------------------- */
PreintKernel::PreintKernel(const Kernel * kernel) :
  kernel(kernel) {

}

/* -------------------------------------------------------------------------- */
PreintKernel::~PreintKernel() {

}

/* -------------------------------------------------------------------------- */
void PreintKernel::preintegrate(double time_factor, double time_step) {

  // time for truncation of this mode
  // Tcut = tcut * q * cs  <-> tcut = Tcut / q / cs
  double trunc = this->kernel->getTruncation() / time_factor;

  // Reduce nb_integration_int by 1 because of i+1 when calling at() and another +1 inside at()
  unsigned int nb_integration_int = (unsigned int)(trunc / time_step) - 1;
  this->values.resize(nb_integration_int);

  // compute trapezoidal integral over time step
  double k_i = this->kernel->at(0.);
  for (unsigned int i=0; i<nb_integration_int; ++i) {
    double k_ii = this->kernel->at( (i+1)*time_step * time_factor); // k after i-th step
    this->values[i] = 0.5 * (k_i + k_ii) * time_step * time_factor;
    k_i = k_ii;
  }

}

/* -------------------------------------------------------------------------- */
void PreintKernel::multiplyBy(double factor) {

  unsigned int nb_integration_int = this->values.size();
  for (unsigned int i=0; i<nb_integration_int; ++i) {
    this->values[i] *= factor;
  }

}

/* -------------------------------------------------------------------------- */
std::complex<double> PreintKernel::convolve(const LimitedHistory * __restrict__ U_r,
					    const LimitedHistory * __restrict__ U_i) {

  unsigned int nb_U = U_r->getNbHistoryPoints();

  double real = 0.;
  double imag = 0.;

  unsigned int index_now = U_r->getIndexNow();
  unsigned int size = U_r->getSize();
  if (nb_U > size) {
    std::cerr << "try to access history value beyond existence" << std::endl;
    throw nb_U;
  }

  const double * __restrict__ Ur_p = U_r->getValues();
  const double * __restrict__ Ui_p = U_i->getValues();


#ifdef UCA_USE_BLAS
  int inca = 1, incb = 1;
  if (size - index_now >= nb_U) {
    real += blas_dot(nb_U, &values[0], inca, &Ur_p[index_now], incb);
    imag += blas_dot(nb_U, &values[0], inca, &Ui_p[index_now], incb);
  }
  else {
    unsigned int length = size - index_now;
    real += blas_dot(length, &values[0], inca, &Ur_p[index_now], incb);
    imag += blas_dot(length, &values[0], inca, &Ui_p[index_now], incb);
    unsigned int offset = size - index_now;
    real += blas_dot(nb_U - offset, &values[offset], inca, &Ur_p[0], incb);
    imag += blas_dot(nb_U - offset, &values[offset], inca, &Ui_p[0], incb);
  }
#else /* VECTORIZED CODE */
  if (size - index_now >= nb_U) {
    for (unsigned int i = 0; i < nb_U; ++i) {
      double K_cum = this->values[i];
      real += K_cum * Ur_p[index_now + i];
      imag += K_cum * Ui_p[index_now + i];
    }
  }
  else {
    for (unsigned int i = 0; i < size - index_now; ++i) {
      double K_cum = this->values[i];
      real += K_cum * Ur_p[index_now + i];
      imag += K_cum * Ui_p[index_now + i];
    }
    unsigned int offset = size - index_now;
    for (unsigned int i = offset; i < nb_U; ++i) {
      double K_cum = this->values[i];
      real += K_cum * Ur_p[i - offset];
      imag += K_cum * Ui_p[i - offset];
    }
  }
#endif /* UCA_USE_BLAS */

  return {real, imag};
}

__END_UGUCA__
