/**
 * @file   test_preint_kernel.cc
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
#include <iomanip> // for std::setprecision
#include <limits> // for std::numeric_limits<>

#include "preint_kernel.hh"
#include "precomputed_kernel.hh"
#include "limited_history.hh"
#include <cmath>

using namespace uguca;


int main(){

  std::cout << "start test: test_preint_kernel" << std::endl;

  // for checking values
  double rel_error = std::numeric_limits<double>::max();
  double expected, found;

  // information about kernel as reference for checks below
  // --------------------------------------------------------------
  std::string kernel_fname = "test_precomputed_kernel.txt";
  unsigned int kernel_length = 2001;
  //double kernel_dt = 0.05;
  double time_factor = 0.1;
  double time_step   = 0.01;
  PrecomputedKernel pk(kernel_fname);

  PreintKernel pik(&pk);

  // --------------------------------------------------------------

  // test preintegrate -- applied trapezoidal rule
  std::cout << "check preintegrate" << std::endl;
  pik.preintegrate(time_factor,time_step);
  int i=kernel_length/15;
  expected = 0.5*(pk.at(i*time_step*time_factor)+pk.at((i+1)*time_step*time_factor))*time_step*time_factor;
  found = pik.getValues()[i];
  rel_error = std::abs(found - expected) / expected;
  if (rel_error > 1e-12) {
    std::cout << "failed" << std::endl;
    return 1; // failure
  }
  std::cout << "kernel preintegrate correctly -> success" << std::endl;

  // --------------------------------------------------------------

  // test getSize
  std::cout << "check getSize" << std::endl;
  if ((unsigned int)(pk.getTruncation() / time_factor / time_step)-1 != pik.getSize()) {
    std::cout << "failed" << std::endl;
    return 1; // failure
  }
  std::cout << "kernel getSize correctly -> success" << std::endl;

  // --------------------------------------------------------------

  // test multiplyBy
  std::cout << "check multiplyBy" << std::endl;
  double factor=2.0;
  expected = factor * pik.getValues()[i];
  pik.multiplyBy(factor);
  found = pik.getValues()[i];
  rel_error = std::abs(found - expected) / expected;
  if (rel_error > 1e-12) {
    std::cout << "failed" << std::endl;
    return 1; // failure
  }
  std::cout << "kernel multiply by correctly -> success" << std::endl;

  // --------------------------------------------------------------

  // test convolve
  std::cout << "check convolve" << std::endl;
  // test case where (lt.getSize()- lt.getIndexNow()>=lt.getNbHistoryPoints())
  LimitedHistory lt(pik.getSize());
  for (double val=0.25; val<2.0;val+=0.25){
    lt.addCurrentValue(val);
  }
  expected = 0.0;
  for (unsigned int n=0; n<lt.getNbHistoryPoints();n++){
    expected+=pik.getValues()[n]*lt.getValues()[lt.getIndexNow()+n];
  }
  found = (pik.convolve(&lt,&lt)).real();
  rel_error = std::abs(found - expected) / expected;
  if (rel_error > 1e-12) {
    std::cout << "failed" << std::endl;
    return 1; // failure
  }
  std::cout << "kernel convolve correctly -> success" << std::endl;

  // test case where (lt.getSize()- lt.getIndexNow()<lt.getNbHistoryPoints())
  for (unsigned int n=0; n<lt.getSize();n++){
    for (double val=0.25; val<.0;val+=0.25){
      lt.addCurrentValue(val);
    }
  }
  expected = 0.0;
  for (unsigned int n=0; n<lt.getSize()-lt.getIndexNow(); n++){
    expected+=pik.getValues()[n]*lt.getValues()[lt.getIndexNow()+n];
  }
  for (unsigned int n=lt.getSize()-lt.getIndexNow(); n<lt.getSize(); n++){
    expected+=pik.getValues()[n]*lt.getValues()[n-(lt.getSize()-lt.getIndexNow())];
  }
  found = (pik.convolve(&lt,&lt)).real();
  rel_error = std::abs(found - expected) / expected;
  if (rel_error > 1e-12) {
    std::cout << "failed" << std::endl;
    return 1; // failure
  }
  std::cout << "kernel convolve loop correctly -> success" << std::endl;


  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
