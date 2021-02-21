/**
 * @file   test_limited_history.cc
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
#include <iostream>

#include "limited_history.hh"

using namespace uguca;

int main(){

  std::cout << "start test: limited_history" << std::endl;

  unsigned int size = 4;
  LimitedHistory lh1(size);

  std::cout << "check size" << std::endl;
  if (lh1.getSize() != size) {
    std::cerr << "size incorrect: " << lh1.getSize() << std::endl;
    return 1; // failure
  }
  std::cout << "size correct -> success" << std::endl;

  double add1 = 2;
  for (int i=0; i<add1; ++i)
    lh1.addCurrentValue(i+1);

  std::cout << "check nb history" << std::endl;
  if (lh1.getNbHistoryPoints() != 2) {
    std::cerr << "nb history wrong: " << lh1.getNbHistoryPoints() << std::endl;
    return 1; // failure
  }
  std::cout << "nb history correct -> success" << std::endl;

  double add2 = 3;
  for (int i=0; i<add2; ++i)
    lh1.addCurrentValue(add1+i+1);

  std::cout << "check nb history" << std::endl;
  if (lh1.getNbHistoryPoints() != size) {
    std::cerr << "nb history wrong: " << lh1.getNbHistoryPoints() << std::endl;
    return 1; // failure
  }
  std::cout << "nb history correct -> success" << std::endl;

  std::cout << "check index position" << std::endl;
  if (lh1.getIndexNow() != 2*size-add1-add2) {
    std::cerr << "index position wrong: " << lh1.getIndexNow() << std::endl;
    return 1; //failure
  }
  std::cout << "index position correct -> success" << std::endl;

  std::cout << "check history" << std::endl;
  for (unsigned int i=0; i<lh1.getNbHistoryPoints(); ++i) {
    double val = add1+add2-i;
    if (lh1.at(i) != val) {
      std::cerr << "wrong history value (" << val << "): "
		<< lh1.at(i) << std::endl;
      return 1; // failure
    }
  }
  std::cout << "history correct -> success" << std::endl;

  std::cout << "check change current value" << std::endl;
  double val = 102;
  lh1.changeCurrentValue(val);
  if (lh1.at(0) != val) {
    std::cerr << "change current value wrong: " << lh1.at(0) << std::endl;
    return 1; // failure
  }
  std::cout << "change current value correct -> success" << std::endl;

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
