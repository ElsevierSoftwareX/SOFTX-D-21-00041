/**
 * @file   test_parameter_reader.cc
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

#include "uca_parameter_reader.hh"

using namespace uguca;

int main(){
  std::cout << "start test: parameter_reader" << std::endl;

  ParameterReader pr1;

  std::cout << "check reading file" << std::endl;
  pr1.readInputFile("test_parameter_reader.txt");
  std::cout << "reading works -> success" << std::endl;

  std::cout << "check writing file" << std::endl;
  pr1.writeInputFile("test_parameter_reader.out");
  std::cout << "writing works -> success" << std::endl;

  // check content in re-read file. thus verify indirectly writing
  ParameterReader pr2;
  pr2.readInputFile("test_parameter_reader.out");

  // check 'has' function
  std::cout << "test 'has' function" << std::endl;
  if (!pr2.has<std::string>("string_1")) {
    std::cerr << "should have been true." << std::endl;
    return 1; // failure
  }
  if (pr2.has<std::string>("string_2")) {
    std::cerr << "should have been false." << std::endl;
    return 1; // failure
  }
  std::cout << "'has' function works -> success" << std::endl;

  // check content
  std::cout << "check content" << std::endl;
  if (pr2.get<std::string>("string_1") != "burger") {
    std::cerr << "wrong string" << std::endl;
    return 1; // failure
  }
  if (!pr2.get<bool>("bool_1")) {
    std::cerr << "wrong bool" << std::endl;
    return 1; // failure
  }
  if (pr2.get<double>("double_1") != 2.2) {
    std::cerr << "wrong double" << std::endl;
    return 1; // failure
  }
  if (pr2.get<int>("int_1") != 3) {
    std::cerr << "wrong int" << std::endl;
    return 1; // failure
  }
  if (pr2.get<unsigned int>("uint_1") != 4) {
    std::cerr << "wrong uint" << std::endl;
    return 1; // failure
  }
  std::cout << "content correct -> success" << std::endl;

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
