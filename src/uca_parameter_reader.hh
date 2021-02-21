/**
 * @file   uca_parameter_reader.hh
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
#ifndef __PARAMETER_READER_HH__
#define __PARAMETER_READER_HH__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"

// std
#include <set>
#include <map>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class ParameterReader {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ParameterReader();
  virtual ~ParameterReader() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// read input file
  void readInputFile(std::string file_name);

  /// write input file
  void writeInputFile(std::string file_name) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  ///
  template<typename T>
  T get(std::string key) const;

  template<typename T>
  bool has(std::string key) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// type of data available
  std::set<std::string> data_types;

  /// data
  std::map<std::string,std::string> string_data;
  std::map<std::string,int> int_data;
  std::map<std::string,unsigned int> uint_data;
  std::map<std::string,double> double_data;
  std::map<std::string,bool> bool_data;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

__END_UGUCA__

//#include "parameter_reader_inline_impl.cc"

#endif /* __PARAMETER_READER_HH__ */
