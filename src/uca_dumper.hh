/**
 * @file   uca_dumper.hh
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
#ifndef __DUMPER_H__
#define __DUMPER_H__
/* -------------------------------------------------------------------------- */

#include <fstream>
#include <map>
#include <sstream>
#include <string>

#include "nodal_field.hh"
#include "uca_common.hh"
#include "uca_mesh.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class Dumper {
public:
  enum class Format { ASCII, CSV, Binary };
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  typedef std::map<std::ofstream *, const NodalField *> FileToFieldMap;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Dumper(Mesh & mesh);
  virtual ~Dumper();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void initDump(const std::string &bname, const std::string &path,
                        const Format format = Format::ASCII);

  virtual void registerDumpField(const std::string &field_name) = 0;

  void registerDumpFields(const std::string & field_names,
			  char delim = ',');

  void registerForDump(const std::string &field_name,
                       const NodalField *nodal_field);

  void dump(unsigned int step, double time);

 protected:
  void setBaseName(const std::string & bname);

  void setCoords();

  void dumpField(std::ofstream * dump_file,
		 const NodalField * nodal_field);

protected:
  void closeFiles(bool release_memory);

  /* ------------------------------------------------------------------------ */
  /* File system related methods                                              */
  /* ------------------------------------------------------------------------ */
public:
  static std::string directorySeparator();
  static void createDirectory(std::string path);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

 /* ------------------------------------------------------------------------ */
 /* Class Members                                                            */
 /* ------------------------------------------------------------------------ */
protected:
  Mesh & mesh;
  Format dump_format;

private:
  // base name
  std::string base_name;

  // path to dumped files
  std::string path;

  // has dump been initiated?
  bool initiated;

  // information based on base_name
  std::string info_file_name;
  std::string time_file_name;
  std::string coord_file_name;
  std::string field_file_name;
  std::string folder_name;

  // files corresponding to field
  FileToFieldMap files_and_fields;

  // file with time stamps
  std::ofstream * time_file;

  // file with coord
  std::ofstream * coord_file;

  // file with field infos
  std::ofstream * field_file;

  // characteristics of dumper
  std::string separator = " ";

};

__END_UGUCA__

//#include "uca_dumper_impl.cc"

#endif /* __DUMPER_H__ */
