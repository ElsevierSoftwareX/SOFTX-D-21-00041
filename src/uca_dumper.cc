/**
 * @file   uca_dumper.cc
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
#include "uca_dumper.hh"

#include <cstdio>
#include <iomanip>
#include <iostream>

#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
Dumper::Dumper(Mesh & mesh) :
  mesh(mesh) {

  // default name and path
  this->setBaseName("standard-bname");
  this->path = ".";

  this->time_file = NULL;
  this->field_file = NULL;

  this->initiated = false;
}

/* -------------------------------------------------------------------------- */
Dumper::~Dumper() {
  this->closeFiles(true);
}

void Dumper::closeFiles(bool release_memory) {
  FileToFieldMap::iterator it = this->files_and_fields.begin();
  FileToFieldMap::iterator end = this->files_and_fields.end();

  for (; it != end; ++it) {
    it->first->close();
    if (release_memory) delete it->first;
  }

  if (this->initiated) {
    this->time_file->close();
    if (release_memory) delete this->time_file;

    this->coord_file->close();
    if (release_memory) delete this->coord_file;

    this->field_file->close();
    if (release_memory) delete this->field_file;
  }
}

/* -------------------------------------------------------------------------- */
void Dumper::setBaseName(const std::string & bname) {

  this->base_name = bname;

  this->info_file_name  = this->base_name + ".info";
  this->time_file_name  = this->base_name + ".time";
  this->coord_file_name = this->base_name + ".coord";
  this->field_file_name = this->base_name + ".fields";
  this->folder_name     = this->base_name + "-DataFiles";
}

/* -------------------------------------------------------------------------- */
void Dumper::initDump(const std::string & bname,
		      const std::string & path,
                      const Format format) {
  this->initiated = true;

  this->dump_format = format;

  this->setBaseName(bname);
  this->path = path;

  // create folder for files (works only on linux)
  // read/write/search permission for owner and group
  // read/search permissions for others
  std::string full_path_to_folder =
      this->path + Dumper::directorySeparator() + this->folder_name;
  Dumper::createDirectory(full_path_to_folder);

  // info file
  std::string path_to_info_file =
      this->path + Dumper::directorySeparator() + this->info_file_name;
  std::ofstream info_file(path_to_info_file, std::ios::out);
  info_file << "field_description " << this->field_file_name << std::endl;
  info_file << "time_description " << this->time_file_name << std::endl;
  info_file << "coord_description " << this->coord_file_name << std::endl;
  info_file << "folder_name " << this->folder_name << std::endl;
  switch (this->dump_format) {
    case Format::ASCII: {
      info_file << "output_format ascii" << std::endl;
      break;
    }
    case Format::CSV: {
      info_file << "output_format csv" << std::endl;
      this->separator = ",";
      break;
    }
    case Format::Binary: {
      info_file << "output_format binary" << std::endl;
      break;
    }
    default:
      throw std::runtime_error("Unsupported output format.");
  }
  info_file.close();

  // time file
  std::string path_to_time_file =
      this->path + Dumper::directorySeparator() + this->time_file_name;

  this->time_file = new std::ofstream(path_to_time_file, std::ios::out);

  (*this->time_file) << std::scientific << std::setprecision(10);

  // coord file
  std::string path_to_coord_file =
      this->path + Dumper::directorySeparator() + this->coord_file_name;

  this->coord_file = new std::ofstream(path_to_coord_file, std::ios::out);

  (*this->coord_file) << std::scientific << std::setprecision(10);

  // field file
  std::string path_to_field_file =
      this->path + Dumper::directorySeparator() + this->field_file_name;

  this->field_file = new std::ofstream(path_to_field_file, std::ios::out);

  // write coords file
  this->setCoords();
}

/* -------------------------------------------------------------------------- */
void Dumper::registerForDump(const std::string & field_name,
			     const NodalField * nodal_field) {

  // name and path
  std::string file_extension = ".out";
  if (this->dump_format == Format::CSV) {
    file_extension = ".csv";
  }
  std::string file_name = field_name + file_extension;
  std::string path_to_file = this->path + Dumper::directorySeparator() +
                             this->folder_name + Dumper::directorySeparator() +
                             file_name;

  // open file
  std::ofstream * new_file = new std::ofstream();

  switch (this->dump_format) {
    case Format::ASCII:
    case Format::CSV: {
      new_file->open(path_to_file, std::ios::out);
      (*new_file) << std::scientific << std::setprecision(10);
      break;
    }
    case Format::Binary: {
      new_file->open(path_to_file, std::ios::out | std::ios::binary);  // open as binary file
      break;
    }
    default:
      throw std::runtime_error("Unsupported output format.");
  }

  // keep reference to file
  this->files_and_fields[new_file] = nodal_field;

  // put info into field file
  (*this->field_file) << field_name << " " << file_name << std::endl;
}

/* -------------------------------------------------------------------------- */
// coords has all coordinates of all points
// prints coords to file
void Dumper::setCoords() {

  if (!this->initiated) return;
  std::vector<NodalField *> coords = this->mesh.getCoords();

  for (int n=0; n<this->mesh.getNbNodes(); ++n) {
    for (int d=0; d<this->mesh.getDim(); ++d) {
      if (d != 0) {
	(*this->coord_file) << this->separator;
      }
      (*this->coord_file) << (*coords[d])(n);
    }
    (*this->coord_file) << std::endl;
  }

}

/* -------------------------------------------------------------------------- */
void Dumper::registerDumpFields(const std::string & field_names,
				char delim) {

  // separate string by delim
  std::vector<std::string> field_name_list;
  std::stringstream ss(field_names);
  std::string substring;
  while (std::getline(ss, substring, delim)) {
    field_name_list.push_back(substring);
  }

  // register dump fields
  for (std::vector<std::string>::iterator it = field_name_list.begin();
       it != field_name_list.end();
       ++it) {
    this->registerDumpField(*it);
  }
}

/* -------------------------------------------------------------------------- */
void Dumper::dump(unsigned int step, double time) {

  if (!this->initiated) return;

  FileToFieldMap::iterator it = this->files_and_fields.begin();
  FileToFieldMap::iterator end = this->files_and_fields.end();

  for (; it!=end; ++it) {
    this->dumpField(it->first, it->second);
  }

  (*this->time_file) << step << this->separator << time << std::endl;
}

/* -------------------------------------------------------------------------- */
void Dumper::dumpField(std::ofstream * dump_file,
		       const NodalField * nodal_field) {
  if (!this->initiated) return;

  switch (this->dump_format) {
    case Format::ASCII:
    case Format::CSV: {
      for (int n = 0; n < this->mesh.getNbNodes(); ++n) {
        if (n != 0) (*dump_file) << this->separator;
        (*dump_file) << nodal_field->at(n);
      }
      (*dump_file) << std::endl;
      break;
    }
    case Format::Binary: {
      float temp = 0.0;
      for (int n = 0; n < this->mesh.getNbNodes(); ++n) {
        temp = (float)(nodal_field->at(n));
        (*dump_file).write((char *)&temp, sizeof(float));
      }
      break;
    }
    default:
      throw std::runtime_error("Unsupported output format.");
  }
}

std::string Dumper::directorySeparator() {
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
  return "\\";
#else
  return "/";
#endif
}

void Dumper::createDirectory(std::string path) {
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
  _mkdir(path.c_str());
#else
  mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
}
__END_UGUCA__
