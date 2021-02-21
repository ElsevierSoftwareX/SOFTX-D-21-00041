/**
 * @file   test_dumper.cc
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
#include <stdio.h>

#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "uca_dumper.hh"
#include "uca_mesh.hh"

using namespace uguca;

class TestDumper : public Dumper {
public:
 TestDumper(Mesh& mesh)
     : Dumper(mesh), field1(mesh.getNbNodes()), field2(mesh.getNbNodes()) {}
 virtual void registerDumpField(const std::string& field_name) {
   (void)field_name;  // unused parameter
   this->registerForDump("field1", &(this->field1));
   this->registerForDump("field2", &(this->field2));
  }
  void closeAllFiles() { this->closeFiles(false); }
  NodalField field1;
  NodalField field2;
};

void removeFile(std::string path) { remove(path.c_str()); }

void cleanupFiles(std::string bname, std::string path, std::string file_ext){
  std::string sep = Dumper::directorySeparator();
  removeFile(path + sep + bname + ".coord");
  removeFile(path + sep + bname + ".fields");
  removeFile(path + sep + bname + ".info");
  removeFile(path + sep + bname + ".time");
  removeFile(path + sep + bname + "-DataFiles/field1" + file_ext);
  removeFile(path + sep + bname + "-DataFiles/field2" + file_ext);
  removeFile(path + sep + bname + "-DataFiles");
}

std::string random_string(size_t length) {
  static std::string chars =
      "0123456789"
      "abcdefghijklmnopqrstuvwxyz"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  thread_local static std::mt19937 rg{std::random_device{}()};
  thread_local static std::uniform_int_distribution<std::string::size_type>
      pick(0, chars.length() - 2);
  std::string str;
  str.reserve(length);
  while (length--) str += chars[pick(rg)];
  return str;
}

int checkInfo(std::string bname, std::string path, std::string format) {
  std::ifstream file(path + Dumper::directorySeparator() + bname + ".info");
  std::string line;
  std::vector<std::string> lines;
  while (std::getline(file, line)) {
    lines.push_back(line);
  }
  file.close();
  if (lines.size() != 5) {
    std::cout << "wrong # of lines in .info dump file" << std::endl;
    return 1;
  }
  if (lines[0] != ("field_description " + bname + ".fields") ||
      lines[1] != ("time_description " + bname + ".time")||
      lines[2] != ("coord_description " + bname + ".coord") ||
      lines[3] != ("folder_name " + bname + "-DataFiles") ||
      lines[4] != ("output_format " + format)) {
    std::cout << "wrong content found in .info dump file" << std::endl;
    return 1;
  }
  std::cout << "*.info correct" << std::endl;
  return 0; // success
}

int checkCoords(Mesh& mesh, std::string bname, std::string path,
                Dumper::Format format) {
  std::ifstream file(path + Dumper::directorySeparator() + bname + ".coord");
  std::vector<std::vector<double>> coords;
  if (!file.is_open()) {
    std::cout << "cannot open *.coord dump file" << std::endl;
    return 1;  // failure
  }
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss;
    iss.str(line);
    std::vector<double> row;
    double value;
    while (!iss.eof()) {
      iss >> value;
      row.push_back(value);
      if (format == Dumper::Format::CSV) {
        char c;
        iss >> c;
        if (c != ',') {
          std::cout << "wrong separator in *.coord" << std::endl;
        }
      }
    }
    if (row.size() == 3) {
      coords.push_back(row);
    } else {
      std::cout << "wrong dimensions in *.coord" << std::endl;
      return 1;  // failure
    }
  }
  file.close();
  if (coords.size() != (size_t)mesh.getNbNodes()) {
    std::cout << "wrong # of nodes in *.coord" << std::endl;
    return 1;  // failure
  }
  const std::vector<NodalField*> coords_ref = mesh.getCoords();
  double tol = 1e-10;
  for (int i = 0; i < mesh.getNbNodes(); ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      if (std::abs(coords[i][j] - (*coords_ref[j])(i)) > tol) {
        std::cout << "discrepancy found in *.coord" << std::endl;
        return 1;  // failure
      }
    }
  }
  std::cout << "*.coord correct" << std::endl;
  return 0; // success
}

int checkTimes(std::string bname, std::string path, Dumper::Format format) {
  std::ifstream file(path + Dumper::directorySeparator() + bname + ".time");
  std::vector<std::vector<double>> times;
  if (!file.is_open()) {
    std::cout << "cannot open *.time dump file" << std::endl;
    return 1;  // failure
  }
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss;
    iss.str(line);
    std::vector<double> row;
    double value;
    while (!iss.eof()) {
      iss >> value;
      row.push_back(value);
      if (format == Dumper::Format::CSV) {
        char c;
        iss >> c;
        if (c != ',') {
          std::cout << "wrong separator in *.coord" << std::endl;
        }
      }
    }
    if (row.size() == 2) {
      times.push_back(row);
    } else {
      std::cout << "wrong dimensions in *.time" << std::endl;
      return 1;  // failure
    }
  }
  file.close();
  if (times.size() != 2) {
    std::cout << "wrong # of steps in *.time" << std::endl;
    return 1;  // failure
  }
  double tol = 1e-10;
  if (times[0][0] != 0 || std::abs(times[0][1] - 0.0) > tol ||
      times[1][0] != 1 || std::abs(times[1][1] - 0.1) > tol) {
    std::cout << "discrepancy found in *.time" << std::endl;
    return 1;  // failure
  }
  std::cout << "*.time correct" << std::endl;
  return 0;  // success
}

int checkFields(std::string bname, std::string path,
                std::string file_ext) {
  std::ifstream file(path + Dumper::directorySeparator() + bname + ".fields");
  std::string line;
  std::vector<std::string> lines;
  while (std::getline(file, line)) {
    lines.push_back(line);
  }
  file.close();
  if (lines.size() != 2) {
    std::cout << "wrong # of lines in .fields dump file" << std::endl;
    return 1;
  }
  if (lines[0].compare("field1 field1" + file_ext) != 0 ||
      lines[1].compare("field2 field2" + file_ext) != 0) {
    std::cout << "wrong content found in .fields dump file" << std::endl;
    return 1;
  }
  std::cout << "*.fields correct" << std::endl;
  return 0;  // success
}

int checkField(Mesh& mesh, std::string bname, std::string path,
               NodalField& field, std::string name, Dumper::Format format,
               std::string file_ext) {
  std::string sep = Dumper::directorySeparator();
  std::string file_name = name + file_ext;
  std::string file_path = path + sep + bname + "-DataFiles" + sep + file_name;
  switch (format) {
    case Dumper::Format::ASCII:
    case Dumper::Format::CSV: {
      // read data
      std::ifstream file(file_path);
      std::vector<std::vector<double>> data;
      if (!file.is_open()) {
        std::cout << "cannot open " << file_name << " dump file" << std::endl;
        return 1;  // failure
      }
      std::string line;
      while (std::getline(file, line)) {
        std::istringstream iss;
        iss.str(line);
        std::vector<double> row;
        double value;
        while (!iss.eof()) {
          iss >> value;
          row.push_back(value);
          if (format == Dumper::Format::CSV) {
            char c;
            iss >> c;
            if (c != ',') {
              std::cout << "wrong separator in *.coord" << std::endl;
              return 1;  // failure
            }
          }
        }
        if (row.size() == (size_t)mesh.getNbNodes()) {
          data.push_back(row);
        } else {
          std::cout << "wrong dimensions in \"" << file_name << "\" dump file"
                    << std::endl;
          return 1;  // failure
        }
      }
      file.close();

      // inspect data
      if (data.size() != 2) {
        std::cout << "wrong # of steps in " << file_name << std::endl;
        return 1;  // failure
      }
      for (size_t i = 0; i < (size_t)mesh.getNbNodes(); ++i) {
        if (data[0][i] != 0) {
          std::cout << "discrepancy found in " << file_name << " at step 0"
                    << std::endl;
          return 1;  // failure
        }
      }
      double tol = 1.0e-10;
      for (size_t i = 0; i < (size_t)mesh.getNbNodes(); ++i) {
        if (std::abs(data[1][i] - field((int)i)) > tol) {
          std::cout << "discrepancy found in " << file_name << " at step 1"
                    << std::endl;
          return 1;  // failure
        }
      }
      break;
    }
    case Dumper::Format::Binary: {
      // read data
      std::ifstream file(file_path, std::ios::in | std::ios::binary);
      if (!file.is_open()) {
        std::cout << "cannot open " << file_name << " dump file" << std::endl;
        return 1;  // failure
      }
      file.seekg(0, std::ios::end);
      const size_t num_elements = file.tellg() / sizeof(float);
      file.seekg(0, std::ios::beg);
      std::vector<float> data(num_elements);
      file.read(reinterpret_cast<char*>(&data[0]), num_elements * sizeof(float));
      file.close();

      // inspect data
      if (data.size() != 2 * (size_t)mesh.getNbNodes()) {
        std::cout << "wrong # of dumped values in " << file_name << std::endl;
        return 1;  // failure
      }
      float tol = 1e-7; // machine epsilon of float = 1.19e-7
      for (size_t i = 0; i < (size_t)mesh.getNbNodes(); ++i) {
        if (data[i] != 0) {
          std::cout << "discrepancy found in " << file_name << " at step 0"
                    << std::endl;
          return 1;  // failure
        }
        if (std::abs(data[i + (size_t)mesh.getNbNodes()] - field(i)) > tol) {
          std::cout << "discrepancy found in " << file_name << " at step 1"
                    << std::endl;
          return 1;  // failure
        }
      }
      break;
    }
    default:
      std::cout << "unsupported dump format appears in this test" << std::endl;
      return 1;  // failure
  }

  std::cout << file_name << " correct" << std::endl;
  return 0;
}

int testDumper(Mesh& mesh, std::string bname, std::string path,
               Dumper::Format format, std::string format_str,
               std::string file_ext) {
  std::string _bname = bname + "_" + format_str;
  std::cout << "initialize " << format_str << " dumper" << std::endl;
  TestDumper* dumper = new TestDumper(mesh);
  dumper->initDump(_bname, path, format);
  dumper->registerDumpField("");

  std::cout << "dump initial values" << std::endl;
  dumper->dump(0, 0);

  std::cout << "assign and dump random values at step 1 (t = 0.1)" << std::endl;
  std::uniform_real_distribution<double> unif(0, 1);
  std::default_random_engine re;
  for (size_t i = 0; i < (size_t)mesh.getNbNodes(); ++i) {
    dumper->field1(i) = unif(re);
    dumper->field2(i) = unif(re);
  }
  dumper->dump(1, 0.1);

  dumper->closeAllFiles(); // close files

  std::cout << "check dumped files" << std::endl;
  if (checkInfo(_bname, path, format_str) ||
      checkFields(_bname, path, file_ext) ||
      checkCoords(mesh, _bname, path, format) ||
      checkTimes(_bname, path, format) ||
      checkField(mesh, _bname, path, dumper->field1, "field1", format,
                 file_ext) ||
      checkField(mesh, _bname, path, dumper->field2, "field2", format,
                 file_ext)) {
    cleanupFiles(_bname, path, file_ext);
    return 1; // failure
  }

  std::cout << "remove testing files" << std::endl;
  delete dumper;
  cleanupFiles(_bname, path, file_ext);

  std::cout << format_str << " dumper correct -> success" << std::endl;
  return 0;  // success
}

int main(){
  std::cout << "start test: test_dumper" << std::endl;

  // information about the test grid as reference for checks below
  // --------------------------------------------------------------
  double Lx = 0.3;
  int Nx = 2;
  double Lz = 0.5;
  int Nz = 3;
  Mesh mesh(Lx, Nx, Lz, Nz);
  std::string bname = "test_dump";
  std::string path = ".";
  // --------------------------------------------------------------
  std::string rnd_str = "_" + random_string(10);
  if (testDumper(mesh, bname + rnd_str, path, Dumper::Format::ASCII, "ascii", ".out") ||
      testDumper(mesh, bname + rnd_str, path, Dumper::Format::CSV, "csv", ".csv") ||
      testDumper(mesh, bname + rnd_str, path, Dumper::Format::Binary, "binary", ".out"))
    return 1;  // failure

  std::cout << "all checks passed -> overall success" << std::endl;
  return 0; // success
}
