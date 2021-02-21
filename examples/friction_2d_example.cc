/**
 * @file   friction_2d_example.cc
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

#include "cmath" // for pi
#include <limits> // for max of double

#include "static_communicator_mpi.hh"
#include "uca_parameter_reader.hh"
#include "material.hh"

#include "unimat_shear_interface.hh"
#include "linear_coulomb_friction_law.hh"

using namespace uguca;

int main(int argc, char *argv[]) {

  // communicator for parallel simulation
  StaticCommunicatorMPI * comm = StaticCommunicatorMPI::getInstance();
  int world_rank = comm->whoAmI();

  // get input file name
  std::string fname;
  if(argc<2) {
    if (world_rank==0) {
      std::cerr << "Not enough arguments:"
		<< " ./friction_2d_example <input_file>" << std::endl;
    }
    return 1;
  }
  else {
    fname = argv[1];
  }

  // read input file
  ParameterReader data;
  data.readInputFile(fname);

  // mesh
  double length   = data.get<double>("length");
  int nb_elements = data.get<int>("nb_elements");
  Mesh mesh(length, nb_elements);

  // constitutive interface law
  LinearCoulombFrictionLaw law(mesh,
			       data.get<double>("mus"),
			       data.get<double>("muk"),
			       data.get<double>("dc"));

  // materials
  Material top_mat = Material(data.get<double>("E_top"),
			      data.get<double>("nu_top"),
			      data.get<double>("rho_top"));
  top_mat.readPrecomputedKernels();

  // interface
  UnimatShearInterface interface(mesh, top_mat, law);

  // time step
  double duration = data.get<double>("duration");
  double time_step = data.get<double>("tsf") * interface.getStableTimeStep();
  int nb_time_steps = duration/time_step;
  interface.setTimeStep(time_step);

  // heterogeneity in system: "random" strength
  NodalField * X = mesh.getCoords()[0];
  NodalField * mus = law.getMuS();
  double mus_min = std::numeric_limits<double>::max();
  for (int i=0;i<mesh.getNbNodes(); ++i) {
    double x = (*X)(i);
    (*mus)(i) = data.get<double>("mus")
      + data.get<double>("mus_ampl") * (1.0*std::sin(2*x*M_PI)
					+ 0.5*std::cos(12*x*M_PI)
					+ 0.7*std::sin(18*x*M_PI));
    mus_min = std::min(mus_min,(*mus)(i));
  }
  if (world_rank==0) std::cout << "min mu_s: " << mus_min << std::endl;

  // external loading
  double shear_load_rate = data.get<double>("shear_load_rate");
  double normal_load = data.get<double>("normal_load");
  NodalField * ext_shear = interface.getShearLoad();
  NodalField * ext_normal = interface.getNormalLoad();
  ext_normal->setAllValuesTo(normal_load);
  ext_shear->setAllValuesTo(mus_min*std::abs(normal_load));

  // initialization
  interface.init();

  // dumping
  interface.initDump("friction_2d_example",".");
  interface.registerDumpFields(data.get<std::string>("dump_fields"));
  interface.dump(0,0);
  unsigned int dump_int = std::max(1, nb_time_steps/data.get<int>("nb_dumps"));

  // time stepping
  for (int s=1; s<=nb_time_steps; ++s) {
    if (world_rank==0) {
      std::cout << "s=" << s << "/" << nb_time_steps << "\r";
      std::cout.flush();
    }

    ext_shear->setAllValuesTo(mus_min*std::abs(normal_load) + s*time_step*shear_load_rate);

    // time integration
    interface.advanceTimeStep();

    // dump
    if (s % dump_int == 0)
      interface.dump(s,s*time_step);
  }

  delete comm;
  return 0;
}
