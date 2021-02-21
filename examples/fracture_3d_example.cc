/**
 * @file   fracture_3d_example.cc
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

#include "static_communicator_mpi.hh"
#include "uca_parameter_reader.hh"
#include "material.hh"

#include "defrig_interface.hh"
#include "barras_law.hh"

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
		<< " ./fracture_3d_example <input_file>" << std::endl;
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
  double x_length   = data.get<double>("x_length");
  double z_length   = data.get<double>("z_length");
  int nb_x_elements = data.get<int>("nb_x_elements");
  int nb_z_elements = data.get<int>("nb_z_elements");
  Mesh mesh(x_length, nb_x_elements,
	    z_length, nb_z_elements);

  // constitutive interface law
  BarrasLaw law(mesh,
		data.get<double>("tauc"),
		data.get<double>("dc"));

  // materials
  Material top_mat = Material(data.get<double>("E_top"),
			      data.get<double>("nu_top"),
			      data.get<double>("rho_top"));
  top_mat.readPrecomputedKernels();

  // interface
  DefRigInterface interface(mesh, top_mat, law);

  // external loading
  NodalField * ext_normal = interface.getNormalLoad();
  ext_normal->setAllValuesTo(data.get<double>("normal_load"));

  // time step
  double duration = data.get<double>("duration");
  double time_step = data.get<double>("tsf") * interface.getStableTimeStep();
  int nb_time_steps = duration/time_step;
  interface.setTimeStep(time_step);

  // initialization
  interface.init();

  // heterogeneity for nucleation: decreased strength
  NodalField * X = mesh.getCoords()[0];
  NodalField * Z = mesh.getCoords()[2];
  NodalField * tau_max = law.getTauMax();
  double a0 = data.get<double>("a0");
  for (int i=0;i<mesh.getNbNodes(); ++i)
    if ((std::abs((*X)(i) - x_length/2.) < a0/2.) &&
	(std::abs((*Z)(i) - z_length/2.) < a0/2.))
      (*tau_max)(i) = 0.;

  // dumping
  interface.initDump("fracture_3d_example",".");
  interface.registerDumpFields(data.get<std::string>("dump_fields"));
  interface.dump(0,0);
  unsigned int dump_int = std::max(1, nb_time_steps/data.get<int>("nb_dumps"));

  // time stepping
  for (int s=1; s<=nb_time_steps; ++s) {
    if (world_rank==0) {
      std::cout << "s=" << s << "/" << nb_time_steps << "\r";
      std::cout.flush();
    }

    // time integration
    interface.advanceTimeStep();

    // dump
    if (s % dump_int == 0)
      interface.dump(s,s*time_step);
  }

  delete comm;
  return 0;
}
