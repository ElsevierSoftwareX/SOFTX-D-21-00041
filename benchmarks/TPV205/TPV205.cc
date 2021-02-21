/**
 * @file   TPV205.cc
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
#include <math.h>
#include <sys/time.h>
#include <unistd.h>

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "bimat_interface.hh"
#include "interface.hh"
#include "linear_shear_cohesive_law.hh"
#include "material.hh"
#include "static_communicator_mpi.hh"
#include "unimat_shear_interface.hh"
#include "uca_mesh.hh"

using namespace uguca;

int main(int argc, char *argv[]) {
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  double length_x_rpt = 30e3;
  double length_z_rpt = 30e3;

  double domain_factor=2.0;

  double a0 = 3e3; // nucleation domain
  double duration = 12.;
  double dump_int = 0.01;

  int nb_nodes_x = (length_x_rpt*domain_factor)/100;
  double time_step_factor = 0.35;//1;//35;

  int s_dump = 0;
  int nb_time_steps = 0;

  int nb_pc=1;

  // Argument processing

  int c;

  bool is_unimat_interface=false;
  bool no_residual=false;

  extern char* optarg;
  while ((c = getopt(argc, argv, "hurN:T:t:s:f:p:")) != -1) {
    switch (c) {
    case 'h':
      fprintf(stderr,
	      "%s\n"
	      "\t-h: print this message\n"
	      "\t-N: number of elements power of 2 (%d)\n"
	      "\t-T: number of timesteps (%d)\n"
	      "\t-t: time between dumps (%f)\n"
	      "\t-s: factor of domain size (%f)\n"
	      "\t-f: time step factor (%g)\n"
	      "\t-p: number of predictor corrector steps >=0 (%d)\n"
	      "\t-u: use unimaterial interface (default bimaterial interface with same materials)\n"
	      "\t-r: remove residual friciton \n",

	      argv[0], nb_nodes_x, nb_time_steps, dump_int, domain_factor,time_step_factor,nb_pc);
      return -1;
    case 'N': nb_nodes_x     = atoi(optarg); break;
    case 'T': nb_time_steps  = atof(optarg); break;
    case 't': dump_int       = atof(optarg); break;
    case 's': domain_factor  = atof(optarg); break;
    case 'f': time_step_factor = atof(optarg); break;
    case 'p': nb_pc          = atoi(optarg); break;
    case 'u': is_unimat_interface=true; break;
    case 'r': no_residual=true; break;

    default:
      fprintf(stderr, "Unknown option (-%c)\n", c);
      return -1;
    }
  }

  double length_x = domain_factor*length_x_rpt;
  double length_z = domain_factor*length_z_rpt;

  int nb_nodes_z = nb_nodes_x/length_x*length_z;

  // loading and properties outside nucleation

  double normal_load = -120e6;

  double f_c = 0.677;
  double f_r = 0.525;
  double dc = 0.4;

  double tau_c = -f_c * normal_load;
  double tau_r = -f_r * normal_load;

  double shear_load  = 70e6;

  // loading and properties nucleation

  double nuc_shear_load = 81.6e6;

  double weak_patch_shear_load = 78.0e6;

  double strong_patch_shear_load = 62.0e6;

  double strength_barrier_center = 7.5e3;

  // take away tau_r for simplicity
  if (no_residual){
    shear_load -= tau_r;
    nuc_shear_load -= tau_r;
    weak_patch_shear_load   -= tau_r;
    strong_patch_shear_load -= tau_r;

    tau_c -= tau_r;
    tau_r = 0.0;
  }
  double Gc = 0.5*(tau_c-tau_r)*dc;

  // mesh
  Mesh mesh(length_x,nb_nodes_x,
	    length_z,nb_nodes_z);

  // constitutive interface law
  LinearShearCohesiveLaw law(mesh,
			     Gc,tau_c,tau_r);

  // materials

  // cp = sqrt((lambda + 2*mu)/rho)
  // cs = sqrt(mu/rho)
  // lambda = nu*E/(1+nu)/(1-2*nu)
  // mu = 0.5*E/(1+nu)

  double Cp=6000.0;
  double Cs=3464.0;
  double rho = 2670.0;
  double mu = Cs*Cs*rho;
  double lambda = Cp*Cp*rho - 2.0*mu;
  double nu = 0.5*(lambda / (lambda + mu));
  double E = mu*(3.0*lambda + 2.0*mu)/(lambda + mu);
  printf("E=%g\nnu=%g\n",E,nu);

  Material top_mat = Material(E,nu,rho);
  top_mat.readPrecomputedKernels();
  Material bot_mat = Material(E,nu,rho);
  bot_mat.readPrecomputedKernels();
  if ((std::abs((top_mat.getCp()-Cp))>1e-15*Cp) ||
      (std::abs((top_mat.getCs()-Cs))>1e-15*Cs))
    return -1;

  // ---------------------------------------------------------------------------

  // weak interface

  //time_step
  Interface * interface;
  if (is_unimat_interface)
    interface = new UnimatShearInterface(mesh, top_mat, law);
  else
    interface = new BimatInterface(mesh, top_mat, bot_mat, law);

  interface->initPredictorCorrector(nb_pc);

  // time step
  double time_step = time_step_factor * interface->getStableTimeStep();
  interface->setTimeStep(time_step);

  // external loading
  interface->getLoad(0)->setAllValuesTo(shear_load);
  interface->getLoad(1)->setAllValuesTo(0.);//normal_load);
  interface->getLoad(2)->setAllValuesTo(0.);

  interface->init();

  if (! nb_time_steps)
    nb_time_steps = duration/time_step;


  //--------------
  // nucleation
  NodalField * load_0 = interface->getLoad(0);
  NodalField * tauc = law.getTauc();
  NodalField * Gamma_c = law.getGc();
  double tol = 0.1*length_x/nb_nodes_x/2.0;

  std::vector<NodalField *>  coords = mesh.getCoords();
  for (int i=0; i<mesh.getNbNodes(); i++) {
    if (std::abs( (*coords[0])(i) - length_x/2.0) < a0/2.0+tol &&
	std::abs( (*coords[2])(i) - length_z/2.0 + length_z_rpt/4.0) < a0/2.0+tol) {
      (*load_0)(i) = 0.5*(nuc_shear_load+shear_load);
    }
  }

  for (int i=0; i<mesh.getNbNodes(); i++) {
    if (std::abs( (*coords[0])(i) - length_x/2.0) < a0/2.0-tol &&
	std::abs( (*coords[2])(i) - length_z/2.0 + length_z_rpt/4.0) < a0/2.0-tol) {
	(*load_0)(i) = nuc_shear_load;
    }
  }

  //--------------
  // infinite strength zone;

  for (int i=0; i<mesh.getNbNodes(); i++) {
    if (std::abs( (*coords[0])(i) - length_x/2.0)  > length_x_rpt/2.0-tol ||
	std::abs( (*coords[2])(i) - length_z/2.0) > length_z_rpt/2.0-tol) {
	(*Gamma_c)(i) = 1e24*Gc;
	(*tauc)(i) = 1e24*tau_c;
    }
  }

  //-------------
  // strengh barrier higher
  for (int i=0; i<mesh.getNbNodes(); i++) {
    if (std::abs( (*coords[0])(i) - (length_x/2.0+strength_barrier_center)) < a0/2.0+tol &&
	std::abs( (*coords[2])(i) - length_z/2.0 + length_z_rpt/4.0) < a0/2.0+tol)
      {
	(*load_0)(i) = 0.5*(strong_patch_shear_load+shear_load);
      }
  }
  for (int i=0; i<mesh.getNbNodes(); i++) {
    if (std::abs( (*coords[0])(i) - (length_x/2.0+strength_barrier_center)) < a0/2.0-tol &&
	std::abs( (*coords[2])(i) - length_z/2.0 + length_z_rpt/4.0) < a0/2.0-tol)
      {
	(*load_0)(i) = strong_patch_shear_load;
      }
  }

  //-------------
  // strengh barrier lower
  for (int i=0; i<mesh.getNbNodes(); i++) {
    if (std::abs( (*coords[0])(i) - (length_x/2.0-strength_barrier_center)) < a0/2.0+tol &&
	std::abs( (*coords[2])(i) - length_z/2.0 + length_z_rpt/4.0) < a0/2.0+tol)
      {
	(*load_0)(i) = 0.5*(weak_patch_shear_load+shear_load);
      }
  }

  for (int i=0; i<mesh.getNbNodes(); i++) {
    if (std::abs( (*coords[0])(i) - (length_x/2.0-strength_barrier_center)) < a0/2.0-tol &&
	std::abs( (*coords[2])(i) - length_z/2.0 + length_z_rpt/4.0) < a0/2.0-tol)
      {
	(*load_0)(i) = weak_patch_shear_load;
      }
  }

  //--------------

  // dumping
  std::string bname = "TPV205_N"+std::to_string(nb_nodes_x)
    +"_s"+std::to_string((int)domain_factor)+"."+std::to_string((int)((fmod(domain_factor,1))*100))
    +"_tf0."+std::to_string((int)(time_step_factor*100)/10)+std::to_string((int)(time_step_factor*100)%10)
    +"_pc"+std::to_string(nb_pc);

  if (is_unimat_interface)
    bname +="_unimat";
  if (no_residual)
    bname +="_nores";

  std::cout<<bname<<std::endl;

  interface->initDump(bname,".",Dumper::Format::Binary);

  interface->registerDumpField("cohesion_0");
  interface->registerDumpField("cohesion_1");
  interface->registerDumpField("cohesion_2");

  interface->registerDumpField("top_disp_0");
  interface->registerDumpField("top_disp_1");
  interface->registerDumpField("top_disp_2");
  interface->registerDumpField("bot_disp_0");
  interface->registerDumpField("bot_disp_1");
  interface->registerDumpField("bot_disp_2");

  interface->registerDumpField("top_velo_0");
  interface->registerDumpField("top_velo_1");
  interface->registerDumpField("top_velo_2");
  interface->registerDumpField("bot_velo_0");
  interface->registerDumpField("bot_velo_1");
  interface->registerDumpField("bot_velo_2");
  /*
  interface->registerDumpField("top_internal_0");
  interface->registerDumpField("top_internal_1");
  interface->registerDumpField("top_internal_2");
  interface->registerDumpField("bot_internal_0");
  interface->registerDumpField("bot_internal_1");
  interface->registerDumpField("bot_internal_2");

  interface->registerDumpField("load_0");
  interface->registerDumpField("load_1");
  interface->registerDumpField("load_2");
  */
  interface->registerDumpField("tau_c");

  interface->dump(0,0);
  if (!s_dump)
    s_dump = dump_int / time_step + 1;

  if (world_rank==0) {
    std::cout << "time step     = " << time_step << std::endl;
    std::cout << "nb time steps = " << nb_time_steps << std::endl;
    std::cout << "s_dump = " << s_dump << std::endl;
  }

  if (world_rank==0) std::cout << "dumping..."<< std::endl;

  // free surface

  HalfSpace * top = &interface->getTop();
  HalfSpace * bot = &interface->getBot();

  double * u0_top_p = top->getDisp(0)->storage();
  double * v0_top_p = top->getVelo(0)->storage();
  double * u2_top_p = top->getDisp(2)->storage();
  double * v2_top_p = top->getVelo(2)->storage();

  double * u0_bot_p = bot->getDisp(0)->storage();
  double * v0_bot_p = bot->getVelo(0)->storage();
  double * u2_bot_p = bot->getDisp(2)->storage();
  double * v2_bot_p = bot->getVelo(2)->storage();

  int ww=log(nb_time_steps)/log(10);
  if (world_rank==0) {
    std::cout<<"profiling:\n";
    std::cout<<std::setw(ww)<<"s"<<"/"
	     <<std::setw(ww)<<"smax"
	     <<std::setw(10)<<"free surf"
	     <<std::setw(10)<<"time step"
	     <<std::setw(10)<<"dumping"<<"\n";
  }
  // time stepping

  for (int s=1; s<=nb_time_steps; ++s) {
    struct timeval t0, t1, t2, t3;

    if (false){//world_rank==0 && s%10==0) {
      std::cout << "s=" << s << "/" << nb_time_steps << "\r";
      std::cout<<std::flush;
    }

    gettimeofday(&t0, NULL);

    // free surface
    int nb_nodes_x = mesh.getGlobalNbNodesX();
    int nb_nodes_z =  mesh.getGlobalNbNodesZ();

    for (int i = 0; i < nb_nodes_x; ++i) {
      for (int j = 1; j < nb_nodes_z / 2; ++j) {
	int ij = i * nb_nodes_z + j;
	int ijsym = i * nb_nodes_z + nb_nodes_z - j;

	u0_top_p[ijsym] =  u0_top_p[ij];
	v0_top_p[ijsym] =  v0_top_p[ij];
	u2_top_p[ijsym] = -u2_top_p[ij];
	v2_top_p[ijsym] = -v2_top_p[ij];

	if (!is_unimat_interface) {
	  u0_bot_p[ijsym] =  u0_bot_p[ij];
	  v0_bot_p[ijsym] =  v0_bot_p[ij];
	  u2_bot_p[ijsym] = -u2_bot_p[ij];
	  v2_bot_p[ijsym] = -v2_bot_p[ij];
	}
      }
    }
    gettimeofday(&t1, NULL);

    // time integration
    interface->advanceTimeStep();

    gettimeofday(&t2, NULL);

    if (world_rank==0) //dump
      if (s % s_dump == 0)
	interface->dump(s,s*time_step);

    gettimeofday(&t3, NULL);

    if (world_rank == 0 && s % s_dump == 0) {
      std::cout<<std::setw(ww)<<s<<"/"
	       <<std::setw(ww)<<nb_time_steps
	       <<std::setw(10)<< (t1.tv_sec-t0.tv_sec) + (t1.tv_usec-t0.tv_usec) *1e-6
	       <<std::setw(10)<< (t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec) *1e-6
	       <<std::setw(10)<< (t3.tv_sec-t2.tv_sec) + (t3.tv_usec-t2.tv_usec) *1e-6 << "\n";
      std::cout.flush();
    }
  }

  if (world_rank==0)
    std::cout << "Cleaning up..." << std::endl;

  delete interface;

  StaticCommunicatorMPI::getInstance()->finalize();
  return 0;
}
