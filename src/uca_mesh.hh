/**
 * @file   uca_mesh.hh
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
#ifndef __MESH_H__
#define __MESH_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "nodal_field.hh"
#include <cmath>

#include <fftw3.h>

/* -------------------------------------------------------------------------- */

__BEGIN_UGUCA__

/* --------------------------------------------------------------------------
 * Mesh has information about both representation:
 *   - spatial
 *   - spectral
 */
class Mesh {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Mesh(double Lx, int Nx,
       double Lz, int Nz);

  Mesh(double Lx, int Nx,
       double Lz, int Nz,
       std::vector<NodalField *> & coords_costum);

  Mesh(double Lx, int Nx);

  Mesh(double Lx, int Nx,
       std::vector<NodalField *> & coords_costum);

  ~Mesh();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // for load balance in parallel convolution
  template <typename T> void sortFFT  (T * un_sorted, T * sorted, int root_rank=0);
  template <typename T> void unsortFFT(T * sorted,    T * un_sorted, int root_rank=0);

  void sortAndScatterFFTModes(int * Uglobal, int * Ulocal, int root_rank=0);
  void sortAndScatterFFTModes(int * U, int root_rank=0){sortAndScatterFFTModes(U,U,root_rank);};

  void sortAndScatterFFTModes(double * Uglobal, double * Ulocal, int root_rank=0); //used for distributing Modes
  void sortAndScatterFFTModes(double * U, int root_rank=0){sortAndScatterFFTModes(U,U,root_rank);};

  void sortAndScatterFFTModes(fftw_complex * Uglobal, fftw_complex * Ulocal, int root_rank=0);
  void sortAndScatterFFTModes(fftw_complex * U, int root_rank=0){sortAndScatterFFTModes(U,U,root_rank);};

  void gatherAndSortFFTModes (fftw_complex * Uglobal, fftw_complex * Ulocal, int root_rank=0);
  void gatherAndSortFFTModes (fftw_complex * U, int root_rank=0){gatherAndSortFFTModes(U,U,root_rank);};

  // for costum mesh
  template <typename T> void sortCostumNodes  (T * un_sorted, T * sorted, int root_rank=0);
  template <typename T> void unsortCostumNodes(T * sorted,    T * un_sorted, int root_rank=0);

  void sortAndScatterCostumNodes(int * Uglobal, int * Ulocal, int root_rank=0);
  void sortAndScatterCostumNodes(int * U, int root_rank=0){sortAndScatterCostumNodes(U,U,root_rank);};

  void sortAndScatterCostumNodes(double * Uglobal, double * Ulocal, int root_rank=0);//used for distributing Nodes
  void sortAndScatterCostumNodes(double * U, int root_rank=0){sortAndScatterCostumNodes(U,U,root_rank);};

  void gatherAndSortCostumNodes(int * Uglobal, int * Ulocal, int root_rank=0);
  void gatherAndSortCostumNodes(int * U, int root_rank=0){gatherAndSortCostumNodes(U,U,root_rank);};

  void gatherAndSortCostumNodes(double * Uglobal, double * Ulocal, int root_rank=0);
  void gatherAndSortCostumNodes(double * U, int root_rank=0){gatherAndSortCostumNodes(U,U,root_rank);};

private:
  // for load balance in parallel convolution
  template <typename T> void sortAndScatterFFTModes(T * Uglobal, T * Ulocal, T * buffer, int size, int root_rank=0);
  template <typename T> void gatherAndSortFFTModes (T * Uglobal, T * Ulocal, T * buffer, int size, int root_rank=0);
  // for costum mesh
  template <typename T> void sortAndScatterCostumNodes(T * Uglobal, T * Ulocal, T * buffer, int size, int root_rank=0);
  template <typename T> void gatherAndSortCostumNodes (T * Uglobal, T * Ulocal, T * buffer, int size, int root_rank=0);

protected:
  void initGlobal();
  void initParallel();

private:
  // for FFTW SERIAL - datastructure real: N0 x N1
  void initCoordsGlobal(std::vector<NodalField*> & coords_global);
  void initWaveNumbersGlobal(std::vector<NodalField*> & wave_numbers_global);

  // for parallel implementation of computeStressFourierCoeff()
  void assignFFTModes(std::vector<NodalField*> & wave_numbers_global);
  void printMaximumNbProc(std::vector<double> & work_per_mode);
  int getMaximumNbProc(std::vector<double> & work_per_mode);
  void computeWorkPerMode(std::vector<NodalField*> & wave_numbers_global,
			  std::vector<double> & work_per_mode);

  // for serial implementation of nodal operations
  void assignGlobalNodes(std::vector<NodalField*> & coords_global);

  // for parallel implementation of all nodal operations

  // costum local coordinates - to match a FEM code.
  void initCoordsLocalCostum(std::vector<NodalField*> &  coords );

  void initSortCostumNodesMap();

  void checkCostumCoords(std::vector<NodalField*> & coords_global);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // accessors
  const std::vector<NodalField *> & getCoords(){return this->coords;}
  const std::vector<NodalField *> & getWaveNumbers(){return this->wave_numbers;}

  int getDim(){return this->dim;}
  int getNbNodes(){return this->nb_nodes_local;}// for fftw_mpi includes padding

  int getGlobalNbNodes(){return this->nb_nodes_global;} // for fftw_mpi includes padding N0x2*(N1/2+1)
  int getGlobalNbNodesX(){return this->nb_nodes_x_global;}
  int getGlobalNbNodesZ(){return this->nb_nodes_z_global;}
  int getMaxNbNodesPp(){return this->max_nb_nodes_pp;}
  int getMaxNbCostumNodesPp(){return this->max_nb_costum_nodes_pp;}
  int getMaxNbFFTPp(){return this->max_nb_fft_pp;}

  double getDeltaX() {return this->length_x / this->nb_nodes_x_global;};
  double getDeltaZ() {return this->length_z / this->nb_nodes_z_global;};

  int getGlobalNbFFTX(){return this->nb_fft_x_global;}
  int getGlobalNbFFTZ(){return this->nb_fft_z_global;}
  int getNbFFT(){return this->nb_fft_local;}
  double getLengthX(){return this->length_x;}
  double getLengthZ(){return this->length_z;}

  bool isSpatialDomainParallel(){return this->spatial_domain_parallel;}
  bool isCostum(){return this->is_costum;}

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  // global information ( initialized in contructor )

  int dim;

  // length of domain / replication length
  double length_x;
  double length_z;

  // number of nodes
  int nb_nodes_x_global; // global
  int nb_nodes_z_global; // global

  // additional info /switches
  bool spatial_domain_parallel;
  bool is_costum=false;

  // additional global and local information (initialized in init)
  int nb_nodes_global;   // global
  int nb_nodes_local = 0;    //local
  int max_nb_nodes_pp = 0;   //local for fftw mpi datastrucutre
  int max_nb_costum_nodes_pp = 0; // for costum mesh

  // Fourier modes
  int nb_fft_x_global; // global
  int nb_fft_z_global; // global

  int nb_fft_global;   // global
  int nb_fft_local = 0;       //local
  int max_nb_fft_pp = 0;   //local

  std::vector<NodalField *> wave_numbers;  // local {k,-,m}
  std::vector<NodalField *> coords; // local {x,y,z}

  int * sort_fft_modes_map;  // for load balance sort modes
  int * sort_costum_nodes_map;  // for sorting nodes to processes

  fftw_complex * fftw_complex_buffer; //to sort modes
  double * double_buffer; //to sort nodes

};

__END_UGUCA__

//#include "mesh_impl.cc"

#endif /* __MESH_H__ */
