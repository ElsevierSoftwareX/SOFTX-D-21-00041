/**
 * @file   uca_mesh.cc
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
#include "uca_mesh.hh"
#include "static_communicator_mpi.hh"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <numeric>      // std::iota
#include <algorithm>
#include <unistd.h> //sleep

__BEGIN_UGUCA__

/* --------------------------------------------------------------------------
 * 2D
 *
 * All operations in spatial domain are performed serially:
 * local==global for spatial domain
 * operations in fourier domain can be performed in parallel by default
 */
Mesh::Mesh(double Lx, int Nx) :
  dim(2),
  length_x(Lx),
  length_z(0.0),
  nb_nodes_x_global(Nx),
  nb_nodes_z_global(1),
  spatial_domain_parallel(false)
{
  // init global variables // allocate sorting maps
  this->initGlobal();

  std::vector<NodalField *> coords_global_tmp;
  std::vector<NodalField *> wavenb_global_tmp;

  coords_global_tmp.resize(this->dim);
  wavenb_global_tmp.resize(this->dim);

  for (int d=0; d<this->dim; ++d) {
    coords_global_tmp[d] = new NodalField(this->nb_nodes_global);
    wavenb_global_tmp[d] = new NodalField(this->nb_nodes_global);
  }
 
  // init coords default and wave number
  this->initCoordsGlobal(coords_global_tmp);
  this->initWaveNumbersGlobal(wavenb_global_tmp);

  // for parallel implementation
  // sort and assing modes for load balance
  this->assignFFTModes(wavenb_global_tmp);
  // spatial domain serial -> local=global
  this->assignGlobalNodes(coords_global_tmp);

  this->initParallel();

  for (int d=0; d<this->dim; ++d) {
    delete coords_global_tmp[d];
    delete wavenb_global_tmp[d];
  }
}
/* --------------------------------------------------------------------------
 * 2D costum
 *
 * All operations in spatial domain are performed in parallel
 * operations in fourier domain can be performed in parallel by default
 */
Mesh::Mesh(double Lx, int Nx,
	   std::vector<NodalField *> & coords_costum) :
  dim(2),
  length_x(Lx),
  length_z(0.0),
  nb_nodes_x_global(Nx),
  nb_nodes_z_global(1),
  spatial_domain_parallel(true),
  is_costum(true)
{

  // init global variables // allocate sorting maps
  this->initGlobal();

  std::vector<NodalField *> wavenb_global_tmp;

  wavenb_global_tmp.resize(this->dim);

  for (int d=0; d<this->dim; ++d) {
    wavenb_global_tmp[d] = new NodalField(this->nb_nodes_global);
  }
  
  // fill wave number arrays
  this->initWaveNumbersGlobal(wavenb_global_tmp);

  // for parallel implementation

  // save local coords provided by the user (e.g. for coupling with FEM)
  this->initCoordsLocalCostum(coords_costum);

  // sort and assing modes for load balance
  this->assignFFTModes(wavenb_global_tmp);
  // sort and assign nodes
  this->initSortCostumNodesMap();

  this->initParallel();

  for (int d=0; d<this->dim; ++d) {
    delete wavenb_global_tmp[d];
  }
}

/* --------------------------------------------------------------------------
 * 3D with default coords
 *
 * All operations in spatial domain are performed serially:
 * local==global for spatial domain
 * operations in fourier domain can be performed in parallel by default
 *
 */
Mesh::Mesh(double Lx, int Nx,
	   double Lz, int Nz) :
  dim(3),
  length_x(Lx),
  length_z(Lz),
  nb_nodes_x_global(Nx),
  nb_nodes_z_global(Nz),
  spatial_domain_parallel(false)
{
  // init global variables // allocate sorting maps
  this->initGlobal();

  std::vector<NodalField *> coords_global_tmp;
  std::vector<NodalField *> wavenb_global_tmp;

  coords_global_tmp.resize(this->dim);
  wavenb_global_tmp.resize(this->dim);

  for (int d=0; d<this->dim; ++d) {
    coords_global_tmp[d] = new NodalField(this->nb_nodes_global);
    wavenb_global_tmp[d] = new NodalField(this->nb_nodes_global);
  }

  // init coords wave number default
  this->initCoordsGlobal(coords_global_tmp);
  this->initWaveNumbersGlobal(wavenb_global_tmp);

  // for parallel implementation
  // in assignNodes() //spatial domain is serial

  // sort and assing modes for load balance
  this->assignFFTModes(wavenb_global_tmp);
  // spatial domain serial -> local=global
  this->assignGlobalNodes(coords_global_tmp);

  this->initParallel();

  for (int d=0; d<this->dim; ++d) {
    delete coords_global_tmp[d];
    delete wavenb_global_tmp[d];
  }
}

/* --------------------------------------------------------------------------
 * 3D with custom coords
 *
 * All operations is spatial domain are performed in parallel:
 * FFTW is performed in serial (for now)
 * operations in fourier domain can be performed in parallel by default
 *
 */

Mesh::Mesh(double Lx, int Nx,
	   double Lz, int Nz,
	   std::vector<NodalField *> & coords_costum) :
  dim(3),
  length_x(Lx),
  length_z(Lz),
  nb_nodes_x_global(Nx),
  nb_nodes_z_global(Nz),
  spatial_domain_parallel(true),
  is_costum(true)
{

  // init global variables // allocate sorting maps
  this->initGlobal();

  std::vector<NodalField *> wavenb_global_tmp;
  wavenb_global_tmp.resize(this->dim);
  
  for (int d=0; d<this->dim; ++d) {
      wavenb_global_tmp[d] = new NodalField(this->nb_nodes_global);
  }

  // fill wave number arrays
  this->initWaveNumbersGlobal(wavenb_global_tmp);

  // for parallel implementation

  // save local coords provided by the user (e.g. for coupling with FEM)
  this->initCoordsLocalCostum(coords_costum);

  // sort and assing modes for load balance
  this->assignFFTModes(wavenb_global_tmp);
  // sort and assign nodes
  this->initSortCostumNodesMap();

  this->initParallel();

  for (int d=0; d<this->dim; ++d) {
    delete wavenb_global_tmp[d];
  }
}

/* -------------------------------------------------------------------------- */
Mesh::~Mesh() {
  for (int d=0; d<this->dim; ++d) {
    delete this->coords[d];
  }
  for (int d=0; d<this->dim; ++d) {
    delete this->wave_numbers[d];
  }

  delete[] this->double_buffer;
  delete[] this->fftw_complex_buffer;

  if (this->is_costum)
    delete[] this->sort_costum_nodes_map;

  delete[] this->sort_fft_modes_map;
}

/* -------------------------------------------------------------------------- */
/* define global values resize and allocate global members
 */
void Mesh::initGlobal() {

  // global information
  if (this->dim==2) {
    this->nb_fft_x_global = this->nb_nodes_x_global / 2 + 1;
    this->nb_fft_z_global = 1;
  }
  else if (this->dim==3) {
    this->nb_fft_x_global = this->nb_nodes_x_global;
    this->nb_fft_z_global = this->nb_nodes_z_global / 2 + 1;
  }

  this->nb_fft_global=this->nb_fft_x_global*this->nb_fft_z_global;
  this->nb_nodes_global=this->nb_nodes_x_global*this->nb_nodes_z_global;

  this->coords.resize(this->dim);
  this->wave_numbers.resize(this->dim);

  //serial case
  this->nb_nodes_local = this->nb_nodes_global;

  this->nb_fft_local = this->nb_fft_global;
  this->max_nb_fft_pp = this->nb_fft_global;

}


/* --------------------------------------------------------------------------
 */
void Mesh::initParallel() {

  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();

  if (world_size%2!=0 && world_size>1){
    std::cerr<<"ERROR: number of mpi process needs to be even number\n"<<std::flush;
    throw;
  }

  this->fftw_complex_buffer = new fftw_complex[std::max(this->nb_fft_global,this->max_nb_fft_pp*world_size)];
  this->double_buffer       = new double[this->max_nb_costum_nodes_pp*world_size];

#ifdef UCA_VERBOSE
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (world_rank==0)
    printf("Init Parallel: nb_procs %d, proc_id %d, nb_fft_local %d, "
	   "nb_nodes_local %d\n",
	   world_size,
	   world_rank,
	   this->nb_fft_local,
	   this->nb_nodes_local);
#endif /* UCA_VERBOSE */

}


/* --------------------------------------------------------------------------
 * FFTW SERAIL METHODS
 * -------------------------------------------------------------------------- */
/*
 * default coordinate initiation as needed for FFTW serial
 *
 *              real data          complex data
 * array size   N0 x N1            2 x N0 x (N1/2+1)
 * stored as    N0 x N1            2 x N0 x (N1/2+1)
 *
 */

void Mesh::initCoordsGlobal(std::vector<NodalField*> & coords_global) {

  // determine element size
  double dx = this->getDeltaX();
  double dz = this->getDeltaZ();

  // fill coords for this mesh
  if (true) { // both are equivalent - verified
      for (int i=0; i<this->nb_nodes_x_global; ++i)
	for (int j=0; j<this->nb_nodes_z_global; ++j) {
	  int ij = i*this->nb_nodes_z_global +j;

	  (*coords_global[0])(ij) = i*dx;
	  (*coords_global[1])(ij) = 0.0; // we are on the xz plane

	  if (this->dim==3)
	    (*coords_global[2])(ij) = j*dz;
	}
    }
  else {// get rid of beause unintuitive
    for (int n=0; n<this->nb_nodes_global; ++n) {
      (*coords_global[0])(n) = n/this->nb_nodes_z_global*dx;
      (*coords_global[1])(n) = 0.0; // we are on the xz plane

      if (this->dim==3)
	(*coords_global[2])(n) = n%this->nb_nodes_z_global*dz;
    }
  }
}

/* -------------------------------------------------------------------------- */
/*
 * Serial version includes mode 0 (not used in computeStressFourierCoeff)
 * to avoid sorting before adding the zero mode back in
 * computeStressFourierCoeff()
 *
 */

void Mesh::initWaveNumbersGlobal(std::vector<NodalField*> & wave_numbers_global) {

  // fundamental mode
  double k1 = 2*M_PI / this->length_x;
  double m1 =0.0;
  if (this->dim==3)
    m1 = 2*M_PI / this->length_z;

  // compute nyquest frequency
  int f_ny_x;
  if (this->dim==2)
    f_ny_x = this->nb_fft_x_global;
  else
    f_ny_x = this->nb_fft_x_global/2+1;

  // init wave numbers global
  if (true) { //both are equivalent - verified
    for (int i=0; i<this->nb_fft_x_global; i++) for (int j=0; j<this->nb_fft_z_global; j++) {
	int ij =  i*this->nb_fft_z_global +j;

	if (this-> dim==2) {
	  (*wave_numbers_global[0])(ij) = k1 * i;
	}
	else {
	  // after nyquest frequncy modes are negative
	  // e.g. 0, 1, 2, ... ny, -ny, -ny+1, ... -2, -1
	  (*wave_numbers_global[0])(ij) = k1 * (i - (i/f_ny_x)*this->nb_fft_x_global);
	  (*wave_numbers_global[2])(ij) = m1 * j;
	}
	(*wave_numbers_global[1])(ij) = 0.0; // we re on the xz plane
      }
  }
  else { // get rid of beause unintuitive
    for (int j=0; j<this->nb_fft_global; j++) {

      (*wave_numbers_global[0])(j) = k1 * (j/this->nb_fft_z_global - ((j/this->nb_fft_z_global)/(f_ny_x))*this->nb_fft_x_global);
      (*wave_numbers_global[1])(j) = 0.0; // we re on the xz plane
      if (this->dim==3)
	(*wave_numbers_global[2])(j) = m1 * (j%this->nb_fft_z_global);
    }
  }
}

/* --------------------------------------------------------------------------
 * COSTUM MESH METHODS
 * -------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------
 * save user provided local coords int class member coords
 */
void Mesh::initCoordsLocalCostum(std::vector<NodalField *> & coords_costum) {
  this->nb_nodes_local = coords_costum[0]->getNbNodes();

#ifdef UCA_VERBOSE
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  std::cout<<"rank "<<world_rank<<" local nb nodes = "<<this->nb_nodes_local<<std::endl;
#endif /* UCA_VERBOSE */

  // find max nb nodes per process
  StaticCommunicatorMPI::getInstance()->allReduce(
	       &this->nb_nodes_local,
	       &this->max_nb_costum_nodes_pp,
	       1,
	       MPI_MAX);

#ifdef UCA_VERBOSE
  if (world_rank==0)
    std::cout<<"max local nb nodes = "<<this->max_nb_costum_nodes_pp<<std::endl;
#endif /* UCA_VERBOSE */

  for (int d=0; d<this->dim; ++d){
    this->coords[d] = new NodalField(this->max_nb_costum_nodes_pp);
    this->coords[d]->setAllValuesTo(NAN);
    double * costum = coords_costum[d]->storage();
    double * local  = this->coords[d]->storage();
    for (int n=0; n<this->nb_nodes_local; ++n){
      local[n]=costum[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
struct meshcompare {
  meshcompare(double * coord_x, double * coord_z) :
    coord_x(coord_x),
    coord_z(coord_z)
  { }

  bool operator() (int i, int j ) {
    if (std::isnan(this->coord_x[i]))
      return false;

    if (std::isnan(this->coord_x[j]))
      return true;

    if ((float)this->coord_x[i]==(float)this->coord_x[j]) // to avoid problems with numerical error in the mesh
      return this->coord_z[i]<this->coord_z[j];
    else
      return this->coord_x[i]<this->coord_x[j];
  }

  double * coord_x;
  double * coord_z;
};

/* --------------------------------------------------------------------------
 * provide array with nodal indexes - costum domain decomposition
 * needs reordering for compatibility with FFTW
 */
void Mesh::initSortCostumNodesMap() {
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  std::vector<NodalField *> coords_global_tmp;
  
  coords_global_tmp.resize(this->dim);

  // allocate global coords
  for (int d=0; d<this->dim; ++d){
    coords_global_tmp[d] = new NodalField(this->max_nb_costum_nodes_pp*world_size);
  }
  
  // allocate sort nodes now that you know needed size
  this->sort_costum_nodes_map = new int [this->max_nb_costum_nodes_pp*world_size];

  // gather all local mesh to master
  for (int d=0; d<this->dim; ++d) {
    StaticCommunicatorMPI::getInstance()->gather(
		coords_global_tmp[d]->storage(),
		this->coords[d]->storage(),
		this->max_nb_costum_nodes_pp);
  }

  if (world_rank==0) {

    meshcompare mymeshcompare(coords_global_tmp[0]->storage(),
			      coords_global_tmp[this->dim-1]->storage());

    // argsort by increasing x, z coords

    std::vector<int> sort_map_vec;

    // init vector with indexes from 0 to max_nb_nodespp*world_size
    for (int i=0; i<this->max_nb_costum_nodes_pp*world_size; ++i)
      sort_map_vec.push_back(i);

    // reorder indexes -> argsort
    std::sort(sort_map_vec.begin(), sort_map_vec.end(), mymeshcompare);

    // create sorting map
    for (int i=0; i<this->max_nb_costum_nodes_pp*world_size; ++i)
      this->sort_costum_nodes_map[i]=sort_map_vec[i];

    // verify coords
    this->checkCostumCoords(coords_global_tmp);

  }
  
  for (int d=0; d<this->dim; ++d)
    delete coords_global_tmp[d];

#ifdef UCA_USE_MPI
  StaticCommunicatorMPI::getInstance()->broadcast(
			   this->sort_costum_nodes_map,
			   this->max_nb_costum_nodes_pp*world_size);
#endif
}

/* -------------------------------------------------------------------------- */
void Mesh::checkCostumCoords(std::vector<NodalField*> & coords_global) {
  double dx = this->getDeltaX();

  double tolerance = dx*1e-6;

  // declare arrays and alloc memory
  std::vector<NodalField *> coords_global_sorted;
  std::vector<NodalField *> coords_global_ref;

  coords_global_sorted.resize(this->dim);
  coords_global_ref.resize(this->dim);

  for (int d=0; d<this->dim; ++d) {
    coords_global_sorted[d] = new NodalField(this->nb_nodes_global);
    coords_global_ref[d] = new NodalField(this->nb_nodes_global);
  }
 
  //----------------------------------------------------
  // sort coords using sorting map

  for (int d=0; d<this->dim; ++d)
    this->sortCostumNodes(coords_global[d]->storage(),coords_global_sorted[d]->storage());

  // find origin
  double x0 = (*coords_global_sorted[0])(0);
  double z0 = 0;
  if (this->dim==3)
    z0 = (*coords_global_sorted[2])(0);
  std::vector<double> origin = {x0,0,z0};

#ifdef UCA_VERBOSE
  std::cout<<"origin "<<x0<<", 0";
  if (this->dim==3)
    std::cout<<", "<<z0;
  std::cout<<std::endl;
#endif /* UCA_VERBOSE */

  //----------------------------------------------------
  // compare with default coord

  // get default coords
  this->initCoordsGlobal(coords_global_ref);

  // compare and raise error is difference exceeds tolerance
  bool test_passed = true;

  for (int n=0; n<this->nb_nodes_global; ++n) {
    for (int d=0; d<this->dim; d+=2) {
      double error = std::abs((*coords_global_sorted[d])(n)
			      - origin[d]
			      - (*coords_global_ref[d])(n));
      if (error>tolerance) {
	test_passed = false;
	if (d>0)
	  std::cerr<<n<<" error in costum mesh : is "
		   <<std::setw(6)<<      (*coords_global_sorted[0])(n) - origin[0]
		   <<std::setw(6)<<", "<<(*coords_global_sorted[d])(n) - origin[d]
		   <<" should "
		   <<std::setw(6)<<      (*coords_global_ref[0])(n)
		   <<std::setw(6)<<", "<<(*coords_global_ref[d])(n)
		   <<std::endl;
      }
    }
  }
  if (!test_passed)
    throw std::runtime_error("Error; costum mesh is not a regular grid\n");

  // free memory
  for (int d=0; d<this->dim; ++d) {
    delete coords_global_sorted[d];
    delete coords_global_ref[d];
   }
}

// END COSTUM MESH METHODS


/* --------------------------------------------------------------------------
 * SPATIAL DOMAIN Serial
 */
void Mesh::assignGlobalNodes(std::vector<NodalField*> & coords_global) {


  for (int d=0; d<this->dim; ++d){
    this->coords[d] = new NodalField(this->nb_nodes_local);

    double * global_p = coords_global[d]->storage();
    double * local_p = this->coords[d]->storage();

    for (int i=0; i<this->nb_nodes_local; ++i)
      local_p[i] = global_p[i];
  }
}


/* --------------------------------------------------------------------------
 * FFT LOAD_BALANCE
 * for small mesh sizes or very big world_size not optimal because 1st mode is very big
 */
void Mesh::assignFFTModes(std::vector<NodalField*> & wave_numbers_global) {

  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  int nb_jobs = this->nb_fft_global-1;//take away zero frequency

  //split modes_km between processes

  if (world_size>1) {

    // adapt for 3D case order modes with respect to their length

    std::vector<double> totwork(world_size, 0);
    std::vector<int> mode_assigned(this->nb_fft_global, 0); //mode_assinged[job_id] = rank
    std::vector<int> nb_modes_per_rank(world_size, 0);

    if (world_rank==0) {
      std::vector<double> work_per_mode(this->nb_fft_global,0.0);
      this->computeWorkPerMode(wave_numbers_global,work_per_mode);

      // argsort by decreasing work_per_mode
      //
      std::vector<int> sorted_modes_idx(this->nb_fft_global,0);
      // init indexes
      std::iota(sorted_modes_idx.begin(), sorted_modes_idx.end(),0);
      // sort indexes based on comparing values in work
      std::sort(sorted_modes_idx.begin(), sorted_modes_idx.end(),
		[&work_per_mode](int i1, int i2){return work_per_mode[i1]>work_per_mode[i2];});
      // end argsort

      this->printMaximumNbProc(work_per_mode);

#ifdef UCA_VERBOSE
      std::cout<<"loop over modes by decreasing work load\n";
#endif /* UCA_VERBOSE */

      for (int i=0; i<nb_jobs; ++i) {
	// loop over modes by decreasing work load
	int job_id = sorted_modes_idx[i];

	// assign job to process rank with minimum work
	int min_work_rank = std::distance(totwork.begin(),
					  std::min_element(totwork.begin(),
							   totwork.end()));
	mode_assigned[job_id] = min_work_rank;
	nb_modes_per_rank[min_work_rank]+=1;
	totwork[min_work_rank] += work_per_mode[job_id];
      }

#ifdef UCA_VERBOSE
      if (world_rank==0) {
	std::cout<<"mode_assigned, rank"<<std::endl;
	for (int i=0; i<this->nb_fft_global; ++i) 	  
	  std::cout<< i<<", "<<mode_assigned[i]<<std::endl;
      }
      for (int rank=0; rank<world_size; ++rank)
      	printf("proc_id %d, totwork %g, nb_fft_local %d\n",
	       rank, totwork[rank],nb_modes_per_rank[rank]);
#endif /* UCA_VERBOSE */

      // find max nb jobs per proc and set nb modes per proc
      this->max_nb_fft_pp = *std::max_element(nb_modes_per_rank.begin(),
					      nb_modes_per_rank.end());
      this->max_nb_fft_pp++;
      
    } // world_rank==0

    // done by everybody
    StaticCommunicatorMPI::getInstance()->broadcast(&this->max_nb_fft_pp,1);
    StaticCommunicatorMPI::getInstance()->scatter(nb_modes_per_rank.data(), &this->nb_fft_local,1);

    // build sorting array from mode_assigned
    // allocate memory now that we know the size
    this->sort_fft_modes_map = new int [std::max(this->nb_fft_global,this->max_nb_fft_pp*world_size)];

    if (world_rank==0) {
      this->sort_fft_modes_map[0]=this->nb_fft_global-1;
      for (int rank=0; rank<world_size; ++rank) {
	int j = 0;
#ifdef UCA_VERBOSE
	std::cout<<"job id, array size, sort fft modes map\n";
#endif /* UCA_VERBOSE */
	for (int i=1; i<this->nb_fft_global; ++i) {
	  int job_id = i;
	  
	  if (rank == mode_assigned[job_id]){
	    this->sort_fft_modes_map[job_id] = j + this->max_nb_fft_pp*rank;
#ifdef UCA_VERBOSE
	    std::cout<<job_id<<", " <<this->max_nb_fft_pp*world_size<<", "<<this->sort_fft_modes_map[job_id]<<std::endl;
#endif /* UCA_VERBOSE */
	    ++j;
	  }
	}
      }
#ifdef UCA_VERBOSE
      for (int i=0; i<this->nb_fft_global; ++i) {
	int job_id = i;
	printf("job_id %d, sort_fft_modes_map %d, rank %d\n",
	       job_id,
	       this->sort_fft_modes_map[job_id],
	       this->sort_fft_modes_map[job_id]/this->max_nb_fft_pp);
      }
#endif /* UCA_VERBOSE */
    } // world_rank==0

#ifdef UCA_VERBOSE
    std::cout << "broadcast max_nb_fft_pp,nb_modes_par_rank,sort_fft_modes_"
	      << std::endl << std::flush;
#endif /* UCA_VERBOSE */

    StaticCommunicatorMPI::getInstance()->broadcast(
      		           this->sort_fft_modes_map,
			   this->max_nb_fft_pp*world_size);
#ifdef UCA_VERBOSE
    printf("proc_id %d, max_nb_fft_pp %d, nb_fft_local %d\n",
	   world_rank, this->max_nb_fft_pp,this->nb_fft_local);
#endif /* UCA_VERBOSE */

    for (int d=0; d<this->dim; ++d){
      // allocate local wave number memory
      this->wave_numbers[d] = new NodalField(this->max_nb_fft_pp);

      // actually sort modes and scatter them
      this->sortAndScatterFFTModes(wave_numbers_global[d]->storage(),
				   this->wave_numbers[d]->storage());
    }
  } // end if world_size>1
  else { //world_size==0
    //serial local = global
    //just to avoid segmentation fault by deleteing an non new pointer
    this->sort_fft_modes_map = new int [1];
    for (int d=0; d<this->dim; ++d) {
      this->wave_numbers[d] = new NodalField(this->nb_fft_local);

      double * global_p = wave_numbers_global[d]->storage();
      double * local_p = this->wave_numbers[d]->storage();

      for (int i=0; i<this->nb_fft_local; ++i) {
	local_p[i] = global_p[i];
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::computeWorkPerMode(std::vector<NodalField*> & wave_numbers_global,
			 std::vector<double> & work_per_mode) {

  // estimate work_per_mode = nb_intrgration_int
  // k,m are components of mode vector q
  for (int i=1; i<this->nb_fft_global; ++i) {
    // took away zero frequency
    double k = (*wave_numbers_global[0])(i);
    double m = 0.0;
    if (this->dim==3)
      m = (*wave_numbers_global[2])(i);
    double q = std::sqrt(k*k + m*m);
    work_per_mode[i] = 1.0 / q;
  }
  work_per_mode[0]=0.0; // no convolution for 0 mode
#ifdef UCA_VERBOSE
  std::cout<<"work per mode"<<std::endl;
  for (int i=0; i<this->nb_fft_global; ++i) {
    std::cout<< work_per_mode[i]<<std::endl;
  }
#endif

}

/* -------------------------------------------------------------------------- */
int Mesh::getMaximumNbProc(std::vector<double> & work_per_mode) {
  /* For 2D:
   * N0 is the number of nodes
   * N0/2+1 is the number of FFT modes
   *
   * The total work during convolution is given by the nth partial sum of the
   * harmonic series,
   * H_n = sum_{k=1}^{n} 1/k
   * where n=N0/2+1
   *
   * For 3D:
   * N0 is the number of nodes along x
   * N1 is the number of nodes along z
   * the number of FFT Modes is N0x(N1/2+1)
   * the total work is give by
   *
   * H = sum_i sum_j 1/q_ij
   * where q_ij = sqrt(k_i^2+m_j^2)
   *
   * for 3D the optimal number of process quickly exceeds the available processes.
   * This is not an issue -> In 3D the are practically no load balance issues.
   * All processes will be assigned the same amount of work.
   *
   * Since our current parallel scheme does not allow us to subdivide the 1st
   * mode between different processes we can prove that the theoretical maximum
   * speedup due to parallelization (neglecting communication overhead) is:
   *
   * Speedup = SetialTime/ParallelTime <= TotalWork/Work1stMode = MaxSpeedup
   *
   * and it is achieved when one process only has the first mode. Therefore, the
   * maximum number of processes is also equivalent to TotalWork/Work1stMode.
   *
   */

  double tot_work = std::accumulate(work_per_mode.begin(),work_per_mode.end(),0.0);
  double max_work_per_mode = *std::max_element(work_per_mode.begin(),work_per_mode.end());

#ifdef UCA_VERBOSE
  std::cout<<"tot work : "<<tot_work<<", max w :"<<max_work_per_mode<<std::endl;
#endif /* UCA_VERBOSE */

  return int(tot_work/max_work_per_mode) + 1;// +1 coz intger division rounds down
}

/* -------------------------------------------------------------------------- */
void Mesh::printMaximumNbProc(std::vector<double> & work_per_mode) {
  if (StaticCommunicatorMPI::getInstance()->whoAmI()==0) {
    std::cout << "---------------------------------------------------------"
	      << std::endl;
    printf("Maximal number of MPI processes for optimal load balance: \n");
    printf("MaxNbProc = %d \n",
	   this->getMaximumNbProc(work_per_mode));
    std::cout << "---------------------------------------------------------"
	      << std::endl;
  }
}

/* --------------------------------------------------------------------------
 * SORT
 * -------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------
 * SORT FFT MODES
 * -------------------------------------------------------------------------- */
template <typename T>
void Mesh::sortFFT(T* un_sorted, T * sorted, int root_rank) {
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (world_rank==root_rank)
    for (int n=0; n<this->nb_fft_global; ++n)
      sorted[n] = un_sorted[this->sort_fft_modes_map[n]];
}
/* -------------------------------------------------------------------------- */
template <typename T>
void Mesh::unsortFFT(T * sorted, T * un_sorted, int root_rank) {
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (world_rank==root_rank)
    for (int n=0; n<this->nb_fft_global; ++n)
      un_sorted[this->sort_fft_modes_map[n]] = sorted[n];
}
/* --------------------------------------------------------------------------
 * SORT COSTUM NODES
 * -------------------------------------------------------------------------- */
template <typename T>
void Mesh::sortCostumNodes(T* un_sorted, T * sorted, int root_rank) {
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (world_rank==root_rank)
    for (int n=0; n<this->nb_nodes_global; ++n)
      sorted[n] = un_sorted[this->sort_costum_nodes_map[n]];
}
/* -------------------------------------------------------------------------- */
template <typename T>
void Mesh::unsortCostumNodes(T * sorted, T * un_sorted, int root_rank) {
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (world_rank==root_rank)
    for (int n=0; n<this->nb_nodes_global; ++n)
      un_sorted[this->sort_costum_nodes_map[n]] = sorted[n];
}

/* --------------------------------------------------------------------------
 * template specialization for fftw complex
 * -------------------------------------------------------------------------- */
template <>
void Mesh::sortFFT(fftw_complex * un_sorted, fftw_complex  * sorted, int root_rank) {
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (world_rank==root_rank)
    for (int n=0; n<this->nb_fft_global; ++n)
      for (int j=0;j<2; ++j)
	sorted[n][j] = un_sorted[this->sort_fft_modes_map[n]][j];
}
/* -------------------------------------------------------------------------- */
template <>
void Mesh::unsortFFT(fftw_complex * sorted, fftw_complex * un_sorted, int root_rank) {
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (world_rank==root_rank)
    for (int n=0; n<this->nb_fft_global; ++n)
      for (int j=0;j<2; ++j)
	un_sorted[this->sort_fft_modes_map[n]][j] = sorted[n][j];
}

/* --------------------------------------------------------------------------
 * SORT SCATTER GATHER FFT MODE TEMPLATE
 * -------------------------------------------------------------------------- */
template <typename T>
void Mesh::sortAndScatterFFTModes(T * Uglobal, T * Ulocal, T * buffer, int size, int root_rank) {
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();

  if (world_size > 1) {
    this->unsortFFT(Uglobal, buffer,root_rank);

    StaticCommunicatorMPI::getInstance()->scatter(
			   buffer, Ulocal, size, root_rank);
  }
}
/* -------------------------------------------------------------------------- */
template <typename T>
void Mesh::gatherAndSortFFTModes(T * Uglobal, T * Ulocal, T * buffer, int size, int root_rank) {
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();

  if (world_size > 1) {
    StaticCommunicatorMPI::getInstance()->gather(
			   buffer, Ulocal, size, root_rank);

    this->sortFFT(buffer,Uglobal,root_rank);
  }
}


/* --------------------------------------------------------------------------
 * SORT SCATTER GATHER COSTUM NODES TEMPLATE
 * -------------------------------------------------------------------------- */
template <typename T>
void Mesh::sortAndScatterCostumNodes(T * Uglobal, T * Ulocal, T * buffer, int size, int root_rank) {
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();

  this->unsortCostumNodes(Uglobal, buffer, root_rank);
    if (world_size > 1) {
    StaticCommunicatorMPI::getInstance()->scatter(
			   buffer, Ulocal, size, root_rank);
  }
  else {
    for (int n=0; n<this->nb_nodes_global; ++n)
      Ulocal[n] = buffer[n];
  }
}
/* -------------------------------------------------------------------------- */
template <typename T>
void Mesh::gatherAndSortCostumNodes(T * Uglobal,
				    T * Ulocal,
				    T * buffer,
				    int size,
				    int root_rank) {
  
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();

  if (world_size > 1) {
    StaticCommunicatorMPI::getInstance()->gather(
  			   buffer, Ulocal, size, root_rank);
  }
  else {
    for (int n=0; n<this->nb_nodes_global; ++n)
      buffer[n]=Ulocal[n];
  }
  this->sortCostumNodes(buffer,Uglobal,root_rank);
}
/* instead of template specialization there are a lot of functions calling
 * the templated one
 */

/* -------------------------------------------------------------------------- */
void Mesh::sortAndScatterCostumNodes(int * Uglobal, int * Ulocal, int root_rank) {
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  int * bufferint=NULL; // we are allocating and freeing memory here just because this function is used only during initialization

  if (world_rank==root_rank)
    bufferint = new int[(this->max_nb_costum_nodes_pp)*(world_size)]();

  this->sortAndScatterCostumNodes(Uglobal,Ulocal,bufferint,this->max_nb_costum_nodes_pp,root_rank);

  if (world_rank==root_rank) delete[] bufferint;
}
/* -------------------------------------------------------------------------- */
void Mesh::gatherAndSortCostumNodes(int * Uglobal, int * Ulocal, int root_rank) {
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  int * bufferint=NULL; // we are allocating and freeing memory here just because this function is used only during initialization

  if (world_rank==root_rank)
    bufferint = new int[(this->max_nb_costum_nodes_pp)*(world_size)]();

  this->gatherAndSortCostumNodes(Uglobal,Ulocal,bufferint,this->max_nb_costum_nodes_pp,root_rank);

  if (world_rank==root_rank) delete[] bufferint;
}
/* -------------------------------------------------------------------------- */
void Mesh::sortAndScatterCostumNodes(double * Uglobal, double * Ulocal, int root_rank) {
  this->sortAndScatterCostumNodes(Uglobal,Ulocal,this->double_buffer,this->max_nb_costum_nodes_pp,root_rank);
}
/* -------------------------------------------------------------------------- */
void Mesh::gatherAndSortCostumNodes(double * Uglobal,
				    double * Ulocal,
				    int root_rank) {

  this->gatherAndSortCostumNodes(Uglobal,
				 Ulocal,
				 this->double_buffer,
				 this->max_nb_costum_nodes_pp,
				 root_rank);
}

/* --------------------------------------------------------------------------
 * SORT AND SCATTER FFTW MODES
 * -------------------------------------------------------------------------- */
void Mesh::sortAndScatterFFTModes(int * Uglobal, int * Ulocal, int root_rank) {
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  int * bufferint=NULL; // we are allocating and freeing memory here just because this function is used only during initialization

  if (world_rank==root_rank)
    bufferint = new int[(this->max_nb_fft_pp)*(world_size)]();

  this->sortAndScatterFFTModes(Uglobal,Ulocal,bufferint,this->max_nb_fft_pp,root_rank);
  if (world_rank==root_rank) delete[] bufferint;
}

/* -------------------------------------------------------------------------- */
void Mesh::sortAndScatterFFTModes(double * Uglobal, double * Ulocal, int root_rank) {
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  double * buffer=NULL;

  if (world_rank==root_rank)
    buffer = new double[(this->max_nb_fft_pp)*(world_size)]();

  this->sortAndScatterFFTModes(Uglobal, Ulocal, buffer, this->max_nb_fft_pp,root_rank);
  if (world_rank==root_rank) delete[] buffer;
}

/* -------------------------------------------------------------------------- */
void Mesh::sortAndScatterFFTModes(fftw_complex * Uglobal, fftw_complex * Ulocal, int root_rank) {
  // 2* accounts for real and imag parts -> fftw_complex aka double[2]
  this->sortAndScatterFFTModes(Uglobal,Ulocal,this->fftw_complex_buffer,(this->max_nb_fft_pp)*2,root_rank);
}

/* --------------------------------------------------------------------------
 * GATHER AND SORT FFTW MODES
 * -------------------------------------------------------------------------- */
void Mesh::gatherAndSortFFTModes(fftw_complex * Uglobal, fftw_complex * Ulocal, int root_rank) {
  // 2* accounts for real and imag parts -> fftw_complex aka double[2]
  this->gatherAndSortFFTModes(Uglobal,Ulocal,this->fftw_complex_buffer,(this->max_nb_fft_pp)*2,root_rank);
}

// end parallel methods

__END_UGUCA__
