/**
 * @file   static_communicator_mpi.cc
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
#include "static_communicator_mpi.hh"
#include <iostream>

__BEGIN_UGUCA__

typedef double fftw_complex[2];

#ifdef UCA_USE_MPI
class MPITypeWrapper {
public:
  template<typename T>
  static inline MPI_Datatype getMPIDatatype();
};
#else
#define UNUSED(...) (void)(__VA_ARGS__)
#endif /* UCA_USE_MPI */

/* -------------------------------------------------------------------------- */
/* implementation (needs to be in same file because of template)              */
/* -------------------------------------------------------------------------- */

bool StaticCommunicatorMPI::is_instantiated = false;
StaticCommunicatorMPI* StaticCommunicatorMPI::static_communicator = NULL;

/* -------------------------------------------------------------------------- */
StaticCommunicatorMPI::StaticCommunicatorMPI() {

  int is_initialized=0;
#ifdef UCA_USE_MPI
  MPI_Initialized(&is_initialized);
  if (!is_initialized)
    MPI_Init(NULL,NULL);
  MPI_communicator=MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_communicator, &world_rank);
  MPI_Comm_size(MPI_communicator, &world_size);

#ifdef UCA_VERBOSE
  if ((!is_initialized) && (world_rank==0)) printf("MPI initialized \n");
#endif /* UCA_VERBOSE */
#endif /* UCA_USE_MPI */

  this->is_externaly_initialized = is_initialized;
  this->is_instantiated = true;
}

/* -------------------------------------------------------------------------- */
StaticCommunicatorMPI::~StaticCommunicatorMPI() {
  this->finalize();
  is_instantiated = false;
  StaticCommunicatorMPI::static_communicator = NULL;
}

/* -------------------------------------------------------------------------- */
void StaticCommunicatorMPI::finalize() {
 #ifdef UCA_USE_MPI
  if(!this->is_externaly_initialized){
#ifdef UCA_VERBOSE
    if (world_rank==0) printf("MPI finalize \n");
#endif /* UCA_VERBOSE */
    MPI_Finalize();
  }
 #endif /* UCA_USE_MPI */
}

/* -------------------------------------------------------------------------- */
StaticCommunicatorMPI * StaticCommunicatorMPI::getInstance() {
  if (!static_communicator)
    static_communicator = new StaticCommunicatorMPI();
  return static_communicator;
}
/* -------------------------------------------------------------------------- */
unsigned int StaticCommunicatorMPI::whoAmI()    const { return this->world_rank; }
unsigned int StaticCommunicatorMPI::getNbProc() const { return this->world_size; }
unsigned int StaticCommunicatorMPI::getRoot()   const { return 0; }
/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::send(T * buffer, int size, int receiver, int tag) {
#ifdef UCA_USE_MPI
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
  MPI_Send(buffer, size, type, receiver, tag, MPI_communicator);
#else
  UNUSED(buffer);
  UNUSED(size);
  UNUSED(receiver);
  UNUSED(tag);
#endif /* UCA_USE_MPI */
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::recv(T * buffer, int size, int sender, int tag) {
#ifdef UCA_USE_MPI
  MPI_Status status;
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
  MPI_Recv(buffer, size, type, sender, tag, MPI_communicator, &status);
#else
  UNUSED(buffer);
  UNUSED(size);
  UNUSED(sender);
  UNUSED(tag);
#endif /* UCA_USE_MPI */
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::scatter(T * root_buffer,
				    T * leaf_buffer,
				    int size_per_proc,
				    int root_rank) {
#ifdef UCA_USE_MPI
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
  MPI_Scatter(root_buffer, size_per_proc, type, leaf_buffer, size_per_proc, type, root_rank,
	      MPI_communicator);
#else
  UNUSED(root_buffer);
  UNUSED(leaf_buffer);
  UNUSED(size_per_proc);
  UNUSED(root_rank);
#endif /* UCA_USE_MPI */
}
template <typename T>
void StaticCommunicatorMPI::scatter(T * buffer,
				    int size_per_proc,
				    int root_rank) {
#ifdef UCA_USE_MPI
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
  MPI_Scatter(buffer, size_per_proc, type, buffer, size_per_proc, type, root_rank,
	      MPI_communicator);
#else
  UNUSED(buffer);
  UNUSED(size_per_proc);
  UNUSED(root_rank);
#endif /* UCA_USE_MPI */
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::gather(T * buffer, int size, int root_rank) {
#ifdef UCA_USE_MPI
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
  MPI_Gather(buffer, size, type, buffer, size, type, root_rank,
	     MPI_communicator);
#else
  UNUSED(buffer);
  UNUSED(size);
  UNUSED(root_rank);
#endif /* UCA_USE_MPI */
}

template <typename T>
void StaticCommunicatorMPI::gather(T * root_buffer,
				   T * leaf_buffer,
				   int size,
				   int root_rank) {
#ifdef UCA_USE_MPI
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
  MPI_Gather(leaf_buffer, size, type, root_buffer, size, type, root_rank,
	     MPI_communicator);
#else
  for (int i=0; i<size; ++i) {
    *root_buffer = *leaf_buffer;
    ++root_buffer;
    ++leaf_buffer;
  }
  UNUSED(root_rank);
#endif /* UCA_USE_MPI */
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::allGather(T * buffer, int size) {
#ifdef UCA_USE_MPI
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
  MPI_Allgather(buffer, size, type, buffer, size, type,
		MPI_communicator);
  if (world_size==2) {
#ifdef UCA_VERBOSE
    std::cout << "Warning: allGather error when using same send and recv buffer "
	      << "and world_size==2" << std::endl;;
    std::cout << "... work around: Broadcast correct result from root" << std::endl;
#endif /* UCA_VERBOSE */
    MPI_Bcast(buffer, size*world_size, type, 0, MPI_communicator);
  }
#else
  UNUSED(buffer);
  UNUSED(size);
#endif /* UCA_USE_MPI */
}
template <typename T>
void StaticCommunicatorMPI::allGather(T * root_buffer, T * leaf_buffer, int size) {
#ifdef UCA_USE_MPI
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
  MPI_Allgather(leaf_buffer, size, type, root_buffer, size, type,
		MPI_communicator);
#else
  UNUSED(root_buffer);
  UNUSED(leaf_buffer);
  UNUSED(size);
#endif /* UCA_USE_MPI */
}
/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::broadcast(T * buffer, int size, int root_rank) {
#ifdef UCA_USE_MPI
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
  MPI_Bcast(buffer, size, type, root_rank, MPI_communicator);
#else
  UNUSED(buffer);
  UNUSED(size);
  UNUSED(root_rank);
#endif /* UCA_USE_MPI */
}

/* -------------------------------------------------------------------------- */
template <typename T>
void StaticCommunicatorMPI::allReduce(T * send_data, T * recv_data,
				      int count,
				      const MPI_Op & op) {
#ifdef UCA_USE_MPI
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<T>();
  MPI_Allreduce(send_data, recv_data, count, type, op, MPI_communicator);
#else
  for (int i=0; i<count; ++i) {
    *recv_data = *send_data;
    ++recv_data;
    ++send_data;
  }
  UNUSED(op);
#endif /* UCA_USE_MPI */
}

/* -------------------------------------------------------------------------- */
/* template specialization                                                    */
/* -------------------------------------------------------------------------- */
#ifdef UCA_USE_MPI
template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<char>() {
  return MPI_CHAR;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<float>() {
  return MPI_FLOAT;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<double>() {
  return MPI_DOUBLE;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<double[2]>() {
  return MPI_DOUBLE;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<long double>() {
  return MPI_LONG_DOUBLE;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<int>() {
  return MPI_INT;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<unsigned int>() {
  return MPI_UNSIGNED;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<signed long int>() {
  return MPI_LONG;
}

template <> MPI_Datatype MPITypeWrapper::getMPIDatatype<unsigned long int>() {
  return MPI_UNSIGNED_LONG;
}

template <>
MPI_Datatype MPITypeWrapper::getMPIDatatype<signed long long int>() {
  return MPI_LONG_LONG;
}

template <>
MPI_Datatype MPITypeWrapper::getMPIDatatype<unsigned long long int>() {
  return MPI_UNSIGNED_LONG_LONG;
}
#endif /* UCA_USE_MPI */


template void StaticCommunicatorMPI::scatter<int>(int*, int, int);
template void StaticCommunicatorMPI::scatter<int>(int*, int*, int, int);
template void StaticCommunicatorMPI::scatter<unsigned int>(unsigned int*, int, int);
template void StaticCommunicatorMPI::scatter<unsigned int>(unsigned int*, unsigned int*, int, int);
template void StaticCommunicatorMPI::scatter<float>(float*, int, int);
template void StaticCommunicatorMPI::scatter<float>(float*, float*, int, int);
template void StaticCommunicatorMPI::scatter<double>(double*, int, int);
template void StaticCommunicatorMPI::scatter<double>(double*, double*, int, int);
template void StaticCommunicatorMPI::scatter<fftw_complex>(fftw_complex*, int, int);
template void StaticCommunicatorMPI::scatter<fftw_complex>(fftw_complex*, fftw_complex*, int, int);

template void StaticCommunicatorMPI::gather<int>(int*, int, int);
template void StaticCommunicatorMPI::gather<int>(int*, int*, int, int);
template void StaticCommunicatorMPI::gather<float>(float*, int, int);
template void StaticCommunicatorMPI::gather<float>(float*, float*, int, int);
template void StaticCommunicatorMPI::gather<double>(double*, int, int);
template void StaticCommunicatorMPI::gather<double>(double*, double*, int, int);
template void StaticCommunicatorMPI::gather<fftw_complex>(fftw_complex*, int, int);
template <> void StaticCommunicatorMPI::gather(fftw_complex * root_buffer,
					       fftw_complex * leaf_buffer,
					       int size,
					       int root_rank) {
#ifdef UCA_USE_MPI
  MPI_Datatype type = MPITypeWrapper::getMPIDatatype<fftw_complex>();
  MPI_Gather(leaf_buffer, size, type, root_buffer, size, type, root_rank,
	     MPI_communicator);
#else
  for (int i=0; i<size; ++i) {
    *root_buffer[0] = *leaf_buffer[0];
    *root_buffer[1] = *leaf_buffer[1];
    ++root_buffer;
    ++leaf_buffer;
  }
  UNUSED(root_rank);
#endif /* UCA_USE_MPI */
}

template void StaticCommunicatorMPI::allGather<float>(float*, int);
template void StaticCommunicatorMPI::allGather<float>(float*, float*, int);
template void StaticCommunicatorMPI::allGather<fftw_complex>(fftw_complex*, int);
template void StaticCommunicatorMPI::allGather<fftw_complex>(fftw_complex*, fftw_complex*, int);
template void StaticCommunicatorMPI::allGather<double>(double*, int);
template void StaticCommunicatorMPI::allGather<double>(double*, double*, int);

template void StaticCommunicatorMPI::recv<float>(float*, int, int, int);
template void StaticCommunicatorMPI::recv<double>(double*, int, int, int);
template void StaticCommunicatorMPI::recv<int>(int*, int, int, int);
template void StaticCommunicatorMPI::send<float>(float*, int, int, int);
template void StaticCommunicatorMPI::send<double>(double*, int, int, int);
template void StaticCommunicatorMPI::send<int>(int*, int, int, int);

template void StaticCommunicatorMPI::broadcast<float>(float*, int, int);
template void StaticCommunicatorMPI::broadcast<double>(double*, int, int);
template void StaticCommunicatorMPI::broadcast<int>(int*, int, int);
template void StaticCommunicatorMPI::broadcast<unsigned int>(unsigned int*, int, int);

template void StaticCommunicatorMPI::allReduce<int>(int*, int*, int, const MPI_Op &);
template void StaticCommunicatorMPI::allReduce<unsigned int>(unsigned int*, unsigned int*, int, const MPI_Op &);

__END_UGUCA__

