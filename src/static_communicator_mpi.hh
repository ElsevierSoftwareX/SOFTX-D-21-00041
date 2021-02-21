/**
 * @file   static_communicator_mpi.hh
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
#ifndef __STATIC_COMMUNICATOR_MPI_H__
#define __STATIC_COMMUNICATOR_MPI_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include <cstddef>
#include <cstdio>

#ifdef UCA_USE_MPI
#include <mpi.h>
#else
typedef int MPI_Op;
static const MPI_Op MPI_MAX=0;
#endif

__BEGIN_UGUCA__


/*
 *
 * Usage in your main file:
 *
 *
 * # Option 1
 *
 * return pointer of the staticCommunicator instance without saving it
 *
 *    int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
 *    StaticCommunicatorMPI::getInstance()->doStuff();
 *
 * call finalize at the end of your main file
 *
 *    StaticCommunicatorMPI::getInstance()->finalize();
 *
 *
 *
 * # Option 2
 *
 * save the pointer in your main file
 *
 *    StaticCommunicatorMPI * comm  = StaticCommunicatorMPI::getInstance();
 *    int world_rank = comm->whoAmI();
 *    comm->doStuff();
 *
 * delete pointer at the end fo your main file
 *
 *    delete comm; //calls finalize()
 *
 */
/* -------------------------------------------------------------------------- */

class StaticCommunicatorMPI {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
protected:
  StaticCommunicatorMPI();

public:

  ~StaticCommunicatorMPI();
  void finalize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void barrier(){
#ifdef UCA_USE_MPI
    MPI_Barrier(MPI_communicator);
#endif
  };

  template <typename T> void send(T * buffer, int size, int receiver, int tag);

  template <typename T> void recv(T * buffer, int size, int sender, int tag);

  template <typename T> void scatter(T * buffer, int size_per_proc, int root_rank=0);

  template <typename T> void scatter(T * root_rank_buffer, T* leaf_buffer, int size_per_proc, int root_rank=0);

  template <typename T> void  gather(T * buffer, int size, int root_rank=0);

  template <typename T> void  gather(T* root_rank_buffer, T * leaf_buffer, int size, int root_rank=0);
  
  template <typename T> void allGather(T * buffer, int size);

  template <typename T> void allGather(T * root_rank_buffer, T* leaf_buffer, int size);

  template <typename T> void broadcast(T * values, int nb_values, int root_rank=0);

  template <typename T> void allReduce(T * send_data, T * recv_data,  int count,
				       const MPI_Op & op);

  /*
  template <typename T> void reduce(T * values, int nb_values,
				    const MPI_Op & op, int root_rank);
  */

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  unsigned int whoAmI() const;
  unsigned int getNbProc() const;
  unsigned int getRoot() const;

  static StaticCommunicatorMPI * getInstance();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // MPI common world
  int world_size=1;
  int world_rank=0;

private:

  static StaticCommunicatorMPI *static_communicator;
  static bool is_instantiated;

#ifdef UCA_USE_MPI
  MPI_Comm MPI_communicator;
#endif

  bool is_externaly_initialized;
};

__END_UGUCA__

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "static_communicator_mpi_inline_impl.hh"


#endif /* __STATIC_COMMUNICATOR_MPI_H__ */

