/**
 * @file   test_static_communicator_mpi.cc
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
#include <iostream>
#include <iomanip>
#include <fstream>

typedef double fftw_complex[2];

#include "static_communicator_mpi.hh"

using namespace uguca;

int main() {

  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();

  int root = StaticCommunicatorMPI::getInstance()->getRoot();
  bool error = false;
  std::cout<<world_rank<<std::endl;

  //---------------------------------------------------
  // SEND RECV

  if (world_rank==root)
    std::cout<<"test send recv \n"<<std::flush;

  int size =16;
  float* buffer = new float[size];

  // init buffer on root only and send to leaf 1
  if (world_rank==root) {
    std::cout<<"rank 0 send:   "<<std::flush;
    for (int i=0; i<size;i++){
      buffer[i]=i;
      std::cout<<std::setw(3)<<buffer[i];
    }
    std::cout << std::endl;
    StaticCommunicatorMPI::getInstance()->send(buffer,size,1,0);
  }
  // receive buffer on leaf 1
  if (world_rank==1) {
    StaticCommunicatorMPI::getInstance()->recv(buffer,size,root,0);
    std::cout << "rank 1 receive:";
    for (int i=0; i<size;i++) {
      std::cout<<std::setw(3)<<buffer[i];
      if (buffer[i]!=i) error=true;
    }
    std::cout<<std::endl;
  }
  if (error) return 1;

  StaticCommunicatorMPI::getInstance()->barrier();

  if (world_rank==root)
    std::cout<<"send recv success!\n"<<std::flush;

  delete buffer;

  //---------------------------------------------------
  //---------------------------------------------------
  // SCATTER in-place
  if (world_rank==root)
    std::cout<<"test scatter in-place\n"<<std::flush;

  StaticCommunicatorMPI::getInstance()->barrier();

  int size_pproc=size/world_size;
  if (world_rank==root && size%world_size!=0)
    printf("WARNING wold_size is not a multiple of problem_size: %d, %d\n",world_size,size);
  std::cout<<std::flush;

  // init buffer on root only
  if (world_rank==root) {
    buffer = new float[size];
    for (int i=0; i<size; i++)
      buffer[i]=i;
  }
  else {
    buffer = new float[size_pproc];
  }

  // scatter from root to leafs
  StaticCommunicatorMPI::getInstance()->scatter(buffer,size_pproc);

  StaticCommunicatorMPI::getInstance()->barrier();

  printf("rank %d receive: ",world_rank);
  for (int i=0; i<size_pproc; i++){
    std::cout<<std::setw(3)<<buffer[i];
    if (buffer[i]!=world_rank*size_pproc+i)
      error=true;
  }
  std::cout<<std::endl<<std::flush;
  if (error) return 1;

  StaticCommunicatorMPI::getInstance()->barrier();

  if (world_rank==root) printf("scatter in-place success!\n") ;


  //---------------------------------------------------
  // GATHER in-place
  StaticCommunicatorMPI::getInstance()->barrier();
  if (world_rank==root)
    std::cout<<"test gather in-place \n"<<std::flush ;
  //---------------------------------------------------
  // change buffer
  for (int i=0; i<size_pproc; i++)
    buffer[i]+=i+world_rank;

  printf("rank %d change:  ",world_rank);
  for (int i=0; i<size_pproc; i++)
    std::cout<<std::setw(3)<<buffer[i];
  std::cout<<std::endl<<std::flush;

  StaticCommunicatorMPI::getInstance()->gather(buffer,size_pproc);

  if (world_rank==root) {
    std::cout<<"rank "<<world_rank<<" gather:  ";
    for (int i=0; i<size;i++) {
      std::cout<<std::setw(3)<<buffer[i];
      if (buffer[i]!=i+i%size_pproc+i/size_pproc){
	error=true;
	std::cout<<"error ";
      }
    }
    std::cout<<std::endl<<std::flush;
  }
  if (error) return 1;

  StaticCommunicatorMPI::getInstance()->barrier();

  if (world_rank==root) printf("gather in-place success!\n") ;

  //---------------------------------------------------
  //---------------------------------------------------
  // SCATTER 2
  StaticCommunicatorMPI::getInstance()->barrier();
  if (world_rank==root) printf("test scatter not-in-place\n") ;
  StaticCommunicatorMPI::getInstance()->barrier();

  size_pproc=size/world_size;
  if (world_rank==root && size%world_size!=0)
    printf("WARNING wold_size is not a multiple of problem_size: %d, %d\n",world_size,size);

  float * root_buffer = new float[size];;
  if (world_rank==root) {
    for (int i=0; i<size; i++)
      root_buffer[i]=i;
  }

  buffer = new float[size_pproc];

  StaticCommunicatorMPI::getInstance()->scatter(root_buffer,buffer,size_pproc);

  StaticCommunicatorMPI::getInstance()->barrier();

  printf("rank %d receive: ",world_rank);
  for (int i=0; i<size_pproc; i++){
    std::cout<<std::setw(3)<<buffer[i];
    if (buffer[i]!=world_rank*size_pproc+i) {
      error=true;
      std::cout<<"error ";
    }
  }
  std::cout<<std::endl<<std::flush;
  if (error) return 1;

  StaticCommunicatorMPI::getInstance()->barrier();

  if (world_rank==root) printf("scatter not-in-place success!\n") ;

  StaticCommunicatorMPI::getInstance()->barrier();

  //---------------------------------------------------
  // change buffer
  for (int i=0; i<size_pproc; i++)
    buffer[i]+=i+world_rank;

  printf("rank %d change:  ",world_rank);
  for (int i=0; i<size_pproc; i++)
    std::cout<<std::setw(3)<<buffer[i];
  std::cout<<std::endl<<std::flush;
  StaticCommunicatorMPI::getInstance()->barrier();

  //---------------------------------------------------
  // GATHER not in place
  if (world_rank==root) printf("test gather not-in-place\n") ;
  StaticCommunicatorMPI::getInstance()->barrier();

  StaticCommunicatorMPI::getInstance()->gather(root_buffer,buffer,size_pproc);

  StaticCommunicatorMPI::getInstance()->barrier();

  if (world_rank==root) {
    std::cout<<"rank "<<world_rank<<" gather:  ";
    for (int i=0; i<size;i++) {
      std::cout<<std::setw(3)<<root_buffer[i];
      if (root_buffer[i]!=i+i%size_pproc+i/size_pproc)
	error=true;
    }
    std::cout<<std::endl<<std::flush;
  }
  if (error) return 1;
  StaticCommunicatorMPI::getInstance()->barrier();
  if (world_rank==root) printf("gather not-in-place success!\n") ;

  //---------------------------------------------------
  //---------------------------------------------------
  // SCATTER in-place - all gather in-place

  StaticCommunicatorMPI::getInstance()->barrier();
  delete buffer;

  if (world_rank==root) printf("test allGather in-place\n") ;
  StaticCommunicatorMPI::getInstance()->barrier();

  size_pproc=size/world_size;
  if (world_rank==root && size%world_size!=0)
    printf("WARNING wold_size is not a multiple of problem_size: %d, %d\n",world_size,size);


  buffer = new float[size]; // all need same length
  if (world_rank==root) {
    for (int i=0; i<size; i++)
      buffer[i]=i;
  }

  StaticCommunicatorMPI::getInstance()->barrier();
  if (world_rank==root) printf("first scatter \n") ;
  StaticCommunicatorMPI::getInstance()->barrier();

  StaticCommunicatorMPI::getInstance()->scatter(buffer,size_pproc);

  printf("rank %d receive:    ",world_rank);
  for (int i=0; i<size_pproc; i++){
    std::cout<<std::setw(3)<<buffer[i];
  }
  std::cout<<std::endl<<std::flush;

  //---------------------------------------------------
  // change buffer
  for (int i=0; i<size_pproc; i++)
    buffer[i]+=i+world_rank;

  printf("rank %d change :    ",world_rank);
  for (int i=0; i<size_pproc; i++)
    std::cout<<std::setw(3)<<buffer[i];
  std::cout<<std::endl<<std::flush;
  StaticCommunicatorMPI::getInstance()->barrier();

  //---------------------------------------------------
  // ALL GATHER in place
  StaticCommunicatorMPI::getInstance()->barrier();

  StaticCommunicatorMPI::getInstance()->allGather(buffer,size_pproc);

  StaticCommunicatorMPI::getInstance()->barrier();


  std::cout<<"rank "<<world_rank<<" allGather : ";
  for (int i=0; i<size;i++) {
    std::cout<<std::setw(3)<<buffer[i];
    if (buffer[i]!=i+i%size_pproc+i/size_pproc){
      error=true;
      std::cout<<"error ";
    }
  }
  //  if (error) return 1;
  std::cout<<std::endl<<std::flush;
  error=false;

  StaticCommunicatorMPI::getInstance()->barrier();
  if (world_rank==root) printf("allGather in-place success!\n") ;

  //---------------------------------------------------
  //---------------------------------------------------
  // SCATTER  allGather not-in-place

  StaticCommunicatorMPI::getInstance()->barrier();
  if (world_rank==root) printf("test allGather not-in-place \n") ;
  StaticCommunicatorMPI::getInstance()->barrier();
  delete buffer;
  delete root_buffer;

  size_pproc=size/world_size;
  if (world_rank==root && size%world_size!=0)
    printf("WARNING wold_size is not a multiple of problem_size: %d, %d\n",world_size,size);


  root_buffer = new float[size]; // all need same length
  if (world_rank==root) {
    for (int i=0; i<size; i++)
      root_buffer[i]=i;
  }

  buffer = new float[size_pproc];
  StaticCommunicatorMPI::getInstance()->barrier();
  if (world_rank==root) printf("scatter first \n") ;
  StaticCommunicatorMPI::getInstance()->scatter(root_buffer,buffer,size_pproc);

  StaticCommunicatorMPI::getInstance()->barrier();

  printf("rank %d receive:   ",world_rank);
  for (int i=0; i<size_pproc; i++){
    std::cout<<std::setw(3)<<buffer[i];
  }
  std::cout<<std::endl<<std::flush;

  StaticCommunicatorMPI::getInstance()->barrier();

  //---------------------------------------------------
  // change buffer
  for (int i=0; i<size_pproc; i++)
    buffer[i]+=i+world_rank;

  StaticCommunicatorMPI::getInstance()->barrier();

  printf("rank %d change :   ",world_rank);
  for (int i=0; i<size_pproc; i++)
    std::cout<<std::setw(3)<<buffer[i];
  std::cout<<std::endl<<std::flush;

  StaticCommunicatorMPI::getInstance()->barrier();

  //---------------------------------------------------
  // ALL GATHER not in place
  StaticCommunicatorMPI::getInstance()->barrier();

  StaticCommunicatorMPI::getInstance()->allGather(root_buffer,buffer,size_pproc);

  StaticCommunicatorMPI::getInstance()->barrier();

  std::cout<<"rank "<<world_rank<<" allGather: ";
  for (int i=0; i<size;i++) {
    std::cout<<std::setw(3)<<root_buffer[i];
    if (root_buffer[i]!=i+i%size_pproc+i/size_pproc){
      error=true;
    }
  }
  if (error) return 1;
  std::cout<<std::endl<<std::flush;


  StaticCommunicatorMPI::getInstance()->barrier();

  if (world_rank==root) printf("allGather not-in-place success!\n") ;


  //---------------------------------------------------
  // 2d array scatter gather
  //---------------------------------------------------
  StaticCommunicatorMPI::getInstance()->barrier();
  if (world_rank==root) printf("test scatter and gather 2D \n") ;
  StaticCommunicatorMPI::getInstance()->barrier();

  fftw_complex * rootB = new fftw_complex[size];
  fftw_complex * leafB = new fftw_complex[size];

  // populate array
  if (world_rank==0) {
    std::cout<<"rank "<<world_rank<<" 2D scatter: ";
    for (int i=0; i<size;i++) {
      for (int j=0;j<2;j++) {
	rootB[i][j]=j+2*i;
	std::cout<<std::setw(3)<<rootB[i][j];
      }
    }
    std::cout<<std::endl<<std::flush;
  }
  StaticCommunicatorMPI::getInstance()->barrier();

  // scatter

  StaticCommunicatorMPI::getInstance()->scatter(rootB,leafB,size_pproc*2);

  StaticCommunicatorMPI::getInstance()->barrier();

  std::cout<<"rank "<<world_rank<<" 2D receive: ";
  for (int i=0; i<size_pproc; i++) {
    for (int j=0; j<2; j++) {
      if (leafB[i][j]!=j+2*i+size_pproc*2*world_rank)
	error=true;
      std::cout<<std::setw(3)<<leafB[i][j];
    }
  }
  std::cout<<std::endl<<std::flush;
  if (error) return 1;

  StaticCommunicatorMPI::getInstance()->barrier();

  if (world_rank==root) printf("scatter 2D success!\n");

  StaticCommunicatorMPI::getInstance()->barrier();

  //---------------------------------------------------
  // change leaf buffer
  for (int i=0; i<size_pproc; i++)
    for (int j=0; j<2; j++)
      leafB[i][j]+=i+world_rank;

  printf("rank %d change 2D: ",world_rank);
  for (int i=0; i<size_pproc; i++)
    for (int j=0; j<2;j++)
      std::cout<<std::setw(3)<<leafB[i][j];
  std::cout<<std::endl<<std::flush;
  StaticCommunicatorMPI::getInstance()->barrier();

  //---------------------------------------------------
  // Gather

  StaticCommunicatorMPI::getInstance()->gather(rootB,leafB,size_pproc*2);

  StaticCommunicatorMPI::getInstance()->barrier();

  if (world_rank==0) {
    std::cout<<"rank "<<world_rank<<" gather 2D: ";
    for (int i=0; i<size;i++) {
      for (int j=0; j<2;j++) {
	std::cout<<std::setw(3)<<rootB[i][j];
	if (rootB[i][j]!=2*i+j+i%size_pproc+i/size_pproc){
	  error=true;
	}
      }
    }

    if (error) return 1;
    std::cout<<std::endl<<std::flush;
  }

  StaticCommunicatorMPI::getInstance()->barrier();

  if (world_rank==root) printf("gather 2D success!\n");
  StaticCommunicatorMPI::getInstance()->barrier();


  //---------------------------------------------------
  // allGather
  if (world_rank==root) printf("test allGather 2D\n");

  StaticCommunicatorMPI::getInstance()->barrier();

  //---------------------------------------------------
  // change leaf buffer
  for (int i=0; i<size_pproc; i++)
    for (int j=0; j<2; j++)
      leafB[i][j]+=i+world_rank;

  printf("rank %d change 2D: ",world_rank);
  for (int i=0; i<size_pproc; i++)
    for (int j=0; j<2;j++)
      std::cout<<std::setw(3)<<leafB[i][j];
  std::cout<<std::endl<<std::flush;
  StaticCommunicatorMPI::getInstance()->barrier();


  StaticCommunicatorMPI::getInstance()->allGather(rootB,leafB,size_pproc*2);

  StaticCommunicatorMPI::getInstance()->barrier();

  if (world_rank==0) {
    std::cout<<"rank "<<world_rank<<" Gather 2D: ";
    for (int i=0; i<size;i++) {
      for (int j=0; j<2;j++) {
	std::cout<<std::setw(3)<<rootB[i][j];
	if (rootB[i][j]!=2*i+j+2*(i%size_pproc+i/size_pproc)){
	  error=true;
	}
      }
    }

    if (error) return 1;
    std::cout<<std::endl<<std::flush;
  }

  StaticCommunicatorMPI::getInstance()->barrier();

  if (world_rank==root) printf("allGather 2D success!\n") ;

  //---------------------------------------------------
  //---------------------------------------------------

  // clean up
  delete root_buffer;
  delete buffer;

  if (world_rank==root) printf("went to the end \n");
  StaticCommunicatorMPI::getInstance()->finalize();

  return 0;
}
