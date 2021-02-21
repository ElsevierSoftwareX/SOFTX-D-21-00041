/**
 * @file   half_space.cc
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
#include "half_space.hh"
#include "static_communicator_mpi.hh"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <algorithm>

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */


HalfSpace::HalfSpace(Mesh & mesh,
		     int side_factor) :
  mesh(mesh),
  time_step(0.),
  side_factor(side_factor) {

  this->disp.resize(this->mesh.getDim());
  this->velo.resize(this->mesh.getDim());
  this->U_r.resize(this->mesh.getDim());
  this->U_i.resize(this->mesh.getDim());
  this->internal.resize(this->mesh.getDim());
  this->residual.resize(this->mesh.getDim());

  for (int d=0; d<this->mesh.getDim(); ++d) {
    this->disp[d]      = new FFTableNodalField(this->mesh);
    this->internal[d]  = new FFTableNodalField(this->mesh);
    this->velo[d]      = new        NodalField(this->mesh.getNbNodes());
    this->residual[d]  = new        NodalField(this->mesh.getNbNodes());
  }

  this->H00_pi.resize(this->mesh.getNbFFT());
  this->H01_pi.resize(this->mesh.getNbFFT());
  this->H11_pi.resize(this->mesh.getNbFFT());
  if (this->mesh.getDim()==3)
    this->H22_pi.resize(this->mesh.getNbFFT());

  for (int d=0; d<this->mesh.getDim(); ++d) {
    this->U_r[d].resize(this->mesh.getNbFFT());
    this->U_i[d].resize(this->mesh.getNbFFT());
  }

#ifdef UCA_VERBOSE
  printf("%dD HalfSpace constructed\n"
	 "nb_nds_local  %d, nb_FFT_local %d\n",
	 this->mesh.getDim(),
	 this->mesh.getNbNodes(),this->mesh.getNbFFT());
#endif /* UCA_VERBOSE */
}

/* -------------------------------------------------------------------------- */
HalfSpace::~HalfSpace() {

  for (int d=0; d<this->mesh.getDim(); ++d) {
    delete this->disp[d];
    delete this->velo[d];

    if (predictor_corrector) {
      delete this->disp_pc[d];
      delete this->velo_pc[d];
    }

    delete this->internal[d];
    delete this->residual[d];

    for (int j=0; j<this->mesh.getNbFFT(); ++j) {
      delete this->U_r[d][j];
      delete this->U_i[d][j];
    }
  }

  for (int j=0; j<this->mesh.getNbFFT(); ++j) {
    delete this->H00_pi[j];
    delete this->H01_pi[j];
    delete this->H11_pi[j];
    if (this->mesh.getDim()==3)
      delete this->H22_pi[j];
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::initPredictorCorrector() {
  this->predictor_corrector = true;
  this->disp_pc.resize(this->mesh.getDim());
  this->velo_pc.resize(this->mesh.getDim());
  for (int d = 0; d < this->mesh.getDim(); ++d) {
    this->disp_pc[d] = new FFTableNodalField(this->mesh);
    this->velo_pc[d] = new NodalField(this->mesh.getNbNodes());
  }
}

/* -------------------------------------------------------------------------- */
double HalfSpace::getStableTimeStep() {
  double delta_x = this->mesh.getDeltaX();
  double delta_z = this->mesh.getDeltaZ();

  if (this->mesh.getDim()==2)
    delta_z = 1e100;
  return std::min(delta_x,delta_z) / this->material->getCs();
}

/* -------------------------------------------------------------------------- */
void HalfSpace::setTimeStep(double time_step) {
  if (this->dynamic) {
    if (((time_step/this->getStableTimeStep()>0.35) && (this->mesh.getDim()==3)) ||
        ((time_step/this->getStableTimeStep()>0.4) && (this->mesh.getDim()==2)))
      throw std::runtime_error("Error: time_step_factor is too large: (<0.4 for 2D and <0.35 for 3D)\n");
  }

  this->time_step = time_step;
}

/* -------------------------------------------------------------------------- */
void HalfSpace::initConvolutions() {
  if (!this->dynamic) return;
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();

  int total_work=0;

  std::vector<NodalField *> wave_numbers = this->mesh.getWaveNumbers();

  // history for q1 is longest q = j*q1

  int j=0; // mode 0 discarted during job assignment
  if (world_size==1) j=1; // discard mode 0
  for (; j<this->mesh.getNbFFT(); ++j) { //parallel loop

    this->H00_pi[j] = new PreintKernel(this->material->getH00());
    this->H01_pi[j] = new PreintKernel(this->material->getH01());
    this->H11_pi[j] = new PreintKernel(this->material->getH11());
    if (this->mesh.getDim()==3)
      this->H22_pi[j] = new PreintKernel(this->material->getH22());

    double qq = 0.0;
    for (int d=0; d<this->mesh.getDim();d+=2)
      qq +=(*wave_numbers[d])(j)*(*wave_numbers[d])(j);

    double qj_cs = std::sqrt(qq) * this->material->getCs();

    this->H00_pi[j]->preintegrate(qj_cs, this->time_step);
    this->H01_pi[j]->preintegrate(qj_cs, this->time_step);
    this->H11_pi[j]->preintegrate(qj_cs, this->time_step);
    if (this->mesh.getDim()==3)
      this->H22_pi[j]->preintegrate(qj_cs, this->time_step);

    std::vector<int> nb_hist={0,0,0};

    nb_hist[0] = std::max(this->H00_pi[j]->getSize(),
			  this->H01_pi[j]->getSize());
    nb_hist[1] = std::max(this->H11_pi[j]->getSize(),
			  this->H01_pi[j]->getSize());
    if (this->mesh.getDim()==3) {
      nb_hist[0] = std::max(std::max(this->H00_pi[j]->getSize(),
				     this->H01_pi[j]->getSize()),
			    this->H22_pi[j]->getSize());
      nb_hist[2] = nb_hist[0];
    }

    for (int d=0; d<this->mesh.getDim(); ++d) {
      this->U_r[d][j] = new LimitedHistory(nb_hist[d]);
      this->U_i[d][j] = new LimitedHistory(nb_hist[d]);
    }

    for (int d=0; d<this->mesh.getDim(); ++d)
      total_work += nb_hist[d];
  }
#ifdef UCA_VERBOSE
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  printf("Rank %d has total work %d \n",world_rank,total_work);
#endif /* UCA_VERBOSE */
}


/* -------------------------------------------------------------------------- */
void HalfSpace::computeDisplacement(bool predicting) {
  for (int d = 0; d < this->mesh.getDim(); ++d) {
    this->computeDisplacement(this->disp[d],
                              predicting ? this->velo_pc[d] : this->velo[d],
                              predicting ? this->disp_pc[d] : this->disp[d]);
  }
}

/* -------------------------------------------------------------------------- */
// u_i+1 = u_i + dt * v_i
void HalfSpace::computeDisplacement(NodalField * _disp,
				    NodalField * _velo, NodalField * target) {

  double * disp_p = _disp->storage();
  double * velo_p = _velo->storage();
  double * target_p = target->storage();

  for (int n=0; n<this->mesh.getNbNodes(); ++n) {
    target_p[n] = disp_p[n] + velo_p[n] * this->time_step;
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::forwardFFT(bool predicting) {
  for (int d = 0; d < this->mesh.getDim(); ++d) {
    // if parallel automatically exectued by rank 0
    (predicting ? this->disp_pc[d] : this->disp[d])->forwardFFT();
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::backwardFFT() {
  for (int d = 0; d < this->mesh.getDim(); ++d) {
    // if parallel automatically exectued by rank 0
    this->internal[d]->backwardFFT();
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::gatherCostumMeshForwardFFT(std::vector<NodalField *> & scratch,
					   bool predicting) {
  std::vector<FFTableNodalField *> & _disp =
    predicting ? this->disp_pc : this->disp;

  for (int d = 0; d < this->mesh.getDim(); ++d) {
    // copy local disp
    double * disp_p = _disp[d]->storage();
    double * disp_copy_p = scratch[d]->storage();

    for (int n=0; n<this->mesh.getNbNodes();++n)
      disp_copy_p[n]=disp_p[n];
  }

  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  if (world_size>this->mesh.getDim()) {
    // for trivially parallel fftw rank i does dimension i for i<dim
    for (int d = 0; d < this->mesh.getDim(); ++d) {
      this->mesh.gatherAndSortCostumNodes(_disp[d]->storage(),d);
    }

    if (world_rank < this->mesh.getDim())
      this->disp[world_rank]->forwardFFT(world_rank); // exectued by world_rank
  }
  else {
    for (int d = 0; d < this->mesh.getDim(); ++d) {
      this->mesh.gatherAndSortCostumNodes(_disp[d]->storage());
      _disp[d]->forwardFFT(); // exectued by world_rank
    }
  }

  for (int d = 0; d < this->mesh.getDim(); ++d) {
    // copy local disp back
    double * disp_p = _disp[d]->storage();
    double * disp_copy_p = scratch[d]->storage();

    for (int n=0; n<this->mesh.getNbNodes();++n)
      disp_p[n]=disp_copy_p[n];
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::backwardFFTscatterCostumMesh() {
#ifdef UCA_USE_MPI
  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();
  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  if (world_size>this->mesh.getDim()) {
    // for trivially parallel fftw rank i does dimension i for i<dim
    if (world_rank < this->mesh.getDim())
      this->internal[world_rank]->backwardFFT(world_rank); // exectued by world_rank
    for (int d = 0; d < this->mesh.getDim(); ++d)
      this->mesh.sortAndScatterCostumNodes(this->internal[d]->storage(),d);
  }
  else {
    for (int d = 0; d < this->mesh.getDim(); ++d) {
      this->internal[d]->backwardFFT(); // exectued by world_rank
      this->mesh.sortAndScatterCostumNodes(this->internal[d]->storage());
    }
  }
#else
  for (int d = 0; d < this->mesh.getDim(); ++d) {
    // if parallel automatically exectued by rank 0
    this->internal[d]->backwardFFT();
    this->mesh.sortAndScatterCostumNodes(this->internal[d]->storage());
  }
#endif
}

/* -------------------------------------------------------------------------- */
void HalfSpace::computeStressFourierCoeff(bool predicting, bool correcting) {
  if (this->dynamic) {
    this->computeStressFourierCoeffDynamic(predicting, correcting);
  } else {
    this->computeStressFourierCoeffStatic(predicting);
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::computeStressFourierCoeffDynamic(bool predicting, bool correcting) {
  std::vector<FFTableNodalField *> &_disp =
      predicting ? this->disp_pc : this->disp;

  int world_size = StaticCommunicatorMPI::getInstance()->getNbProc();

  double mu = this->material->getShearModulus();
  double eta = this->material->getCp() / this->material->getCs();

  std::vector<NodalField *> wave_numbers = this->mesh.getWaveNumbers();

  // imaginary number i
  std::complex<double> imag = {0., 1.};

  // access to fourier coefficients of stresses
  std::vector<fftw_complex *> internal_fd;
  internal_fd.resize(this->mesh.getDim());

  for (int d = 0; d < this->mesh.getDim(); ++d) {
    internal_fd[d] = this->internal[d]->fd_storage();
    // scatter displacements for parallel convolution
    if (this->mesh.isCostum() && world_size>this->mesh.getDim())
      // for trivially parallel fftw rank i does dimension i for i<dim
      this->mesh.sortAndScatterFFTModes(_disp[d]->fd_storage(), d);
    else
      this->mesh.sortAndScatterFFTModes(_disp[d]->fd_storage());
  }

  int j=0;
  if (world_size==1) j=1; //discard mode 0

  for (; j<this->mesh.getNbFFT(); ++j) { // parallel loop over km modes

    std::vector<std::complex<double>> U;
    U.resize(this->mesh.getDim());

    for (int d = 0; d < this->mesh.getDim(); ++d) {
      U[d] = {_disp[d]->fd(j)[0], _disp[d]->fd(j)[1]};

      // store current displacement in history
      if (correcting) {
        this->U_r[d][j]->changeCurrentValue(std::real(U[d]));
        this->U_i[d][j]->changeCurrentValue(std::imag(U[d]));
      } else {
        this->U_r[d][j]->addCurrentValue(std::real(U[d]));
        this->U_i[d][j]->addCurrentValue(std::imag(U[d]));
      }
    }

    int nb_conv = (int)std::pow(2.0, this->mesh.getDim());

    std::vector<std::complex<double>> conv;
    conv.resize(nb_conv);

    std::vector<PreintKernel *> krnl = {
      this->H00_pi[j],
      this->H01_pi[j],
      this->H01_pi[j],
      this->H11_pi[j]};

    if (this->mesh.getDim() == 3) {  // add 3d part
      std::vector<PreintKernel *> krnl3d = {
	this->H00_pi[j],
	this->H01_pi[j],
	this->H22_pi[j],
	this->H22_pi[j]};

      krnl.insert(krnl.end(), krnl3d.begin(), krnl3d.end());
    }

    std::vector<LimitedHistory *> U_real = {
      this->U_r[0][j],
      this->U_r[0][j],
      this->U_r[1][j],
      this->U_r[1][j]};

    if (this->mesh.getDim() == 3) {  // add 3d part
      std::vector<LimitedHistory *> U_real3d = {
	this->U_r[2][j],
	this->U_r[2][j],
	this->U_r[0][j],
	this->U_r[2][j]};

      U_real.insert(U_real.end(), U_real3d.begin(), U_real3d.end());
    }

    std::vector<LimitedHistory *> U_imag = {
      this->U_i[0][j],
      this->U_i[0][j],
      this->U_i[1][j],
      this->U_i[1][j]};

    if (this->mesh.getDim() == 3) {  // add 3d part
      std::vector<LimitedHistory *> U_imag3d = {
	this->U_i[2][j],
	this->U_i[2][j],
	this->U_i[0][j],
	this->U_i[2][j]};

      U_imag.insert(U_imag.end(), U_imag3d.begin(), U_imag3d.end());
    }

#ifdef UCA_USE_OPENMP
#pragma omp parallel for
#endif
    for (int n = 0; n < nb_conv; ++n) {
      conv[n] = krnl[n]->convolve(U_real[n], U_imag[n]);
    }
    // convolutions for both 2d and 3d
    std::complex<double> conv_H00_U0_j = conv[0];
    std::complex<double> conv_H01_U0_j = conv[1];
    std::complex<double> conv_H01_U1_j = conv[2];
    std::complex<double> conv_H11_U1_j = conv[3];

    std::vector<std::complex<double>> F;
    F.resize(this->mesh.getDim());
    if (this->mesh.getDim() == 2) {
      double q = (*wave_numbers[0])(j);

      // - mu * q * int(H00,U0)
      F[0] = -this->side_factor * mu * q * conv_H00_U0_j;

      // + i * (2 - eta) * mu * q * U1
      F[0] += mu * q * (2 - eta) * (imag * U[1]);

      // + i * mu * q * int(H01, U1)
      F[0] += mu * q * imag * conv_H01_U1_j;

      // - mu * q * int(H11, U1)
      F[1] = -this->side_factor * mu * q * conv_H11_U1_j;

      // - i * (2 - etq) * mu * q * U0
      F[1] -= mu * q * (2 - eta) * (imag * U[0]);

      // - i * mu * q * int(H01, U0)
      F[1] -= mu * q * imag * conv_H01_U0_j;
    } else {
      // convolutions for 3d only
      std::complex<double> conv_H00_U2_j = conv[4];
      std::complex<double> conv_H01_U2_j = conv[5];
      std::complex<double> conv_H22_U0_j = conv[6];
      std::complex<double> conv_H22_U2_j = conv[7];
      // q = {k,m} wave number in x,y direction

      double k = (*wave_numbers[0])(j);
      double m = (*wave_numbers[2])(j);

      double q = std::sqrt(k * k + m * m);

      F[0] = imag * mu * (2 - eta) * k * U[1];

      F[0] += imag * mu * k * conv_H01_U1_j;

      F[0] -= this->side_factor * mu *
	((conv_H00_U0_j * ((k * k) / q) + conv_H00_U2_j * ((k * m) / q)) +
	 (conv_H22_U0_j * ((m * m) / q) - conv_H22_U2_j * ((k * m) / q)));

      F[1] = -imag * mu * (2 - eta) * (k * U[0] + m * U[2]);

      F[1] -= mu * imag * (conv_H01_U0_j * k + conv_H01_U2_j * m);

      F[1] -= this->side_factor * mu * q * conv_H11_U1_j;

      F[2] = imag * mu * (2 - eta) * m * U[1];

      F[2] += imag * mu * m * conv_H01_U1_j;

      F[2] -= this->side_factor * mu *
	((conv_H00_U0_j * ((k * m) / q) + conv_H00_U2_j * ((m * m) / q)) +
	 (-conv_H22_U0_j * ((k * m) / q) + conv_H22_U2_j * ((k * k) / q)));
    }
    // set values to internal force
    for (int d = 0; d < this->mesh.getDim(); ++d) {
      internal_fd[d][j][0] = std::real(F[d]);  // real part
      internal_fd[d][j][1] = std::imag(F[d]);  // imag part
    }
  }

  for (int d = 0; d < this->mesh.getDim(); ++d) {
    // gather internal force after parallel convolution
    if (this->mesh.isCostum() && world_size>this->mesh.getDim())
        //for trivially parallel fftw rank i does dimension i for i<dim
      this->mesh.gatherAndSortFFTModes(internal_fd[d],d);
    else
      this->mesh.gatherAndSortFFTModes(internal_fd[d]);

    // set for mode zero fourier coefficients to zero
    internal_fd[d][0][0] = 0.;  // real part
    internal_fd[d][0][1] = 0.;  // imag part
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::computeStressFourierCoeffStatic(bool /*predicting*/) {
  throw std::runtime_error(
      "HalfSpace::computeStressFourierCoeffStatic has not been implemented.");
}

/* -------------------------------------------------------------------------- */
// residual = (internal + external) * side_factor
void HalfSpace::computeResidual(std::vector<NodalField *> & external) {
  for (int d = 0; d < this->mesh.getDim(); ++d) {
    double *int_p = this->internal[d]->storage();
    double *ext_p = external[d]->storage();
    double *res_p = this->residual[d]->storage();

    for (int n = 0; n < this->mesh.getNbNodes(); ++n) {
      res_p[n] = this->side_factor * (int_p[n] + ext_p[n]);
    }
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::computeVelocity(bool predicting) {
  this->computeVelocity(predicting ? this->velo_pc : this->velo);
}

/* -------------------------------------------------------------------------- */
void HalfSpace::updateVelocity() {
  for (int d = 0; d < this->mesh.getDim(); ++d) {
    double *velo_p = this->velo[d]->storage();
    double *velo_pc_p = this->velo_pc[d]->storage();
    for (int n = 0; n < this->mesh.getNbNodes(); ++n) {
      velo_pc_p[n] = velo_p[n];
    }
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::correctVelocity(bool last_step) {
  for (int d = 0; d < this->mesh.getDim(); ++d) {
    this->correctVelocity(this->velo[d], this->velo_pc[d],
                          last_step ? this->velo[d] : this->velo_pc[d]);
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::correctVelocity(NodalField * velo_n,
         NodalField * velo_pc, NodalField * target) {

  double * velo_n_p = velo_n->storage();
  double * velo_pc_p = velo_pc->storage();
  double * target_p = target->storage();

  for (int n = 0; n < this->mesh.getNbNodes(); ++n) {
    target_p[n] = 0.5 * (velo_n_p[n] + velo_pc_p[n]);
  }
}

/* -------------------------------------------------------------------------- */
// velocity = cs / mu       * residual (for in-plane shear components)
// velocity = cs / mu / eta * residual (for normal component)
// velocity = cs / mu       * residual (for out-of-plane shear components)
void HalfSpace::computeVelocity(std::vector<NodalField *> & _velo) {
  double mu = this->material->getShearModulus();
  double Cs = this->material->getCs();
  double Cp = this->material->getCp();
  std::vector<double> eta = {1.0, Cp / Cs, 1.0};

  for (int d=0; d < this->mesh.getDim(); ++d) {
    double * velo_p = _velo[d]->storage();
    double * res_p = this->residual[d]->storage();
    double eta_d = eta[d];

    for (int n=0; n<this->mesh.getNbNodes(); ++n)
      velo_p[n] = Cs / mu / eta_d * res_p[n];
  }
}

/* -------------------------------------------------------------------------- */
bool HalfSpace::registerDumpFieldToDumper(const std::string & field_name,
					  const std::string & dump_name,
					  Dumper * const dumper) {

  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();
  if (world_rank != 0) return true;

  int d = std::atoi(&field_name[field_name.length() - 1]);

  if (d >= this->mesh.getDim())
    throw std::runtime_error("Field "+field_name
			     +" cannot be dumped, too high dimension");
  
  // disp
  if (field_name == "disp_" + std::to_string(d)) {
    dumper->registerForDump(dump_name, this->disp[d]);
    return true;
  }
  // velo
  else if (field_name == "velo_" + std::to_string(d)) {
    dumper->registerForDump(dump_name, this->velo[d]);
    return true;
  }
  // residual
  else if (field_name == "residual_" + std::to_string(d)) {
    dumper->registerForDump(dump_name, this->residual[d]);
    return true;
  }
  // internal
  else if (field_name == "internal_" + std::to_string(d)) {
    dumper->registerForDump(dump_name, this->internal[d]);
    return true;
  }

  return false;
}

__END_UGUCA__
