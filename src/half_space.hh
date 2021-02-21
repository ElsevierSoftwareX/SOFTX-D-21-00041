/**
 * @file   half_space.hh
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
#ifndef __HALF_SPACE_H__
#define __HALF_SPACE_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "material.hh"
#include "preint_kernel.hh"
#include "nodal_field.hh"
#include "fftable_nodal_field.hh"
#include "limited_history.hh"
#include "uca_mesh.hh"
#include "uca_dumper.hh"

#ifdef UCA_USE_OPENMP
#include <omp.h>
#endif /* UCA_USE_OPENMP */

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class HalfSpace {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  // side factor top=1 bot=-1

  HalfSpace(Mesh & mesh, int side_factor);

  virtual ~HalfSpace();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void computeDisplacement(bool predicting = false);
  virtual void computeStressFourierCoeff(bool predicting = false,
					 bool correcting = false);
  virtual void computeResidual(std::vector<NodalField *> &external);
  virtual void computeVelocity(bool predicting = false);

  // init convolutions
  void initConvolutions();

  // apply fft forward on displacement
  virtual void forwardFFT(bool predicting = false);
  // apply fft backward on internal
  virtual void backwardFFT();

  // for costum mesh
  virtual void gatherCostumMeshForwardFFT(std::vector<NodalField *> &scratch,
					  bool predicting = false);
  virtual void backwardFFTscatterCostumMesh();

  // for predictor-corrector implmentation
  void initPredictorCorrector();
  virtual void updateVelocity();
  virtual void correctVelocity(bool last_step);

  // dumper function: returns true if successful
  bool registerDumpFieldToDumper(const std::string & field_name, // what to dump
				 const std::string & dump_name, // how to name it
				 Dumper * const dumper);

protected:
  void computeStressFourierCoeffDynamic(bool predicting, bool correcting);
  void computeStressFourierCoeffStatic(bool predicting);

private:
  void computeDisplacement(NodalField *disp, NodalField *velo, NodalField *target);
  void computeVelocity(std::vector<NodalField *> & _velo);
  void correctVelocity(NodalField *velo_n, NodalField *velo_pc, NodalField *target);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // set a material to half space
  void setMaterial(Material *material) { this->material = material; }
  // get material of half space
  const Material &getMaterial() const { return (*this->material); }

  void setTimeStep(double time_step);

  // accessors
  FFTableNodalField * getDisp(int d, bool predicting = false) {
    return predicting ? this->disp_pc[d] : this->disp[d];
  }
  NodalField * getVelo(int d, bool predicting = false) {
    return predicting ? this->velo_pc[d] : this->velo[d];
  }

  FFTableNodalField * getInternal(int d) { return this->internal[d]; };

  NodalField * getResidual(int d) { return this->residual[d]; };

  // const accessors
  const FFTableNodalField * getDisp(int d, bool predicting = false) const {
    return predicting ? this->disp_pc[d] : this->disp[d];
  }
  const NodalField * getVelo(int d, bool predicting = false) const {
    return predicting ? this->velo_pc[d] : this->velo[d];
  }
  const FFTableNodalField * getInternal(int d) const { return internal[d]; }

  double getStableTimeStep();

  void setDynamic(bool fully_dynamic) { this->dynamic = fully_dynamic; }
  bool getDynamic() { return this->dynamic; }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  Mesh & mesh;

  // time step
  double time_step;

  // used to know in which directions the tractions pull
  int side_factor;

  // material properties
  Material * material;

  // displacement 0 x "in-plane shear" ; 1 y "normal"; 2 z "out-of-plane shear"
  std::vector<FFTableNodalField *> disp;

  // velocity
  std::vector<NodalField *> velo;

  // past values of displacement in frequency domain
  // each LimitedHistory is for a given wave number q
  std::vector<std::vector<LimitedHistory *> > U_r;
  std::vector<std::vector<LimitedHistory *> > U_i;

  // convolutions

  std::vector<PreintKernel *> H00_pi;
  std::vector<PreintKernel *> H01_pi;
  std::vector<PreintKernel *> H11_pi;
  std::vector<PreintKernel *> H22_pi;

  // tractions due to deformation
  std::vector<FFTableNodalField *> internal;

  // all acting forces
  std::vector<NodalField *> residual;

  // for predictor-corrector implmentation
  bool predictor_corrector = false;
  std::vector<FFTableNodalField *> disp_pc;
  std::vector<NodalField *> velo_pc;

  bool dynamic = true;

};

__END_UGUCA__

//#include "half_space_impl.cc"

#endif /* __HALF_SPACE_H__ */
