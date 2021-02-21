/**
 * @file   material.hh
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
#ifndef __MATERIAL_H__
#define __MATERIAL_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "kernel_collection.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */
class Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Material(double E, double nu, double rho, bool pstress=false);
  virtual ~Material() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void readPrecomputedKernels(const std::string & path = global_kernel_path) {
    this->kernels.readPrecomputedKernels(this->nu, this->pstress, path);
  };

private:

  void computeShearModulus();
  void computeFirstLame();
  void computeCp();
  void computeCs();
  void computeCr();
  double computeCr_bisect(double a, double b);
  double computeCr_fcn(double v);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  double getYoungsModulus() const { return this->E; };
  double getShearModulus() const { return this->mu; };
  double getPoissonRatio() const { return this->nu; };
  double getDensity() const { return this->rho; };
  double getCp() const { return this->cp; };
  double getCs() const { return this->cs; };
  double getCr() const { return this->cr; };
  bool getPStress() const { return this->pstress; };

  // access kernels
  Kernel * getH00() { return this->kernels.getH00(); }
  Kernel * getH01() { return this->kernels.getH01(); }
  Kernel * getH11() { return this->kernels.getH11(); }
  Kernel * getH22() { return this->kernels.getH22(); }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // PRIMARY
  // young's modulus
  double E;
  // poisson's ratio
  double nu;
  // density
  double rho;
  // plane stress
  bool pstress;

  // kernel collection
  KernelCollection kernels;

  // SECONDARY
  // first lame
  double lambda;
  // shear modulus
  double mu;
  // p-wave speed
  double cp;
  // s-wave speed
  double cs;
  // Rayleigh wave speed
  double cr;
};

__END_UGUCA__

//#include "material_impl.cc"

#endif /* __MATERIAL_H__ */
