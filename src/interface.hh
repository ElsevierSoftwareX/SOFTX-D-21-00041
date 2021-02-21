/**
 * @file   interface.hh
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
#ifndef __INTERFACE_H__
#define __INTERFACE_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "interface_law.hh"
#include "uca_dumper.hh"
#include "half_space.hh"
#include "uca_mesh.hh"
#include <vector>

__BEGIN_UGUCA__

class Interface : public Dumper {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Interface(Mesh & mesh,
	    InterfaceLaw & law);

  virtual ~Interface();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // get (limiting) stable time step
  virtual double getStableTimeStep();

  // initiate interface model and half spaces
  virtual void init(bool velocity_initial_conditions = false);

  // initiate dump (overload to get coordinates)
  void initDump(const std::string &bname,
		const std::string &path,
                const Dumper::Format format = Dumper::Format::ASCII);

  // initiate predictor-corrector stratch memory
  virtual void initPredictorCorrector(int iterations = 1);

  // iteration of advancing one time step
  virtual void advanceTimeStep();

  // functions used during time stepping for each half-space
  virtual void computeDisplacement(bool predicting = false);
  virtual void forwardFFT(bool predicting = false);
  virtual void computeStressFourierCoeff(bool predicting = false,
					 bool correcting = false);
  virtual void backwardFFT();
  virtual void computeCohesion(bool predicting = false);
  virtual void computeResidual();
  virtual void computeVelocity(bool predicting = false);

  // for costum mesh
  virtual void gatherCostumMeshForwardFFT();
  virtual void backwardFFTscatterCostumMesh();

  // functions used during predictor-corrector time stepping
  virtual void updateVelocity();
  virtual void correctVelocity(bool last_step);

  // function that combine load and cohesionw with correct signs
  void combineLoadAndCohesion(std::vector<NodalField *> & load_and_cohesion);

  // compute force needed to close normal gap
  virtual void closingNormalGapForce(NodalField * close_force,
                                     bool predicting = false) = 0;

  // compute force needed to maintain current shear gap
  virtual void maintainShearGapForce(std::vector<NodalField *> & maintain_force) = 0;

  // compute gap in displacement
  virtual void computeGap(std::vector<NodalField *> & gap,
                          bool predicting = false) = 0;

  // compute gap relative velocity
  virtual void computeGapVelocity(std::vector<NodalField *> & gap_velo,
                                  bool predicting = false) = 0;

  // compute norm of provided field
  virtual void computeNorm(const std::vector<NodalField *> & field,
			   NodalField & norm, bool shear_only = false);

  // multiply fields elementwise with scalar
  virtual void multiplyFieldByScalar(std::vector<NodalField *> & field,
				     const NodalField & scalar,
				     bool shear_only = false);

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  NodalField * getShearLoad()  { return this->load[0]; };
  NodalField * getNormalLoad() { return this->load[1]; };
  NodalField * getLoad(int d) { return this->load[d]; };

  NodalField * getCohesion(int d) { return this->cohesion[d]; };

  std::vector<NodalField *> & getCohesion()    { return this->cohesion; };
  std::vector<NodalField *> & getBufferField() { return this->scratch_field; };

  virtual HalfSpace & getTop() = 0;
  virtual HalfSpace & getBot() = 0;
  double getTimeStep() const { return this->time_step; };

  virtual void setTimeStep(double time_step);
  double getTimeStep() {return this->time_step; };
  void setDynamic(bool fully_dynamic);
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  Mesh & mesh;

  // time step
  double time_step;

  // {top, bot}
  std::vector<HalfSpace *> half_space;

  // external loading
  std::vector<NodalField *> load;

  // interface forces (e.g., cohesion)
  std::vector<NodalField *> cohesion;

  std::vector<NodalField *> scratch_field;

  // predictor-corrector
  int nb_pc = 0;

  // interface law: cohesive law and contact/friction law
  InterfaceLaw & law;
};

__END_UGUCA__

//#include "interface_impl.cc"

#endif /* __INTERFACE_H__ */
