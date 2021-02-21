#ifndef __EXAMPLE_INTERFACE_H__
#define __EXAMPLE_INTERFACE_H__
/* -------------------------------------------------------------------------- */
#include "uca_common.hh"
#include "half_space.hh"
#include "interface_law.hh"
#include "interface.hh"

__BEGIN_UGUCA__

class ExampleInterface : public Interface {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ExampleInterface(Mesh & mesh,
		   Material & top_material,
		   InterfaceLaw & law);

  virtual ~ExampleInterface();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // compute force needed to close normal gap
  virtual void closingNormalGapForce(NodalField * close_force,
				     bool predicting = false);

  // compute force needed to maintain current shear gap
  virtual void maintainShearGapForce(std::vector<NodalField *> & maintain_force);

  // compute gap in displacement
  virtual void computeGap(std::vector<NodalField *> & gap,
			  bool predicting = false);

  // compute gap relative velocity
  virtual void computeGapVelocity(std::vector<NodalField *> & gap_velo,
				  bool predicting = false);

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  virtual HalfSpace & getTop() { return this->top; };
  virtual HalfSpace & getBot() {
#ifdef UCA_VERBOSE
    std::cout << "Warning: UnimatShearInterface::getBot() returns the same "
	      << "HalfSpace as UnimatShearInterface::getTop()" << std::endl;
#endif /* UCA_VERBOSE */
    return this->top; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  // half spaces
  HalfSpace top;
  
  // YOU CAN ADD NEW FIELDS HERE
  std::vector<NodalField *> & example_field
};

__END_UGUCA__

//#include "example_interface_impl.cc"

#endif /* __EXAMPLE_INTERFACE_H__ */
