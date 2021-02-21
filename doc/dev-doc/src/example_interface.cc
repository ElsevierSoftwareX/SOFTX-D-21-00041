#include "example_interface.hh"
#include "static_communicator_mpi.hh"

__BEGIN_UGUCA__

/* -------------------------------------------------------------------------- */

ExampleInterface::ExampleInterface(Mesh & mesh,
				   Material & top_material,
				   InterfaceLaw & law) :
  Interface(mesh, law),
  top(mesh,  1),
{

  this->half_space.resize(1);
  this->half_space[0] = &this->top;

  this->top.setMaterial(&top_material);
}

/* -------------------------------------------------------------------------- */
ExampleInterface::~ExampleInterface() {}

/* -------------------------------------------------------------------------- */
void ExampleInterface::closingNormalGapForce(NodalField *close_force,
					     bool predicting) {
  // YOUR CODE COMES HERE
  // YOU NEED TO FILL THE CLOSE_FORCE FIELD
}

/* -------------------------------------------------------------------------- */
void ExampleInterface::maintainShearGapForce(std::vector<NodalField *> &maintain_force) {
  // YOUR CODE COMES HERE
  // YOU NEED TO FILL THE MAINTAIN_FORCE FIELD
}

/* -------------------------------------------------------------------------- */
void ExampleInterface::computeGap(std::vector<NodalField *> & gap,
                                bool predicting) {
  // YOUR CODE COMES HERE
  // YOU NEED TO FILL THE GAP FIELD
}

/* -------------------------------------------------------------------------- */
void ExampleInterface::computeGapVelocity(std::vector<NodalField *> & gap_velo,
                                        bool predicting) {
  // YOUR CODE COMES HERE
  // YOU NEED TO FILL THE GAP_VELO FIELD
}

/* -------------------------------------------------------------------------- */
void ExampleInterface::registerDumpField(const std::string &field_name) {

  int world_rank = StaticCommunicatorMPI::getInstance()->whoAmI();

  if (world_rank != 0) return;

  int d = std::atoi(&field_name[field_name.length() - 1]);

  if (d >= this->mesh.getDim())
    throw std::runtime_error("Field "+field_name
			     +" cannot be dumped, to high dimension");
  
  bool registered = false;
  // field_name starts with "top"
  if (field_name.rfind("top", 0) == 0) {
    // cut away "top_" from field_name and give interface as dumper
    registered = this->top.registerDumpFieldToDumper(field_name.substr(4),
						     field_name,
						     this);
  }

  // DO SAME WITH BOT IF NEEDED
  
  if (!registered) {
    // YOUR NEW FIELDS COME HERE
    if (field_name == "example_field" + std::to_string(d))
      this->registerForDump(field_name, this->example_field[d]);
    else
      Interface::registerDumpField(field_name);
  }
}

__END_UGUCA__
