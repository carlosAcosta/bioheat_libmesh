#include "initial.h"
#include "femparameters.h"
#include "getpot.h"

FEMParameters *param = NULL;

void read_initial_parameters()
{
  GetPot infile("general.in");
  param = new FEMParameters();
  param->read(infile);
}

void finish_initialization()
{
  delete param;
}

// Initial conditions
Number initial_value(const Point& p,
                     const Parameters&,
                     const std::string&,
                     const std::string&)
{
  // FIXME: we're ignoring per-subdomain initial conditions for now

  libmesh_assert(param->initial_conditions.count(0));
  return (*(param->initial_conditions[0]))(p);
}

Gradient initial_grad(const Point&,
                      const Parameters&,
                      const std::string&,
                      const std::string&)
{
  libmesh_not_implemented();
  return Gradient(0.);
}
