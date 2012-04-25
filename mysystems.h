
// Local includes
#include "bioheatsystem.h"
#include "femparameters.h"
#include "equation_systems.h"

// LibMesh includes
#include "getpot.h"

FEMSystem &build_system(EquationSystems &es, GetPot & /* infile */, FEMParameters &param)
{
  BioheatSystem &bsystem = es.add_system<BioheatSystem> ("BioheatSystem");

  // Use the prescribed FE type
  bsystem.fe_family() = param.fe_family[0];
  bsystem.fe_order() = param.fe_order[0];

  // Use the prescribed material fields
  const unsigned int n_funcs = 
    param.other_interior_functions.size();
  if (n_funcs != 1)
    {
      libMesh::out << "Exactly one interior function should be specified: the liquid fraction" << std::endl;
      libMesh::out << "Error: " << n_funcs << " interior functions were found." << std::endl;
      exit(1);
    }
  bsystem.n_L_values = param.other_interior_functions[0];

  // Use the specified boundary conditions
  bsystem.dirichlet_values = param.dirichlet_conditions;
  bsystem.neumann_values = param.neumann_conditions;
  bsystem.dirichlet_variables = param.dirichlet_condition_variables;
  bsystem.neumann_variables = param.neumann_condition_variables;

  // Use the numerical or analytic jacobians?
  bsystem.analytic_jacobians() = param.analytic_jacobians;

  // Prescribe other parameters in a separate config file
  bsystem.system_config_file() = param.system_config_file;

  return bsystem;
}


// Prototype for calculation of the exact solution.  Useful
// for verification
Number exact_value(const Point&,       // xyz location
                   const Parameters&,  // EquationSystem parameters
                   const std::string&, // sys_name
                   const std::string&) // unknown_name
{libmesh_error(); return 0;}

// Prototype for calculation of the gradient of the exact solution.
Gradient exact_grad(const Point&,       // xyz location
                    const Parameters&,  // EquationSystems parameters
                    const std::string&, // sys_name
                    const std::string&) // unknown_name
{libmesh_error(); return 0;}

