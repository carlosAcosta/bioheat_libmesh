
#ifndef __fem_parameters_h__
#define __fem_parameters_h__

#include "libmesh_common.h"
#include "dof_map.h"
#include "enum_norm_type.h"
#include "function_base.h"
#include "getpot.h"
#include "id_types.h"
#include "periodic_boundaries.h"

#include <limits>
#include <map>
#include <string>
#include <vector>



class FEMParameters
{
public:
    FEMParameters();

    ~FEMParameters();

    void read(GetPot &input);

    // Parameters applicable to entire EquationSystems:

    unsigned int initial_timestep, n_timesteps;
    bool transient;
    unsigned int deltat_reductions;
    std::string timesolver_core;
    Real end_time, deltat, timesolver_theta,
	 timesolver_maxgrowth, timesolver_tolerance,
	 timesolver_upper_tolerance, steadystate_tolerance;
    std::vector<FEMNormType> timesolver_norm;

    //   Mesh generation

    unsigned int dimension;
    std::string domaintype, domainfile, elementtype;
    Real elementorder;
    Real domain_xmin, domain_ymin, domain_zmin;
    Real domain_edge_width, domain_edge_length, domain_edge_height;
    unsigned int coarsegridx, coarsegridy, coarsegridz;
    unsigned int coarserefinements, extrarefinements;

    //   Mesh refinement

    unsigned int nelem_target;
    Real global_tolerance;
    Real refine_fraction, coarsen_fraction, coarsen_threshold;
    unsigned int max_adaptivesteps;
    unsigned int initial_adaptivesteps;

    //   Output

    unsigned int write_interval;
    bool write_gmv_error, write_tecplot_error, output_xda, output_xdr,
	 output_bz2, output_gz, output_gmv, output_tecplot;

    // Types of Systems to create

    std::vector<std::string> system_types;

    // Parameters applicable to each system:

    //   Boundary and initial conditions

    std::vector<PeriodicBoundary> periodic_boundaries;

    std::map<subdomain_id_type, FunctionBase<Number> *> initial_conditions;
    std::map<boundary_id_type, FunctionBase<Number> *>  dirichlet_conditions,
                                                        neumann_conditions;
    std::map<boundary_id_type, std::vector<unsigned int> > dirichlet_condition_variables,
                                                           neumann_condition_variables;
    std::map<int, std::map<subdomain_id_type, FunctionBase<Number> *> >
      other_interior_functions;
    std::map<int, std::map<boundary_id_type, FunctionBase<Number> *> >
      other_boundary_functions;

    //   Execution type

    bool run_simulation, run_postprocess;

    //   Approximation type

    std::vector<std::string> fe_family;
    std::vector<unsigned int> fe_order;
    int extra_quadrature_order;

    //   Assembly options

    bool analytic_jacobians;
    Real verify_analytic_jacobians;
    Real numerical_jacobian_h;

    bool print_solution_norms, print_solutions,
         print_residual_norms, print_residuals,
         print_jacobian_norms, print_jacobians,
         print_element_jacobians;

    //   Solver options

    bool use_petsc_snes;
    bool time_solver_quiet, solver_quiet, solver_verbose, require_residual_reduction;
    Real min_step_length;
    unsigned int max_linear_iterations, max_nonlinear_iterations;
    Real relative_step_tolerance, relative_residual_tolerance,
	 initial_linear_tolerance, minimum_linear_tolerance,
	 linear_tolerance_multiplier;

    // Initialization

    unsigned int initial_sobolev_order;
    unsigned int initial_extra_quadrature;

    //   Error indicators

    std::string indicator_type;
    unsigned int sobolev_order;

    //   System-specific parameters file

    std::string system_config_file;
};

#endif // __fem_parameters_h__
