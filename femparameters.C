
#include "femparameters.h"

#include "factory.h"
#include "parsed_function.h"
#include "zero_function.h"

#define GETPOT_INPUT(A) { A = input(#A, A);\
  const std::string stringval = input(#A, std::string());\
  variable_names.push_back(std::string(#A "=") + stringval); }
#define GETPOT_INT_INPUT(A) { A = input(#A, (int)A);\
  const std::string stringval = input(#A, std::string());\
  variable_names.push_back(std::string(#A "=") + stringval); }

#define GETPOT_FUNCTION_INPUT(A) {\
  const std::string type  = input(#A "_type", "zero");\
  const std::string value = input(#A "_value", "");\
  A = new_function_base(type, value); }
#define GETPOT_REGISTER(A) { \
  std::string stringval = input(#A, std::string());\
  variable_names.push_back(std::string(#A "=") + stringval); }

FEMParameters::FEMParameters() :
      initial_timestep(0), n_timesteps(100),
      transient(true),
      deltat_reductions(0),
      timesolver_core("euler"),
      end_time(std::numeric_limits<Real>::max()),
      deltat(0.0001), timesolver_theta(0.5),
      timesolver_maxgrowth(0.), timesolver_tolerance(0.),
      timesolver_upper_tolerance(0.),
      steadystate_tolerance(0.),
      timesolver_norm(0,L2),

      dimension(2), 
      domaintype("square"), domainfile("mesh.xda"), elementtype("quad"),
      elementorder(2),
      domain_xmin(0.0), domain_ymin(0.0), domain_zmin(0.0),
      domain_edge_width(1.0), domain_edge_length(1.0), domain_edge_height(1.0),
      coarsegridx(1), coarsegridy(1), coarsegridz(1),
      coarserefinements(4), extrarefinements(0),

      nelem_target(8000), global_tolerance(0.0),
      refine_fraction(0.3), coarsen_fraction(0.3), coarsen_threshold(10),
      max_adaptivesteps(1),
      initial_adaptivesteps(0), 

      write_interval(10),
      write_gmv_error(false), write_tecplot_error(false),
      output_xda(false), output_xdr(false),
      output_bz2(true), output_gz(true),
      output_gmv(false), output_tecplot(false),

      system_types(0),

      periodic_boundaries(0),

      run_simulation(true), run_postprocess(false),

      fe_family(1, "LAGRANGE"), fe_order(1, 1),
      extra_quadrature_order(0),

      analytic_jacobians(true), verify_analytic_jacobians(0.0),
      numerical_jacobian_h(TOLERANCE),
      print_solution_norms(false), print_solutions(false),
      print_residual_norms(false), print_residuals(false),
      print_jacobian_norms(false), print_jacobians(false),
      print_element_jacobians(false),

      use_petsc_snes(false),
      time_solver_quiet(true), solver_quiet(true), solver_verbose(false),
      require_residual_reduction(true),
      min_step_length(1e-5),
      max_linear_iterations(200000), max_nonlinear_iterations(20),
      relative_step_tolerance(1.e-7), relative_residual_tolerance(1.e-10),
      initial_linear_tolerance(1.e-3), minimum_linear_tolerance(TOLERANCE*TOLERANCE),
      linear_tolerance_multiplier(1.e-3),

      initial_sobolev_order(1),
      initial_extra_quadrature(0),
      indicator_type("kelly"), sobolev_order(1)
{
}

FEMParameters::~FEMParameters()
{
  for (std::map<subdomain_id_type, FunctionBase<Number> *>::iterator
       i = initial_conditions.begin(); i != initial_conditions.end();
       ++i)
    delete i->second;

  for (std::map<boundary_id_type, FunctionBase<Number> *>::iterator
       i = dirichlet_conditions.begin(); i != dirichlet_conditions.end();
       ++i)
    delete i->second;

  for (std::map<boundary_id_type, FunctionBase<Number> *>::iterator
       i = neumann_conditions.begin(); i != neumann_conditions.end();
       ++i)
    delete i->second;

  for (std::map<int,
	 std::map<boundary_id_type, FunctionBase<Number> *>
       >::iterator
       i = other_boundary_functions.begin(); i != other_boundary_functions.end();
       ++i)
    for (std::map<boundary_id_type, FunctionBase<Number> *>::iterator
         j = i->second.begin(); j != i->second.end();
         ++j)
      delete j->second;

  for (std::map<int,
	 std::map<subdomain_id_type, FunctionBase<Number> *>
       >::iterator
       i = other_interior_functions.begin(); i != other_interior_functions.end();
       ++i)
    for (std::map<subdomain_id_type, FunctionBase<Number> *>::iterator
         j = i->second.begin(); j != i->second.end();
         ++j)
      delete j->second;
}


AutoPtr<FunctionBase<Number> > new_function_base(const std::string& func_type,
                                        const std::string& func_value)
{
  if (func_type == "parsed")
    return AutoPtr<FunctionBase<Number> >
      (new ParsedFunction<Number>(func_value));
  else if (func_type == "factory")
    return AutoPtr<FunctionBase<Number> >
      (Factory<FunctionBase<Number> >::build(func_value).release());
  else if (func_type == "zero")
    return AutoPtr<FunctionBase<Number> >
      (new ZeroFunction<Number>);
  else
    libmesh_not_implemented();

  return AutoPtr<FunctionBase<Number> >();
}


void FEMParameters::read(GetPot &input)
{
    std::vector<std::string> variable_names;

    GETPOT_INT_INPUT(initial_timestep);
    GETPOT_INT_INPUT(n_timesteps);
    GETPOT_INPUT(transient);
    GETPOT_INT_INPUT(deltat_reductions);
    GETPOT_INPUT(timesolver_core);
    GETPOT_INPUT(end_time);
    GETPOT_INPUT(deltat);
    GETPOT_INPUT(timesolver_theta);
    GETPOT_INPUT(timesolver_maxgrowth);
    GETPOT_INPUT(timesolver_tolerance);
    GETPOT_INPUT(timesolver_upper_tolerance);
    GETPOT_INPUT(steadystate_tolerance);

    GETPOT_REGISTER(timesolver_norm);
    const unsigned int n_timesolver_norm = input.vector_variable_size("timesolver_norm");
    timesolver_norm.resize(n_timesolver_norm, L2);
    for (unsigned int i=0; i != n_timesolver_norm; ++i)
      {
        int current_norm = 0; // L2
        if (timesolver_norm[i] == H1)
          current_norm = 1;
        if (timesolver_norm[i] == H2)
          current_norm = 2;
        current_norm = input("timesolver_norm", current_norm, i);
        if (current_norm == 0)
          timesolver_norm[i] = L2;
        else if (current_norm == 1)
          timesolver_norm[i] = H1;
        else if (current_norm == 2)
          timesolver_norm[i] = H2;
        else
          timesolver_norm[i] = DISCRETE_L2;
      }


    GETPOT_INT_INPUT(dimension);
    GETPOT_INPUT(domaintype);
    GETPOT_INPUT(domainfile);
    GETPOT_INPUT(elementtype);
    GETPOT_INPUT(elementorder);
    GETPOT_INPUT(domain_xmin);
    GETPOT_INPUT(domain_ymin);
    GETPOT_INPUT(domain_zmin);
    GETPOT_INPUT(domain_edge_width);
    GETPOT_INPUT(domain_edge_length);
    GETPOT_INPUT(domain_edge_height);
    GETPOT_INT_INPUT(coarsegridx);
    GETPOT_INT_INPUT(coarsegridy);
    GETPOT_INT_INPUT(coarsegridz);
    GETPOT_INT_INPUT(coarserefinements);
    GETPOT_INT_INPUT(extrarefinements);


    GETPOT_INT_INPUT(nelem_target);
    GETPOT_INPUT(global_tolerance);
    GETPOT_INPUT(refine_fraction);
    GETPOT_INPUT(coarsen_fraction);
    GETPOT_INPUT(coarsen_threshold);
    GETPOT_INT_INPUT(max_adaptivesteps);
    GETPOT_INT_INPUT(initial_adaptivesteps);


    GETPOT_INT_INPUT(write_interval);
    GETPOT_INPUT(output_xda);
    GETPOT_INPUT(output_xdr);
    GETPOT_INPUT(output_gz);
#ifndef LIBMESH_HAVE_GZSTREAM
    output_gz                   = false;
#endif
    GETPOT_INPUT(output_bz2);
#ifndef LIBMESH_HAVE_BZ2
    output_bz2                  = false;
#endif
    GETPOT_INPUT(output_gmv);
    GETPOT_INPUT(write_gmv_error);
#ifndef LIBMESH_HAVE_GMV
    output_gmv                  = false;
    write_gmv_error             = false;
#endif
    GETPOT_INPUT(output_tecplot);
    GETPOT_INPUT(write_tecplot_error);
#ifndef LIBMESH_HAVE_TECPLOT_API
    output_tecplot              = false;
    write_tecplot_error         = false;
#endif


    GETPOT_REGISTER(system_types);
    const unsigned int n_system_types =
      input.vector_variable_size("system_types");
    if (n_system_types)
      {
        system_types.resize(n_system_types, "");
        for (unsigned int i=0; i != n_system_types; ++i)
          {
            system_types[i] = input("system_types", system_types[i], i);
          }
      }


    GETPOT_REGISTER(periodic_boundaries);
    const unsigned int n_periodic_bcs =
      input.vector_variable_size("periodic_boundaries");

    if (n_periodic_bcs)
      {
        if (domaintype != "square" && 
            domaintype != "cylinder" && 
            domaintype != "file" && 
            domaintype != "od2")
          {
            libMesh::out << "Periodic boundaries need rectilinear domains" << std::endl;;
            libmesh_error();
          }
        for (unsigned int i=0; i != n_periodic_bcs; ++i)
          {
            unsigned int myboundary =
              input("periodic_boundaries", -1, i);
            unsigned int pairedboundary = 0;
            RealVectorValue translation_vector;
            if (dimension == 2)
              switch (myboundary)
              {
              case 0:
                pairedboundary = 2;
                translation_vector = RealVectorValue(0., domain_edge_length);
                break;
              case 1:
                pairedboundary = 3;
                translation_vector = RealVectorValue(-domain_edge_width, 0);
                break;
              default:
                libMesh::out << "Unrecognized periodic boundary id " <<
                                myboundary << std::endl;;
                libmesh_error();
              }
            else if (dimension == 3)
              switch (myboundary)
              {
              case 0:
                pairedboundary = 5;
                translation_vector = RealVectorValue(0., 0., domain_edge_height);
                break;
              case 1:
                pairedboundary = 3;
                translation_vector = RealVectorValue(0., domain_edge_length, 0.);
                break;
              case 2:
                pairedboundary = 4;
                translation_vector = RealVectorValue(-domain_edge_width, 0., 0.);
                break;
              default:
                libMesh::out << "Unrecognized periodic boundary id " <<
                                myboundary << std::endl;;
                libmesh_error();
              }
            periodic_boundaries.push_back(PeriodicBoundary(translation_vector));
            periodic_boundaries[i].myboundary = myboundary;
            periodic_boundaries[i].pairedboundary = pairedboundary;
          }
      }

    // Use std::string inputs so GetPot doesn't have to make a bunch
    // of internal C string copies
    std::string zero_string = "zero";
    std::string empty_string = "";

    GETPOT_REGISTER(dirichlet_condition_types);
    GETPOT_REGISTER(dirichlet_condition_values);
    GETPOT_REGISTER(dirichlet_condition_boundaries);
    GETPOT_REGISTER(dirichlet_condition_variables);

    const unsigned int n_dirichlet_conditions=
      input.vector_variable_size("dirichlet_condition_types");

    if (n_dirichlet_conditions !=
        input.vector_variable_size("dirichlet_condition_values"))
      {
        libMesh::out << "Error: " << n_dirichlet_conditions
                     << " Dirichlet condition types does not match "
                     << input.vector_variable_size("dirichlet_condition_values")
                     << " Dirichlet condition values." << std::endl;
        
        libmesh_error();
      }

    if (n_dirichlet_conditions !=
        input.vector_variable_size("dirichlet_condition_boundaries"))
      {
        libMesh::out << "Error: " << n_dirichlet_conditions
                     << " Dirichlet condition types does not match "
                     << input.vector_variable_size("dirichlet_condition_boundaries")
                     << " Dirichlet condition boundaries." << std::endl;
        
        libmesh_error();
      }

    if (n_dirichlet_conditions !=
        input.vector_variable_size("dirichlet_condition_variables"))
      {
        libMesh::out << "Error: " << n_dirichlet_conditions
                     << " Dirichlet condition types does not match "
                     << input.vector_variable_size("dirichlet_condition_variables")
                     << " Dirichlet condition variables sets." << std::endl;
        
        libmesh_error();
      }

    for (unsigned int i=0; i != n_dirichlet_conditions; ++i)
      {
        const std::string func_type =
          input("dirichlet_condition_types", zero_string, i);

        const std::string func_value =
          input("dirichlet_condition_values", empty_string, i);

        const boundary_id_type func_boundary =
          input("dirichlet_condition_boundaries", boundary_id_type(0), i);

        dirichlet_conditions[func_boundary] =
          (new_function_base(func_type, func_value).release());

        const std::string variable_set =
          input("dirichlet_condition_variables", empty_string, i);

        for (unsigned int i=0; i != variable_set.size(); ++i)
          {
            if (variable_set[i] == '1')
              dirichlet_condition_variables[func_boundary].push_back(i);
            else if (variable_set[i] != '0')
              {
                libMesh::out << "Unable to understand Dirichlet variable set" 
                             << variable_set << std::endl;
                libmesh_error();
              }
          }
      }

    GETPOT_REGISTER(neumann_condition_types);
    GETPOT_REGISTER(neumann_condition_values);
    GETPOT_REGISTER(neumann_condition_boundaries);
    GETPOT_REGISTER(neumann_condition_variables);

    const unsigned int n_neumann_conditions=
      input.vector_variable_size("neumann_condition_types");

    if (n_neumann_conditions !=
        input.vector_variable_size("neumann_condition_values"))
      {
        libMesh::out << "Error: " << n_neumann_conditions
                     << " Neumann condition types does not match "
                     << input.vector_variable_size("neumann_condition_values")
                     << " Neumann condition values." << std::endl;
        
        libmesh_error();
      }

    if (n_neumann_conditions !=
        input.vector_variable_size("neumann_condition_boundaries"))
      {
        libMesh::out << "Error: " << n_neumann_conditions
                     << " Neumann condition types does not match "
                     << input.vector_variable_size("neumann_condition_boundaries")
                     << " Neumann condition boundaries." << std::endl;
        
        libmesh_error();
      }

    if (n_neumann_conditions !=
        input.vector_variable_size("neumann_condition_variables"))
      {
        libMesh::out << "Error: " << n_neumann_conditions
                     << " Neumann condition types does not match "
                     << input.vector_variable_size("neumann_condition_variables")
                     << " Neumann condition variables sets." << std::endl;
        
        libmesh_error();
      }

    for (unsigned int i=0; i != n_neumann_conditions; ++i)
      {
        const std::string func_type =
          input("neumann_condition_types", zero_string, i);

        const std::string func_value =
          input("neumann_condition_values", empty_string, i);

        const boundary_id_type func_boundary =
          input("neumann_condition_boundaries", boundary_id_type(0), i);

        neumann_conditions[func_boundary] =
          (new_function_base(func_type, func_value).release());

        const std::string variable_set =
          input("neumann_condition_variables", empty_string, i);

        for (unsigned int i=0; i != variable_set.size(); ++i)
          {
            if (variable_set[i] == '1')
              neumann_condition_variables[func_boundary].push_back(i);
            else if (variable_set[i] != '0')
              {
                libMesh::out << "Unable to understand Neumann variable set" 
                             << variable_set << std::endl;
                libmesh_error();
              }
          }
      }

    GETPOT_REGISTER(initial_condition_types);
    GETPOT_REGISTER(initial_condition_values);
    GETPOT_REGISTER(initial_condition_subdomains);

    const unsigned int n_initial_conditions=
      input.vector_variable_size("initial_condition_types");

    if (n_initial_conditions !=
        input.vector_variable_size("initial_condition_values"))
      {
        libMesh::out << "Error: " << n_initial_conditions
                     << " initial condition types does not match "
                     << input.vector_variable_size("initial_condition_values")
                     << " initial condition values." << std::endl;
        
        libmesh_error();
      }

    if (n_initial_conditions !=
        input.vector_variable_size("initial_condition_subdomains"))
      {
        libMesh::out << "Error: " << n_initial_conditions
                     << " initial condition types does not match "
                     << input.vector_variable_size("initial_condition_subdomains")
                     << " initial condition subdomains." << std::endl;
        
        libmesh_error();
      }

    for (unsigned int i=0; i != n_initial_conditions; ++i)
      {
        const std::string func_type =
          input("initial_condition_types", zero_string, i);

        const std::string func_value =
          input("initial_condition_values", empty_string, i);

        const subdomain_id_type func_subdomain =
          input("initial_condition_subdomains", subdomain_id_type(0), i);

        initial_conditions[func_subdomain] =
          (new_function_base(func_type, func_value).release());
      }

    GETPOT_REGISTER(other_interior_function_types);
    GETPOT_REGISTER(other_interior_function_values);
    GETPOT_REGISTER(other_interior_function_subdomains);
    GETPOT_REGISTER(other_interior_function_ids);

    const unsigned int n_other_interior_functions =
      input.vector_variable_size("other_interior_function_types");

    if (n_other_interior_functions !=
        input.vector_variable_size("other_interior_function_values"))
      {
        libMesh::out << "Error: " << n_other_interior_functions
                     << " other interior function types does not match "
                     << input.vector_variable_size("other_interior_function_values")
                     << " other interior function values." << std::endl;
        
        libmesh_error();
      }

    if (n_other_interior_functions !=
        input.vector_variable_size("other_interior_function_subdomains"))
      {
        libMesh::out << "Error: " << n_other_interior_functions
                     << " other interior function types does not match "
                     << input.vector_variable_size("other_interior_function_subdomains")
                     << " other interior function subdomains." << std::endl;
        
        libmesh_error();
      }

    if (n_other_interior_functions !=
        input.vector_variable_size("other_interior_function_ids"))
      {
        libMesh::out << "Error: " << n_other_interior_functions
                     << " other interior function types does not match "
                     << input.vector_variable_size("other_interior_function_ids")
                     << " other interior function ids." << std::endl;
        
        libmesh_error();
      }

    for (unsigned int i=0; i != n_other_interior_functions; ++i)
      {
        const std::string func_type =
          input("other_interior_function_types", zero_string, i);

        const std::string func_value =
          input("other_interior_function_values", empty_string, i);

        const subdomain_id_type func_subdomain =
          input("other_interior_condition_subdomains", subdomain_id_type(0), i);

        const int func_id =
          input("other_interior_condition_ids", int(0), i);

        other_interior_functions[func_id][func_subdomain] =
          (new_function_base(func_type, func_value).release());
      }

    GETPOT_REGISTER(other_boundary_function_types);
    GETPOT_REGISTER(other_boundary_function_values);
    GETPOT_REGISTER(other_boundary_function_boundaries);
    GETPOT_REGISTER(other_boundary_function_ids);

    const unsigned int n_other_boundary_functions =
      input.vector_variable_size("other_boundary_function_types");

    if (n_other_boundary_functions !=
        input.vector_variable_size("other_boundary_function_values"))
      {
        libMesh::out << "Error: " << n_other_boundary_functions
                     << " other boundary function types does not match "
                     << input.vector_variable_size("other_boundary_function_values")
                     << " other boundary function values." << std::endl;
        
        libmesh_error();
      }

    if (n_other_boundary_functions !=
        input.vector_variable_size("other_boundary_function_boundaries"))
      {
        libMesh::out << "Error: " << n_other_boundary_functions
                     << " other boundary function types does not match "
                     << input.vector_variable_size("other_boundary_function_boundaries")
                     << " other boundary function boundaries." << std::endl;
        
        libmesh_error();
      }

    if (n_other_boundary_functions !=
        input.vector_variable_size("other_boundary_function_ids"))
      {
        libMesh::out << "Error: " << n_other_boundary_functions
                     << " other boundary function types does not match "
                     << input.vector_variable_size("other_boundary_function_ids")
                     << " other boundary function ids." << std::endl;
        
        libmesh_error();
      }

    for (unsigned int i=0; i != n_other_boundary_functions; ++i)
      {
        const std::string func_type =
          input("other_boundary_function_types", zero_string, i);

        const std::string func_value =
          input("other_boundary_function_values", empty_string, i);

        const boundary_id_type func_boundary =
          input("other_boundary_function_boundaries", boundary_id_type(0), i);

        const int func_id =
          input("other_boundary_function_ids", int(0), i);

        other_boundary_functions[func_id][func_boundary] =
          (new_function_base(func_type, func_value).release());
      }

    GETPOT_INPUT(run_simulation);
    GETPOT_INPUT(run_postprocess);


    GETPOT_REGISTER(fe_family);
    const unsigned int n_fe_family =
      std::max(1u, input.vector_variable_size("fe_family"));
    fe_family.resize(n_fe_family, "LAGRANGE");
    for (unsigned int i=0; i != n_fe_family; ++i)
      fe_family[i]              = input("fe_family", fe_family[i].c_str(), i);
    GETPOT_REGISTER(fe_order);
    const unsigned int n_fe_order =
      input.vector_variable_size("fe_order");
    fe_order.resize(n_fe_order, 1);
    for (unsigned int i=0; i != n_fe_order; ++i)
      fe_order[i]               = input("fe_order", (int)fe_order[i], i);
    GETPOT_INPUT(extra_quadrature_order);


    GETPOT_INPUT(analytic_jacobians);
    GETPOT_INPUT(verify_analytic_jacobians);
    GETPOT_INPUT(numerical_jacobian_h);
    GETPOT_INPUT(print_solution_norms);
    GETPOT_INPUT(print_solutions);
    GETPOT_INPUT(print_residual_norms);
    GETPOT_INPUT(print_residuals);
    GETPOT_INPUT(print_jacobian_norms);
    GETPOT_INPUT(print_jacobians);
    GETPOT_INPUT(print_element_jacobians);


    GETPOT_INPUT(use_petsc_snes);
    GETPOT_INPUT(time_solver_quiet);
    GETPOT_INPUT(solver_quiet);
    GETPOT_INPUT(solver_verbose);
    GETPOT_INPUT(require_residual_reduction);
    GETPOT_INPUT(min_step_length);
    GETPOT_INT_INPUT(max_linear_iterations);
    GETPOT_INT_INPUT(max_nonlinear_iterations);
    GETPOT_INPUT(relative_step_tolerance);
    GETPOT_INPUT(relative_residual_tolerance);
    GETPOT_INPUT(initial_linear_tolerance);
    GETPOT_INPUT(minimum_linear_tolerance);
    GETPOT_INPUT(linear_tolerance_multiplier);


    GETPOT_INT_INPUT(initial_sobolev_order);
    GETPOT_INT_INPUT(initial_extra_quadrature);
    GETPOT_INPUT(indicator_type);
    GETPOT_INT_INPUT(sobolev_order);

    GETPOT_INPUT(system_config_file);

  std::vector<std::string> bad_variables =
    input.unidentified_arguments(variable_names);

  if (libMesh::processor_id() == 0 && !bad_variables.empty())
    {
      std::cerr << "ERROR: Unrecognized variables:" << std::endl;
      for (unsigned int i = 0; i != bad_variables.size(); ++i)
        std::cerr << bad_variables[i] << std::endl;
      std::cerr << "Not found among recognized variables:" << std::endl;
      for (unsigned int i = 0; i != variable_names.size(); ++i)
        std::cerr << variable_names[i] << std::endl;
      libmesh_error();
    }
}
