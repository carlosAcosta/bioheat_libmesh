
// Local includes
#include "domain.h"
#include "initial.h"
#include "femparameters.h"
#include "mysystems.h"
#include "output.h"

// Libmesh includes
#include "dirichlet_boundaries.h"
#include "equation_systems.h"
#include "error_vector.h"
#include "euler_solver.h"
#include "euler2_solver.h"
#include "exact_error_estimator.h"
#include "fourth_error_estimators.h"
#include "getpot.h"
#include "gmv_io.h"
#include "kelly_error_estimator.h"
#include "mesh.h"
#include "mesh_refinement.h"
#include "newton_solver.h"
#include "numeric_vector.h"
#include "petsc_diff_solver.h"
#include "steady_solver.h"
#include "system_norm.h"
#include "tecplot_io.h"
#include "twostep_time_solver.h"
#include "uniform_refinement_estimator.h"

// Some (older) compilers do not offer full stream
// functionality, OStringStream works around this.
#include "o_string_stream.h"

// C++ includes
#include <iostream>
#include <sys/time.h>

std::string numbered_filename(unsigned int number,
                              std::string type,
                              std::string extension,
                              FEMParameters &param)
{
  OStringStream file_name;
  file_name << "out." << type << '.';
  OSSRealzeroright(file_name, 4, 0, number);
  file_name << '.' << extension;
  if (param.output_bz2)
    file_name << ".bz2";
  else if (param.output_gz)
    file_name << ".gz";
  return file_name.str();
}


void write_output(EquationSystems &es,
                  unsigned int number,
                  FEMParameters &param)
{
  MeshBase &mesh = es.get_mesh();

#ifdef LIBMESH_HAVE_GMV
  if (param.output_gmv)
    {
      OStringStream file_name_gmv;
      file_name_gmv << "out.gmv.";
      OSSRealzeroright(file_name_gmv,4,0,number);

      GMVIO(mesh).write_equation_systems
        (file_name_gmv.str(), es);
    }
#endif

#ifdef LIBMESH_HAVE_TECPLOT_API
  if (param.output_tecplot)
    {
      OStringStream file_name_tecplot;
      file_name_tecplot << "out.";
      OSSRealzeroright(file_name_tecplot,4,0,number);
      file_name_tecplot << ".plt";

      TecplotIO(mesh).write_equation_systems
        (file_name_tecplot.str(), es);
    }
#endif

  if (param.output_xda || param.output_xdr)
    mesh.renumber_nodes_and_elements();
  if (param.output_xda)
    {
      mesh.write(numbered_filename(number, "mesh", "xda", param));
      es.write(numbered_filename(number, "soln", "xda", param),
               EquationSystems::WRITE_DATA |
               EquationSystems::WRITE_ADDITIONAL_DATA);
    }
  if (param.output_xdr)
    {
      mesh.write(numbered_filename(number, "mesh", "xdr", param));
      es.write(numbered_filename(number, "soln", "xdr", param),
               EquationSystems::WRITE_DATA |
               EquationSystems::WRITE_ADDITIONAL_DATA);
    }
}

void write_output_headers(FEMParameters &param)
{
  // Only one processor needs to take care of headers.
  if (libMesh::processor_id() != 0)
    return;

  start_output(param.initial_timestep, "out_clocktime.m", "vector_clocktime");

  if (param.run_simulation)
    {
      start_output(param.initial_timestep, "out_activemesh.m", "table_activemesh");
      start_output(param.initial_timestep, "out_solvesteps.m", "table_solvesteps");

      if (param.timesolver_tolerance)
        {
          start_output(param.initial_timestep, "out_time.m", "vector_time");
          start_output(param.initial_timestep, "out_timesteps.m", "vector_timesteps");
        }
      if (param.steadystate_tolerance)
        start_output(param.initial_timestep, "out_changerate.m", "vector_changerate");
    }
}


void write_output_solvedata(EquationSystems &es,
                            unsigned int a_step,
                            unsigned int newton_steps,
                            unsigned int krylov_steps,
                            unsigned int tv_sec,
                            unsigned int tv_usec)
{
  MeshBase &mesh = es.get_mesh();
  unsigned int n_active_elem = mesh.n_active_elem();
  unsigned int n_active_dofs = es.n_active_dofs();

  if (libMesh::processor_id() == 0)
    {
      // Write out the number of elements/dofs used
      std::ofstream activemesh ("out_activemesh.m",
        std::ios_base::app | std::ios_base::out);
      activemesh.precision(17);
      activemesh << (a_step + 1) << ' '
                 << n_active_elem << ' '
                 << n_active_dofs << std::endl;

      // Write out the number of solver steps used
      std::ofstream solvesteps ("out_solvesteps.m",
        std::ios_base::app | std::ios_base::out);
      solvesteps.precision(17);
      solvesteps << newton_steps << ' '
                 << krylov_steps << std::endl;

      // Write out the clock time
      std::ofstream clocktime ("out_clocktime.m",
        std::ios_base::app | std::ios_base::out);
      clocktime.precision(17);
      clocktime << tv_sec << '.' << tv_usec << std::endl;
    }
}


void write_output_footers(FEMParameters &param)
{
  // Write footers on output .m files
  if (libMesh::processor_id() == 0)
    {
      std::ofstream clocktime ("out_clocktime.m",
        std::ios_base::app | std::ios_base::out);
      clocktime << "];" << std::endl;

      if (param.run_simulation)
        {
          std::ofstream activemesh ("out_activemesh.m",
            std::ios_base::app | std::ios_base::out);
          activemesh << "];" << std::endl;

          std::ofstream solvesteps ("out_solvesteps.m",
            std::ios_base::app | std::ios_base::out);
          solvesteps << "];" << std::endl;

          if (param.timesolver_tolerance)
            {
              std::ofstream times ("out_time.m",
                std::ios_base::app | std::ios_base::out);
              times << "];" << std::endl;
              std::ofstream timesteps ("out_timesteps.m",
                std::ios_base::app | std::ios_base::out);
              timesteps << "];" << std::endl;
            }
          if (param.steadystate_tolerance)
            {
              std::ofstream changerate ("out_changerate.m",
                std::ios_base::app | std::ios_base::out);
              changerate << "];" << std::endl;
            }
        }

      // We're done, so write out a file (for e.g. Make to check)
      std::ofstream complete ("complete");
      complete << "complete" << std::endl;
    }
}

#if defined(LIBMESH_HAVE_GMV) || defined(LIBMESH_HAVE_TECPLOT_API)
void write_error(EquationSystems &es,
                 ErrorVector &error,
                 unsigned int t_number,
                 unsigned int a_number,
                 FEMParameters &param)
#else
void write_error(EquationSystems &,
                 ErrorVector &,
                 unsigned int,
                 unsigned int,
                 FEMParameters &)
#endif
{
#ifdef LIBMESH_HAVE_GMV
  if (param.write_gmv_error)
    {
      OStringStream error_gmv;
      error_gmv << "error.gmv.";
      OSSRealzeroright(error_gmv,4,0, a_number);
      error_gmv << ".";
      OSSRealzeroright(error_gmv,4,0, t_number);
      error.plot_error(error_gmv.str(), es.get_mesh());
    }
#endif

#ifdef LIBMESH_HAVE_TECPLOT_API
  if (param.write_tecplot_error)
    {
      OStringStream error_tecplot;
      error_tecplot << "error.plt.";
      OSSRealzeroright(error_tecplot,4,0, a_number);
      error_tecplot << ".";
      OSSRealzeroright(error_tecplot,4,0, t_number);
      error.plot_error(error_tecplot.str(), es.get_mesh());
    }
#endif
}


void read_output(EquationSystems &es,
                 unsigned int number,
                 FEMParameters &param)
{
  MeshBase &mesh = es.get_mesh();

  std::string file_name_mesh, file_name_soln;
  // Look for ASCII files first
  if (param.output_xda)
    {
      file_name_mesh = numbered_filename(number, "mesh", "xda", param);
      file_name_soln = numbered_filename(number, "soln", "xda", param);
    }
  else if (param.output_xdr)
    {
      file_name_mesh = numbered_filename(number, "mesh", "xdr", param);
      file_name_soln = numbered_filename(number, "soln", "xdr", param);
    }
  else
    libmesh_error();

  // Read in the mesh
  mesh.read(file_name_mesh);

  // And the stored solution
  es.read(file_name_soln,
          EquationSystems::READ_HEADER |
          EquationSystems::READ_DATA |
          EquationSystems::READ_ADDITIONAL_DATA);

  // Put systems in a consistent state
  for (unsigned int i = 0; i != es.n_systems(); ++i)
    es.get_system<FEMSystem>(i).update();

  // Figure out the current time
  Real current_time = 0., current_timestep = 0.;

  if (param.timesolver_tolerance)
    {
      std::ifstream times ("out_time.m");
      std::ifstream timesteps ("out_timesteps.m");
      if (times.is_open() && timesteps.is_open())
        {
          // Read headers
          const unsigned int headersize = 25;
          char header[headersize];
          timesteps.getline (header, headersize);
          if (strcmp(header, "vector_timesteps = [") != 0)
            {
              std::cout << "Bad header in out_timesteps.m:" << std::endl
                        << header
                        << std::endl;
              libmesh_error();
            }

          times.getline (header, headersize);
          if (strcmp(header, "vector_time = [") != 0)
            {
              std::cout << "Bad header in out_time.m:" << std::endl
                        << header
                        << std::endl;
              libmesh_error();
            }

          // Read each timestep
          for (unsigned int i = 0; i != number; ++i)
            {
              if (!times.good())
                libmesh_error();
              times >> current_time;
              timesteps >> current_timestep;
            }
          // Remember to increment the last timestep; out_times.m
          // lists each *start* time
          current_time += current_timestep;
        }
      else
        libmesh_error();
    }
  else
    current_time = number * param.deltat;

  for (unsigned int i = 0; i != es.n_systems(); ++i)
    es.get_system<FEMSystem>(i).time = current_time;
}

void set_system_parameters(FEMSystem &system, FEMParameters &param)
{
  // Verify analytic jacobians against numerical ones?
  system.verify_analytic_jacobians = param.verify_analytic_jacobians;
  system.numerical_jacobian_h = param.numerical_jacobian_h;

  // More desperate debugging options
  system.print_solution_norms    = param.print_solution_norms;
  system.print_solutions         = param.print_solutions;
  system.print_residual_norms    = param.print_residual_norms;
  system.print_residuals         = param.print_residuals;
  system.print_jacobian_norms    = param.print_jacobian_norms;
  system.print_jacobians         = param.print_jacobians;
  system.print_element_jacobians = param.print_element_jacobians;

  // Solve this as a time-dependent or steady system
  if (param.transient)
    {
      UnsteadySolver *innersolver;
      if (param.timesolver_core == "euler2")
        {
          Euler2Solver *euler2solver =
            new Euler2Solver(system);

          euler2solver->theta = param.timesolver_theta;
          innersolver = euler2solver;
        }
      else if (param.timesolver_core == "euler")
        {
          EulerSolver *eulersolver =
            new EulerSolver(system);

          eulersolver->theta = param.timesolver_theta;
          innersolver = eulersolver;
        }
      else
        {
          std::cerr << "Don't recognize core TimeSolver type: "
                    << param.timesolver_core << std::endl;
          libmesh_error();
        }

      if (param.timesolver_tolerance)
        {
          TwostepTimeSolver *timesolver =
            new TwostepTimeSolver(system);

          timesolver->max_growth       = param.timesolver_maxgrowth;
          timesolver->target_tolerance = param.timesolver_tolerance;
          timesolver->upper_tolerance  = param.timesolver_upper_tolerance;
          timesolver->component_norm   = SystemNorm(param.timesolver_norm);

          timesolver->core_time_solver =
            AutoPtr<UnsteadySolver>(innersolver);
          system.time_solver =
            AutoPtr<UnsteadySolver>(timesolver);
        }
      else
        system.time_solver =
          AutoPtr<TimeSolver>(innersolver);
    }
  else
    system.time_solver =
      AutoPtr<TimeSolver>(new SteadySolver(system));

  system.time_solver->reduce_deltat_on_diffsolver_failure =
                                        param.deltat_reductions;
  system.time_solver->quiet           = param.time_solver_quiet;

  // Create any periodic boundary conditions
  for (unsigned int i=0; i != param.periodic_boundaries.size(); ++i)
    system.get_dof_map().add_periodic_boundary(param.periodic_boundaries[i]);

  // Create any Dirichlet boundary conditions
  typedef 
    std::map<boundary_id_type, FunctionBase<Number> *>::
    const_iterator Iter;

  for (Iter i = param.dirichlet_conditions.begin(); 
       i != param.dirichlet_conditions.end(); ++i)
    {
      boundary_id_type b = i->first;
      FunctionBase<Number> *f = i->second;
      std::set<boundary_id_type> bdys; bdys.insert(b);
      system.get_dof_map().add_dirichlet_boundary
        (DirichletBoundary
           (bdys, param.dirichlet_condition_variables[b], f));
      std::cout << "Added Dirichlet boundary " << b << " for variables ";
      for (unsigned int vi=0; vi !=
           param.dirichlet_condition_variables[b].size(); ++vi)
        std::cout << param.dirichlet_condition_variables[b][vi];
      std::cout << std::endl;
    }

  // Set the time stepping options
  system.deltat = param.deltat;

  // And the integration options
  system.extra_quadrature_order = param.extra_quadrature_order;

  // And the nonlinear solver options
  if (param.use_petsc_snes)
    {
#ifdef LIBMESH_HAVE_PETSC
      PetscDiffSolver *solver = new PetscDiffSolver(system);
      system.time_solver->diff_solver() = AutoPtr<DiffSolver>(solver);

      solver->quiet                       = param.solver_quiet;
      solver->verbose                     = param.solver_verbose;
      solver->max_nonlinear_iterations    = param.max_nonlinear_iterations;
      //solver->minsteplength               = param.min_step_length;
      solver->relative_step_tolerance     = param.relative_step_tolerance;
      solver->relative_residual_tolerance = param.relative_residual_tolerance;
      //solver->require_residual_reduction  = param.require_residual_reduction;
      //solver->linear_tolerance_multiplier = param.linear_tolerance_multiplier;
      if (system.time_solver->reduce_deltat_on_diffsolver_failure)
        {
          solver->continue_after_max_iterations = true;
          solver->continue_after_backtrack_failure = true;
        }

      // And the linear solver options
      solver->max_linear_iterations       = param.max_linear_iterations;
      solver->initial_linear_tolerance    = param.initial_linear_tolerance;
      solver->minimum_linear_tolerance    = param.minimum_linear_tolerance;

#else
      libmesh_error();
#endif
    }
  else
    {
      NewtonSolver *solver = new NewtonSolver(system);
      system.time_solver->diff_solver() = AutoPtr<DiffSolver>(solver);

      solver->quiet                       = param.solver_quiet;
      solver->verbose                     = param.solver_verbose;
      solver->max_nonlinear_iterations    = param.max_nonlinear_iterations;
      solver->minsteplength               = param.min_step_length;
      solver->relative_step_tolerance     = param.relative_step_tolerance;
      solver->relative_residual_tolerance = param.relative_residual_tolerance;
      solver->require_residual_reduction  = param.require_residual_reduction;
      solver->linear_tolerance_multiplier = param.linear_tolerance_multiplier;
      if (system.time_solver->reduce_deltat_on_diffsolver_failure)
        {
          solver->continue_after_max_iterations = true;
          solver->continue_after_backtrack_failure = true;
        }

      // And the linear solver options
      solver->max_linear_iterations       = param.max_linear_iterations;
      solver->initial_linear_tolerance    = param.initial_linear_tolerance;
      solver->minimum_linear_tolerance    = param.minimum_linear_tolerance;
    }

/*
  NewtonSolver *adjoint_solver = new NewtonSolver(system);
  system.time_solver->linear_solver() = AutoPtr<DiffSolver>(adjoint_solver);

  // We'll assume some consistent parameter settings for adjoint
  // solves for now
  adjoint_solver->quiet                       = param.solver_quiet;
  adjoint_solver->relative_residual_tolerance = param.relative_residual_tolerance;
  adjoint_solver->max_linear_iterations       = param.max_linear_iterations;
*/
}



AutoPtr<MeshRefinement> build_mesh_refinement(MeshBase &mesh,
                                              FEMParameters &param)
{
  AutoPtr<MeshRefinement> mesh_refinement(new MeshRefinement(mesh));
  mesh_refinement->coarsen_by_parents() = true;
  mesh_refinement->absolute_global_tolerance() = param.global_tolerance;
  mesh_refinement->nelem_target()      = param.nelem_target;
  mesh_refinement->refine_fraction()   = param.refine_fraction;
  mesh_refinement->coarsen_fraction()  = param.coarsen_fraction;
  mesh_refinement->coarsen_threshold() = param.coarsen_threshold;

  return mesh_refinement;
}



AutoPtr<ErrorEstimator> build_error_estimator(FEMParameters &param)
{
  AutoPtr<ErrorEstimator> error_estimator;

  if (param.indicator_type == "uniform")
    {
      UniformRefinementEstimator *u =
        new UniformRefinementEstimator;

      if (param.sobolev_order == 0)
        u->error_norm = L2;
      else if (param.sobolev_order == 1)
        u->error_norm = H1;
      else if (param.sobolev_order == 2)
        u->error_norm = H2;

      error_estimator.reset(u);
    }
  else if (param.indicator_type == "kelly")
    {
      error_estimator.reset(new KellyErrorEstimator);
    }
  else if (param.indicator_type == "laplacianjump")
    {
      error_estimator.reset(new LaplacianErrorEstimator);
    }
  else if (param.indicator_type == "exact")
    {
      ExactErrorEstimator *exact_estimator =
        new ExactErrorEstimator;
      error_estimator.reset(exact_estimator);

      if (param.sobolev_order == 0)
      exact_estimator->error_norm = L2;
      else if (param.sobolev_order == 1)
      exact_estimator->error_norm = H1;
      else if (param.sobolev_order == 2)
      exact_estimator->error_norm = H2;

      exact_estimator->extra_quadrature_order(param.extra_quadrature_order);
      exact_estimator->attach_exact_value(exact_value);
      exact_estimator->attach_exact_deriv(exact_grad);
    }
  else
    {
      std::cerr << "Unknown indicator_type" << std::endl;
      libmesh_error();
    }
  return error_estimator;
}



// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  std::cout << "Started " << argv[0] << std::endl;

  // Make sure the general input file exists, and parse it
  {
    std::ifstream i("general.in");
    if (!i)
      {
        std::cerr << '[' << libMesh::processor_id() 
                  << "] Can't find general.in; exiting early." 
                  << std::endl; 
        libmesh_error();
      }
  }
  GetPot infile("general.in");

  // Read in parameters from the input file
  FEMParameters param;
  param.read(infile);

  Mesh mesh (param.dimension);

  // And an object to refine it
  AutoPtr<MeshRefinement> mesh_refinement =
    build_mesh_refinement(mesh, param);

  // And an EquationSystems to run on it
  EquationSystems equation_systems (mesh);

  std::cout << "Building mesh" << std::endl;

  // If our initial timestep is zero, we need to generate a mesh and
  // system
  if (!param.initial_timestep)
    build_domain(mesh, param);

  std::cout << "Building system" << std::endl;

  FEMSystem &system = build_system(equation_systems, infile, param);

  set_system_parameters(system, param);

  std::cout << "Initializing systems" << std::endl;

  // If our initial timestep is non-zero, or if we're only
  // postprocessing, we want to read a restart file.
  if (param.initial_timestep || !param.run_simulation)
    read_output(equation_systems, param.initial_timestep, param);
  else
    {
      // Just initialize the system
      equation_systems.init ();
    }

  std::cout << "Setting initial conditions" << std::endl;

  // Set the initial conditions for transient problems
  if (param.transient && !param.initial_timestep &&
      param.run_simulation)
    {
      read_initial_parameters();

      // Adaptively find a mesh for the initial conditions
      unsigned int a_step = 0;
      for (; a_step != param.initial_adaptivesteps; ++a_step)
        {
          // We can't adapt to both a tolerance and a
          // target mesh size
          if (param.global_tolerance != 0.)
            libmesh_assert (param.nelem_target == 0);
          // If we aren't adapting to a tolerance we need a
          // target mesh size
          else
            libmesh_assert (param.nelem_target > 0);

          system.project_solution(initial_value, initial_grad,
                                  equation_systems.parameters);

          typedef 
            std::map<boundary_id_type, FunctionBase<Number> *>::
            const_iterator Iter;

          for (Iter i = param.dirichlet_conditions.begin(); 
               i != param.dirichlet_conditions.end(); ++i)
            {
              boundary_id_type b = i->first;
              FunctionBase<Number> *f = i->second;
              std::set<boundary_id_type> bdys; bdys.insert(b);
              system.boundary_project_solution
                (bdys, param.dirichlet_condition_variables[b], f);
            }

          ErrorVector error;
          ExactErrorEstimator initial_estimator;
          if (param.sobolev_order == 0)
            initial_estimator.error_norm = L2;
          else if (param.sobolev_order == 1)
            initial_estimator.error_norm = H1;
          else if (param.sobolev_order == 2)
            initial_estimator.error_norm = H2;

          initial_estimator.extra_quadrature_order(param.initial_extra_quadrature);
          initial_estimator.attach_exact_value(initial_value);
          initial_estimator.attach_exact_deriv(initial_grad);
          initial_estimator.estimate_error(system, error, NULL, true);

          // Write error values if requested
          if (param.write_gmv_error || param.write_tecplot_error)
            write_error(equation_systems, error, 0, a_step+1, param);

          // Print out status at each initial adaptive step.
          Real global_error = error.l2_norm();
          std::cout << "Initial adaptive step " << a_step << ": ";
          if (param.global_tolerance != 0.)
            std::cout << "global_error = " << global_error
                      << " with ";
          std::cout << mesh.n_active_elem()
                    << " active elements and "
                    << equation_systems.n_active_dofs()
                    << " active dofs." << std::endl;
          if (param.global_tolerance != 0.)
            std::cout << "worst element error = " << error.maximum()
                      << ", mean = " << error.mean() << std::endl;

          if (param.global_tolerance != 0.)
            {
              // If we've reached our desired tolerance, we
              // don't need any more adaptive steps
              if (global_error < param.global_tolerance)
                break;
              mesh_refinement->flag_elements_by_error_tolerance (error);
            }
          else
            {
              // If flag_elements_to_nelem_target returns true, this
              // should be our last adaptive step.
              if (mesh_refinement->flag_elements_by_nelem_target (error))
                {
                  mesh_refinement->refine_and_coarsen_elements();
                  equation_systems.reinit();
                  break;
                }
            }

          // Carry out the adaptive mesh refinement/coarsening
          mesh_refinement->refine_and_coarsen_elements();
          equation_systems.reinit();
        }

      if (param.initial_adaptivesteps)
        {
          std::cout << "Initial adaptivity done: ";
          std::cout << mesh.n_active_elem()
                    << " active elements and "
                    << equation_systems.n_active_dofs()
                    << " active dofs." << std::endl;
        }

      // Do one last projection if necessary
      if (a_step == param.initial_adaptivesteps)
        {
          system.project_solution(initial_value, initial_grad,
                                  equation_systems.parameters);
        }

      // Close up any resources initial.C needed
      finish_initialization();

      // Refine the grid again if requested
      for (unsigned int i=0; i != param.extrarefinements; ++i)
        {
          mesh_refinement->uniformly_refine(1);
          equation_systems.reinit();
        }

      // Plot the initial conditions
      write_output(equation_systems, 0, param);

      // Postprocess the initial conditions
      if (param.run_postprocess)
        system.postprocess();
    }

  // Print information about the mesh and system to the screen.
  mesh.print_info();
  equation_systems.print_info();

  write_output_headers(param);

  // In optimized mode we catch any solver errors, so that we can
  // write the proper footers before closing.  In debug mode we just
  // let the exception throw so that gdb can grab it.
#ifdef NDEBUG
  try
  {
#endif
    // Now we begin the timestep loop to compute the time-accurate
    // solution of the equations.
    for (unsigned int t_step=param.initial_timestep;
         t_step != param.initial_timestep + param.n_timesteps; ++t_step)
      {
        if (param.run_simulation)
        {
        // A pretty update message
        std::cout << " Solving time step " << t_step << ", time = "
                  << system.time << std::endl;

        // Time ourselves
        struct timeval tv_start;
        gettimeofday(&tv_start, NULL);

        // Count the number of stups used
        unsigned int newton_steps = 0;
        unsigned int krylov_steps = 0;
        unsigned int a_step = 0;

        // Adaptively solve the timestep
        for (; a_step != param.max_adaptivesteps; ++a_step)
          {
            // We can't adapt to both a tolerance and a
            // target mesh size
            if (param.global_tolerance != 0.)
              libmesh_assert (param.nelem_target == 0);
            // If we aren't adapting to a tolerance we need a
            // target mesh size
            else
              libmesh_assert (param.nelem_target > 0);

            system.solve();
            newton_steps +=
              system.time_solver->diff_solver()->total_outer_iterations();
            krylov_steps +=
              system.time_solver->diff_solver()->total_inner_iterations();

            ErrorVector error;

            AutoPtr<ErrorEstimator> error_estimator =
              build_error_estimator(param);

            error_estimator->estimate_error(system, error);

            // Print out status at each adaptive step.
            Real global_error = error.l2_norm();
            std::cout << "Time step " << t_step 
                      << ", adaptive step " << a_step << ": ";
            if (param.global_tolerance != 0.)
              std::cout << "global_error = " << global_error
                        << " with ";
            std::cout << mesh.n_active_elem()
                      << " active elements and "
                      << equation_systems.n_active_dofs()
                      << " active dofs." << std::endl;
            if (param.global_tolerance != 0.)
              std::cout << "worst element error = " << error.maximum()
                        << ", mean = " << error.mean() << std::endl;

            // Write error values if requested
            if ((t_step+1)%param.write_interval == 0 &&
                (param.write_gmv_error || param.write_tecplot_error))
              write_error(equation_systems, error, t_step+1, a_step+1, param);

            if (param.global_tolerance != 0.)
              {
                // If we've reached our desired tolerance, we
                // don't need any more adaptive steps
                if (global_error < param.global_tolerance)
                  break;
                mesh_refinement->flag_elements_by_error_tolerance (error);
              }
            else
              {
                // If flag_elements_to_nelem_target returns true, this
                // should be our last adaptive step.
                if (mesh_refinement->flag_elements_by_nelem_target (error))
                  {
                    mesh_refinement->refine_and_coarsen_elements();
                    equation_systems.reinit();
                    break;
                  }
              }

            // Carry out the adaptive mesh refinement/coarsening
            mesh_refinement->refine_and_coarsen_elements();
            equation_systems.reinit();
          }

        if (param.max_adaptivesteps > 1)
          {
            std::cout << "Time step " << t_step << " adaptivity done: ";
            std::cout << mesh.n_active_elem()
                      << " active elements and "
                      << equation_systems.n_active_dofs()
                      << " active dofs." << std::endl;
          }

        // Do one last solve if necessary
        if (a_step == param.max_adaptivesteps)
          {
            system.solve();
            newton_steps +=
              system.time_solver->diff_solver()->total_outer_iterations();
            krylov_steps +=
              system.time_solver->diff_solver()->total_inner_iterations();
          }

        struct timeval tv_end;
        gettimeofday(&tv_end, NULL);
        write_output_solvedata(equation_systems,
                               a_step, newton_steps, krylov_steps,
                               tv_end.tv_sec - tv_start.tv_sec,
                               tv_end.tv_usec - tv_start.tv_usec);

        // Write out the time of this timestep if that's adaptive
        if (param.timesolver_tolerance && libMesh::processor_id() == 0)
          {
            std::ofstream times ("out_time.m",
              std::ios_base::app | std::ios_base::out);
            times.precision(17);
            times << system.time << std::endl;
            std::ofstream timesteps ("out_timesteps.m",
              std::ios_base::app | std::ios_base::out);
            timesteps.precision(17);
            timesteps << system.deltat << std::endl;
          }

        // Double check that we aren't approaching steady state
        if (param.steadystate_tolerance)
          {
            Real dunorm = system.time_solver->du(param.timesolver_norm);
            if (dunorm/system.deltat < param.steadystate_tolerance)
              {
                std::cout << " Reached STEADY STATE at time step "
                          << t_step << ", time = " << system.time
                          << ", rate of change = " 
                          << (dunorm/system.deltat) << std::endl;

                // Write any final output
                write_output(equation_systems, t_step+1, param);

                // A pretty update message
                std::cout << " Postprocessing time step " << t_step <<
                             ", time = " << system.time << std::endl;

                if (param.run_postprocess)
                  system.postprocess();

                // Close our output, then exit
                // by breaking out of for(t_step) loop
                break;
              }
            else if (libMesh::processor_id() == 0)
              {
                std::ofstream changerate ("out_changerate.m",
                  std::ios_base::app | std::ios_base::out);
                changerate.precision(17);
                changerate << (dunorm/system.deltat) << std::endl;
              }
          }

        // Advance to the next timestep in a transient problem
        system.time_solver->advance_timestep();

        // Double check that we aren't at our final time
        if (system.time + (param.end_time * TOLERANCE) >
            param.end_time)
          {
            // Write any final output
            write_output(equation_systems, t_step+1, param);

            // A pretty update message
            std::cout << " Postprocessing time step " << t_step << 
                         ", time = " << system.time << std::endl;

            if (param.run_postprocess)
              system.postprocess();

            // Close our output, then exit
            // by breaking out of for(t_step) loop
            break;
          }

        // Write out this timestep if we're requested to
        if ((t_step+1)%param.write_interval == 0)
          write_output(equation_systems, t_step+1, param);

        // Double check that we aren't about to overshoot our final
        // time
        if (system.time + system.deltat > param.end_time)
          {
            system.deltat = param.end_time - system.time;
          }
        }
        else // !run_simulation
        {
        // Figure out the current time
        Real current_time = 0.;

        if (param.timesolver_tolerance)
          {
            std::ifstream times ("out_time.m");
            if (times.is_open())
              {
                for (unsigned int i = 0; i != t_step + 1; ++i)
                  {
                    if (!times.good())
                      libmesh_error();
                    times >> current_time;
                  }
              }
            system.time = current_time;
          }
        else if (system.time + system.deltat > param.end_time)
          system.time = param.end_time;
        else
          system.time += system.deltat;
        }

        // A pretty update message
        std::cout << " Postprocessing time step " << t_step << ", time = "
                  << system.time << std::endl;

        // Post process this timestep if we're requested to
        if (param.run_postprocess)
          system.postprocess();
      }
#ifdef NDEBUG
  }
  catch (...)
  {
    std::cerr << '[' << libMesh::processor_id()
              << "] Caught exception; exiting early." << std::endl; 
  }
#endif

  std::cerr << '[' << libMesh::processor_id() 
            << "] Completing output." << std::endl; 

  write_output_footers(param);
  
  // All done.  
  return 0;
}
