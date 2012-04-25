/* $Id: naviersystem.h,v 1.1 2007-01-15 03:20:07 roystgnr Exp $ */

/* Copyright (C) 2012  Roy H. Stogner */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

// DiffSystem framework files
#include "fem_system.h"

// The Bio-Heat-Transfer system class.
// FEMSystem, TimeSolver and NewtonSolver will handle most tasks,
// but we must specify element residuals
class BioheatSystem : public FEMSystem
{
public:
  // Constructor
  BioheatSystem(EquationSystems& es,
               const std::string& name,
               const unsigned int number)
  : FEMSystem(es, name, number),
    _rho_L(1.), _rho_S(1.),
    _c_L(1.), _c_S(1.),
    _k_L(1.), _k_S(1.),
    _r_L(0.), _r_S(0.),
    _D_c(1.),
    _mu(1.),
    _evaluate_T(true),
    _pressure_penalty(1e9),
    _pressure_point(0),
    _fe_family("LAGRANGE"), _fe_order(1),
    _analytic_jacobians(true),
    _system_config_file("bioheat.in") {}

  // Discretization to use
  std::string & fe_family() { return _fe_family;  }
  unsigned int & fe_order() { return _fe_order;  }

  // Config file from which to read parameters
  std::string& system_config_file() { return _system_config_file; }

  // Constant material parameters to use
  Real & rho_L() { return _rho_L;  }
  Real & rho_S() { return _rho_S;  }
  Real & c_L() { return _c_L;  }
  Real & c_S() { return _c_S;  }
  Real & k_L() { return _k_L;  }
  Real & k_S() { return _k_S;  }
  Real & r_L() { return _r_L;  }
  Real & r_S() { return _r_S;  }
  Real & D_c() { return _D_c;  }
  Real & mu() { return _mu;  }
  bool & analytic_jacobians() { return _analytic_jacobians; }

  // Parameters for pinned pressure
  Real & pressure_penalty() { return _pressure_penalty; }
  Point & pressure_point() { return _pressure_point; }

  // Functions returning liquid fraction values
  std::map<subdomain_id_type, FunctionBase<Number> *> n_L_values;

  // Functions returning boundary condition values and fluxes
  std::map<boundary_id_type, FunctionBase<Number> *> dirichlet_values,
                                                     neumann_values;

  // Vectors specifying which variables to subject to each boundary
  // condition
  std::map<boundary_id_type, std::vector<unsigned int> > dirichlet_variables,
                                                         neumann_variables;

protected:
  // System initialization
  virtual void init_data ();

  // Context initialization
  virtual void init_context (DiffContext &context);

  // Element residual and jacobian calculations
  // Time dependent parts
  virtual bool element_time_derivative (bool request_jacobian,
                                        DiffContext& context);
  virtual bool side_time_derivative (bool request_jacobian,
                                     DiffContext& context);

  // Constraint parts
  virtual bool element_constraint (bool request_jacobian,
                                   DiffContext& context);
  virtual bool side_constraint (bool request_jacobian,
                                DiffContext& context);

  // Indices for each variable;
  unsigned int p_var, u_var, v_var, w_var, T_var;

  // Config file name from which to read parameters
  std::string _system_config_file;

  // The densities of liquid and solid
  Real _rho_L, _rho_S;

  // The specific heat capacities of liquid and solid
  Real _c_L, _c_S;

  // The heat conductivities of liquid and solid
  Real _k_L, _k_S;

  // The external heating supplied to liquid and solid
  Real _r_L, _r_S;

  // The characteristic particle diameter of the solid
  Real _D_c;

  // The dynamic viscosity of the liquid
  Real _mu;

  // Whether to evaluate temperature
  bool _evaluate_T;

  // How strongly (and whether) to pin the pressure at a point
  Real _pressure_penalty;

  // What point to pin the pressure at
  Point _pressure_point;

  // The FE type to use
  std::string _fe_family;
  unsigned int _fe_order;

  // Calculate Jacobians analytically or not?
  bool _analytic_jacobians;
};
