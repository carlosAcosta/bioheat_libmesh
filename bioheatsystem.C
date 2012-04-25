/* $Id: naviersystem.C,v 1.3 2008-04-13 03:35:57 roystgnr Exp $ */

/* Copyright (C) 2012  Roy H. Stogner */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

#include "getpot.h"

#include "bioheatsystem.h"

#include "boundary_info.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "fem_context.h"
#include "function_base.h"
#include "mesh.h"
#include "quadrature.h"
#include "string_to_enum.h"

#include <algorithm>


void BioheatSystem::init_data ()
{
  // Make sure the system config file exists, and parse it.
  {
    std::ifstream i(_system_config_file.c_str());
    if (!i)
      {
        libMesh::out << "Can't find physics config file \"" 
		     << _system_config_file << "\"; exiting early." 
                     << std::endl; 
        libmesh_error();
      }
  }
  GetPot infile(_system_config_file);
  _rho_L = infile("rho_L", 1.);
  _rho_S = infile("rho_S", 1.);
  _c_L = infile("c_L", 1.);
  _c_S = infile("c_S", 1.);
  _k_L = infile("k_L", 1.);
  _k_S = infile("k_S", 1.);
  _r_L = infile("r_L", 0.);
  _r_S = infile("r_S", 0.);
  _D_c = infile("D_c", 1.);
  _mu = infile("mu", 1.);
  _evaluate_T = infile("evaluate_T", true);
  _pressure_penalty = infile("pressure_penalty", 1e9);

  const int pp_components = infile.vector_variable_size("pressure_point");
  for (unsigned int i=0; i != std::max(pp_components,3); ++i)
    _pressure_point(i) = infile("pressure_point", 0., i);

  const unsigned int dim = this->get_mesh().mesh_dimension();

  // Add the pressure variable "p". This will
  // be approximated with a first-order basis,
  // providing an LBB-stable pressure-velocity pair.
  // Add the velocity components "u" & "v".  They
  // will be approximated using second-order approximation.
  u_var = this->add_variable ("u", static_cast<Order>(_fe_order+1),
                              Utility::string_to_enum<FEFamily>(_fe_family));
  if (dim > 1)
    v_var = this->add_variable ("v", static_cast<Order>(_fe_order+1),
                                Utility::string_to_enum<FEFamily>(_fe_family));
  else
    v_var = 0;

  if (dim > 2)
    w_var = this->add_variable ("w", static_cast<Order>(_fe_order+1),
                                Utility::string_to_enum<FEFamily>(_fe_family));
  else
    w_var = 0;

  p_var = this->add_variable ("p", static_cast<Order>(_fe_order),
                              Utility::string_to_enum<FEFamily>(_fe_family));

  if (_evaluate_T)
    T_var = this->add_variable ("T", static_cast<Order>(_fe_order+1),
                                Utility::string_to_enum<FEFamily>(_fe_family));
  else
    T_var = 0;

  // Do the parent's initialization after variables are defined
  FEMSystem::init_data();

  // Tell the system to march velocity forward in time, but 
  // leave p as a constraint only
  this->time_evolving(u_var);
  if (dim > 1)
    this->time_evolving(v_var);
  if (dim > 2)
    this->time_evolving(w_var);
  if (_evaluate_T)
    this->time_evolving(T_var);
}



void BioheatSystem::init_context (DiffContext &context)
{
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // To enable FE optimizations, we should prerequest all the data
  // we will need to build the linear system.
  c.element_fe_var[u_var]->get_JxW();
  c.element_fe_var[u_var]->get_xyz();
  c.element_fe_var[u_var]->get_phi();
  c.element_fe_var[u_var]->get_dphi();

  c.element_fe_var[p_var]->get_dphi();
  c.element_fe_var[p_var]->get_phi();

  c.side_fe_var[u_var]->get_JxW();
  c.side_fe_var[u_var]->get_xyz();
  c.side_fe_var[u_var]->get_phi();

  c.side_fe_var[p_var]->get_JxW();
  c.side_fe_var[p_var]->get_phi();

  if (_evaluate_T)
    {
      c.element_fe_var[T_var]->get_phi();
      c.element_fe_var[T_var]->get_dphi();
      c.side_fe_var[T_var]->get_phi();
    }
}


bool BioheatSystem::element_time_derivative (bool request_jacobian,
                                            DiffContext &context)
{
  const bool compute_jacobians = request_jacobian && _analytic_jacobians;

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.

  // Subdomain id (for liquid fraction field)
  const subdomain_id_type sbdid = c.elem->subdomain_id();

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real> &JxW = c.element_fe_var[u_var]->get_JxW();

  // Physical location of quadrature points
  const std::vector<Point> &xyz = c.element_fe_var[u_var]->get_xyz();

  // The velocity shape functions at interior quadrature points.
  const std::vector<std::vector<Real> >& phi = c.element_fe_var[u_var]->get_phi();

  // The velocity shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi =
    c.element_fe_var[u_var]->get_dphi();

  // The temperature shape functions at interior
  // quadrature points.
  const std::vector<std::vector<Real> >& Tphi =
    c.element_fe_var[T_var]->get_phi();

  // The pressure shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<Real> >& psi =
    c.element_fe_var[p_var]->get_phi();

//  const std::vector<std::vector<RealGradient> >& dpsi =
//    c.element_fe_var[p_var]->get_dphi();

  // The temperature shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<RealGradient> >& dTphi =
    c.element_fe_var[T_var]->get_dphi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size(); 
  libmesh_assert (n_u_dofs == c.dof_indices_var[v_var].size()); 
  libmesh_assert (n_u_dofs == c.dof_indices_var[w_var].size()); 
  const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();
  const unsigned int n_T_dofs = c.dof_indices_var[T_var].size();

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[u_var][u_var];
  DenseSubMatrix<Number> &Kup = *c.elem_subjacobians[u_var][p_var];

  DenseSubVector<Number> &Fu = *c.elem_subresiduals[u_var];

// We're currently optimizing for the dim==3, _evaluate_T==true case.
//  if (_evaluate_T)
//    {
      DenseSubMatrix<Number> &KTu = *c.elem_subjacobians[T_var][u_var];

      DenseSubVector<Number> &FT = *c.elem_subresiduals[T_var];
//    }

//  if (dim > 1)
//    {
      DenseSubMatrix<Number> &Kuv = *c.elem_subjacobians[u_var][v_var];
      DenseSubMatrix<Number> &Kvu = *c.elem_subjacobians[v_var][u_var];
      DenseSubMatrix<Number> &Kvv = *c.elem_subjacobians[v_var][v_var];
      DenseSubMatrix<Number> &Kvp = *c.elem_subjacobians[v_var][p_var];

      DenseSubVector<Number> &Fv = *c.elem_subresiduals[v_var];

//      if (_evaluate_T)
        DenseSubMatrix<Number> &KTv = *c.elem_subjacobians[T_var][v_var];

//      if (dim > 2)
//        {
          DenseSubMatrix<Number> &Kuw = *c.elem_subjacobians[u_var][w_var];
          DenseSubMatrix<Number> &Kvw = *c.elem_subjacobians[v_var][w_var];
          DenseSubMatrix<Number> &Kwu = *c.elem_subjacobians[w_var][u_var];
          DenseSubMatrix<Number> &Kwv = *c.elem_subjacobians[w_var][v_var];
          DenseSubMatrix<Number> &Kww = *c.elem_subjacobians[w_var][w_var];
          DenseSubMatrix<Number> &Kwp = *c.elem_subjacobians[w_var][p_var];

          DenseSubVector<Number> &Fw = *c.elem_subresiduals[w_var];

//          if (_evaluate_T)
            DenseSubMatrix<Number> &KTw = *c.elem_subjacobians[T_var][w_var];
            DenseSubMatrix<Number> &KTT = *c.elem_subjacobians[T_var][T_var];
//        }
//    }

  // Get the liquid fraction on this subdomain
  libmesh_assert(n_L_values[sbdid]);
  FunctionBase<Number> &n_L_func = *n_L_values[sbdid];
      
  // Now we will build the element Jacobian and residual.
  // Constructing the residual requires the solution and its
  // gradient from the previous timestep.  This must be
  // calculated at each quadrature point by summing the
  // solution degree-of-freedom values by the appropriate
  // weight functions.
  unsigned int n_qpoints = c.element_qrule->n_points();

  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the solution & its gradient at the old Newton iterate
      const Number 
p = c.interior_value(p_var, qp),
             u = c.interior_value(u_var, qp),
             v = (dim > 1) ? c.interior_value(v_var, qp) : 0,
             w = (dim > 2) ? c.interior_value(w_var, qp) : 0;
      const Gradient 
               grad_p = c.interior_gradient(p_var, qp),
               grad_u = c.interior_gradient(u_var, qp),
               grad_v = (dim > 1) ? c.interior_gradient(v_var, qp) : 0,
               grad_w = (dim > 2) ? c.interior_gradient(w_var, qp) : 0,
               grad_T = c.interior_gradient(T_var, qp);

      // Definitions for convenience.  It is sometimes simpler to do a
      // dot product if you have the full vector at your disposal.
      NumberVectorValue U     (u, v, w);
      const Number  p_x = grad_p(0), p_y = grad_p(1), p_z = grad_p(2),
                    u_x = grad_u(0), u_y = grad_u(1), u_z = grad_u(2),
                    v_x = grad_v(0), v_y = grad_v(1), v_z = grad_v(2),
                    w_x = grad_w(0), w_y = grad_w(1), w_z = grad_w(2);

      const Number div_U = u_x + v_y + w_z;

      const Number n_L = n_L_func(xyz[qp], time);
      const Number n_S = (1 - n_L);

      const Number K_P = _D_c*_D_c / 180 * (n_L*n_L*n_L) / (1 - n_L) / (1 - n_L);
      const Number k_v_L = n_L*n_L * _mu / K_P;

//      if (_evaluate_T)
//        {
          const Number K_eff = n_L*_k_L + n_S*_k_S;
          const Gradient omega_eff = _rho_L*_c_L*U;
          const Number R_eff = _rho_L*_r_L + _rho_S*_r_S + k_v_L*(U*U);
//        }
          
      // First, an i-loop over the velocity degrees of freedom.
      // We know that n_u_dofs == n_v_dofs etc. so we can compute
      // contributions for all at the same time.
      for (unsigned int i=0; i != n_u_dofs; i++)
        {
          const Real phi_i = phi[i][qp];
          const RealGradient dphi_i = dphi[i][qp];

          Fu(i) += JxW[qp] * (
                     -_rho_L*(U*grad_u)*phi_i         // convection term
//                     - n_L*(p_x*phi_i)                // pressure term
                     + n_L*(p*dphi_i(0))              // integrated pressure term
                     + (2*n_L*_mu/3.)*div_U*dphi_i(0) // diffusion term 1
                     - (n_L*_mu) * (                  // diffusion term 2
                       2*u_x*dphi_i(0)
                       + (u_y+v_x)*dphi_i(1)
                       + (u_z+w_x)*dphi_i(2)
                     )
                     - k_v_L * u * phi_i              // drag term
                   );

          if (dim > 1)
          Fv(i) += JxW[qp] * (
                     -_rho_L*(U*grad_v)*phi_i         // convection term
//                     - n_L*(p_y*phi_i)                // pressure term
                     + n_L*(p*dphi_i(1))              // integrated pressure term
                     + (2*n_L*_mu/3.)*div_U*dphi_i(1) // diffusion term 1
                     - (n_L*_mu) * (                  // diffusion term 2
                       2*v_y*dphi_i(1)
                       + (u_y+v_x)*dphi_i(0)
                       + (v_z+w_y)*dphi_i(2)
                     )
                     - k_v_L * v * phi_i              // drag term
                   );

          if (dim > 2)
          Fw(i) += JxW[qp] * (
                     -_rho_L*(U*grad_w)*phi_i         // convection term
//                     - n_L*(p_z*phi_i)                // pressure term
                     + n_L*(p*dphi_i(2))              // integrated pressure term
                     + (2*n_L*_mu/3.)*div_U*dphi_i(2) // diffusion term 1
                     - (n_L*_mu) * (                  // diffusion term 2
                       2*w_z*dphi_i(2)
                       + (u_z+w_x)*dphi_i(0)
                       + (v_z+w_y)*dphi_i(1)
                     )
                     - k_v_L * w * phi_i              // drag term
                   );

          // Note that the Fp block is identically zero unless we are using
          // some kind of artificial compressibility scheme...

          if (compute_jacobians) {

          // Matrix contributions for the velocity couplings.
            for (unsigned int j=0; j != n_u_dofs; j++)
              {
                const Real phi_j = phi[j][qp];
                const RealGradient dphi_j = dphi[j][qp];

                Kuu(i,j) += JxW[qp] * (
                     - _rho_L*(U*dphi_j+phi_j*u_x)*phi_i  // convection term
                     + (2*n_L*_mu/3.)*dphi_j(0)*dphi_i(0) // diffusion term 1
                     - (n_L*_mu) * (                      // diffusion term 2
                       2*dphi_j(0)*dphi_i(0)
                       + dphi_j(1)*dphi_i(1)
                       + dphi_j(2)*dphi_i(2)
                     )
                     - k_v_L * phi_j * phi_i              // drag term
                   );

                if (dim > 1) {
                Kuv(i,j) += JxW[qp] * (
                     -_rho_L*(phi_j*u_y)*phi_i            // convection term
                     + (2*n_L*_mu/3.)*dphi_j(1)*dphi_i(0) // diffusion term 1
                     - (n_L*_mu) * (                      // diffusion term 2
                       dphi_j(0)*dphi_i(1)
                     )
                   );

                Kvu(i,j) += JxW[qp] * (
                     -_rho_L*(phi_j*v_x)*phi_i            // convection term
                     + (2*n_L*_mu/3.)*dphi_j(0)*dphi_i(1) // diffusion term 1
                     - (n_L*_mu) * (                      // diffusion term 2
                       dphi_j(1)*dphi_i(0)
                     )
                   );

                Kvv(i,j) += JxW[qp] * (
                     -_rho_L*(U*dphi_j+phi_j*v_y)*phi_i   // convection term
                     + (2*n_L*_mu/3.)*dphi_j(1)*dphi_i(1) // diffusion term 1
                     - (n_L*_mu) * (                      // diffusion term 2
                       2*dphi_j(1)*dphi_i(1)
                       + dphi_j(0)*dphi_i(0)
                       + dphi_j(2)*dphi_i(2)
                     )
                     - k_v_L * phi_j * phi_i              // drag term
                   );
                }

                if (dim > 2) {
                  Kwu(i,j) += JxW[qp] * (
                     -_rho_L*(phi_j*w_x)*phi_i            // convection term
                     + (2*n_L*_mu/3.)*dphi_j(0)*dphi_i(2) // diffusion term 1
                     - (n_L*_mu) * (                        // diffusion term 2
                       dphi_j(2)*dphi_i(0)
                     )
                   );

                  Kuw(i,j) += JxW[qp] * (
                     -_rho_L*(phi_j*u_z)*phi_i            // convection term
                     + (2*n_L*_mu/3.)*dphi_j(2)*dphi_i(0) // diffusion term 1
                     - (n_L*_mu) * (                      // diffusion term 2
                       dphi_j(0)*dphi_i(2)
                     )
                   );

                  Kwv(i,j) += JxW[qp] * (
                     -_rho_L*(phi_j*w_y)*phi_i            // convection term
                     + (2*n_L*_mu/3.)*dphi_j(1)*dphi_i(2) // diffusion term 1
                     - (n_L*_mu) * (                      // diffusion term 2
                       dphi_j(2)*dphi_i(1)
                     )
                   );

                  Kvw(i,j) += JxW[qp] * (
                     -_rho_L*(phi_j*v_z)*phi_i            // convection term
                     + (2*n_L*_mu/3.)*dphi_j(2)*dphi_i(1) // diffusion term 1
                     - (n_L*_mu) * (                      // diffusion term 2
                       dphi_j(1)*dphi_i(2)
                     )
                   );

                  Kww(i,j) += JxW[qp] * (
                     -_rho_L*(U*dphi_j+phi_j*w_z)*phi_i   // convection term
                     + (2*n_L*_mu/3.)*dphi_j(2)*dphi_i(2) // diffusion term 1
                     - (n_L*_mu) * (                      // diffusion term 2
                       2*dphi_j(2)*dphi_i(2)
                       + dphi_j(0)*dphi_i(0)
                       + dphi_j(1)*dphi_i(1)
                     )
                     - k_v_L * phi_j * phi_i              // drag term
                   );
                }
              }

          // Matrix contributions for the pressure couplings.
            for (unsigned int j=0; j != n_p_dofs; j++)
              {
                const Real psi_j = psi[j][qp];
//                const RealGradient dpsi_j = dpsi[j][qp];
  
//                Kup(i,j) += JxW[qp]*-n_L*(dpsi_j(0)*phi_i); // pressure term
                Kup(i,j) += JxW[qp]* n_L*(psi_j*dphi_i(0)); // integrated pressure term
                if (dim > 1)
//                Kvp(i,j) += JxW[qp]*-n_L*(dpsi_j(1)*phi_i); // pressure term
                Kvp(i,j) += JxW[qp]* n_L*(psi_j*dphi_i(1)); // integrated pressure term
                if (dim > 2)
//                Kwp(i,j) += JxW[qp]*-n_L*(dpsi_j(2)*phi_i); // pressure term
                Kwp(i,j) += JxW[qp]* n_L*(psi_j*dphi_i(2)); // integrated pressure term
              }
          }
        }

      if (_evaluate_T)
        for (unsigned int i=0; i != n_T_dofs; i++)
          {
            const Real Tphi_i = Tphi[i][qp];
            const RealGradient dTphi_i = dTphi[i][qp];

            // Matrix contributions for the temperature couplings.
            FT(i) += JxW[qp] * (
                       -K_eff*(grad_T*dTphi_i) // diffusion term
                       - (omega_eff*grad_T     // convection term
                          - R_eff)*Tphi_i);    // production term

            if (compute_jacobians)
              {
                for (unsigned int j=0; j != n_T_dofs; j++)
                  {
                    const RealGradient dTphi_j = dTphi[j][qp];

                    KTT(i,j) += JxW[qp] * (
                         -K_eff*(dTphi_j*dTphi_i)       // diffusion term
                         - (omega_eff*dTphi_j)*Tphi_i); // convection term
                  }

                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    const Real phi_j = phi[j][qp];

                    KTu(i,j) += JxW[qp] * (
                         -(_rho_L*_c_L*phi_j*grad_T(0) // convection term
                           - k_v_L*2*u*phi_j)*Tphi_i); // production term

                    if (dim > 1)
                    KTv(i,j) += JxW[qp] * (
                         -(_rho_L*_c_L*phi_j*grad_T(1) // convection term
                           - k_v_L*2*v*phi_j)*Tphi_i); // production term


                    if (dim > 2)
                    KTw(i,j) += JxW[qp] * (
                         -(_rho_L*_c_L*phi_j*grad_T(2) // convection term
                           - k_v_L*2*w*phi_j)*Tphi_i); // production term
                  }
              }
          }
    } // end of the quadrature point qp-loop
  
  return compute_jacobians;
}



bool BioheatSystem::element_constraint (bool request_jacobian,
                                       DiffContext &context)
{
  const bool compute_jacobians = request_jacobian && _analytic_jacobians;

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.

  // Subdomain id (for liquid fraction field)
  const subdomain_id_type sbdid = c.elem->subdomain_id();

  // Element Jacobian * quadrature weight for interior integration
  const std::vector<Real> &JxW = c.element_fe_var[u_var]->get_JxW();

  // Physical location of quadrature points
  const std::vector<Point> &xyz = c.element_fe_var[u_var]->get_xyz();

  // The velocity shape functions at interior quadrature points.
//  const std::vector<std::vector<Real> >& phi =
//    c.element_fe_var[u_var]->get_phi();

const std::vector<std::vector<RealGradient> >& dphi =
  c.element_fe_var[u_var]->get_dphi();

  // The pressure shape function gradients at interior
  // quadrature points.
  const std::vector<std::vector<Real> >& psi =
    c.element_fe_var[p_var]->get_phi();

  const std::vector<std::vector<RealGradient> >& dpsi =
    c.element_fe_var[p_var]->get_dphi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
  const unsigned int n_p_dofs = c.dof_indices_var[p_var].size();

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kpu = *c.elem_subjacobians[p_var][u_var];
  DenseSubMatrix<Number> &Kpv = *c.elem_subjacobians[p_var][v_var];
  DenseSubMatrix<Number> &Kpw = *c.elem_subjacobians[p_var][w_var];
  DenseSubVector<Number> &Fp = *c.elem_subresiduals[p_var];

  libmesh_assert(n_L_values[sbdid]);
  FunctionBase<Number> &n_L_func = *n_L_values[sbdid];
      
  // Add the constraint given by the continuity equation
  unsigned int n_qpoints = c.element_qrule->n_points();
  for (unsigned int qp=0; qp != n_qpoints; qp++)
    {
      // Compute the velocity at the old Newton iterate
//      const Number
//             u = c.interior_value(u_var, qp),
//             v = (dim > 1) ? c.interior_value(v_var, qp) : 0,
//             w = (dim > 2) ? c.interior_value(w_var, qp) : 0;

      const Gradient grad_u = c.interior_gradient(u_var, qp),
                     grad_v = c.interior_gradient(v_var, qp),
                     grad_w = c.interior_gradient(w_var, qp);
      const Real div_u = grad_u(0) + grad_v(1) + grad_w(2);

      const Number n_L = n_L_func(xyz[qp], time);
      // const Number n_S = (1 - n_L);

      // Now a loop over the pressure degrees of freedom.  This
      // computes the contributions of the continuity equation.
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
//          const RealGradient dpsi_i = dpsi[i][qp];
          const Real psi_i = psi[i][qp];

//          Fp(i) += JxW[qp]*n_L *
//                   (u*dpsi_i(0) + v*dpsi_i(1) + w*dpsi_i(2));
          Fp(i) += JxW[qp]*n_L *
                   div_u*psi_i;

          if (compute_jacobians)
            for (unsigned int j=0; j != n_u_dofs; j++)
              {
//                const Real phi_j = phi[j][qp];
                const RealGradient dphi_j = dphi[j][qp];

//                Kpu(i,j) += JxW[qp]*n_L*phi_j*dpsi_i(0);
                Kpu(i,j) += JxW[qp]*n_L*dphi_j(0)*psi_i;
                if (dim > 1)
//                  Kpv(i,j) += JxW[qp]*n_L*phi_j*dpsi_i(1);
                  Kpv(i,j) += JxW[qp]*n_L*dphi_j(1)*psi_i;
                if (dim > 2)
//                  Kpw(i,j) += JxW[qp]*n_L*phi_j*dpsi_i(2);
                  Kpw(i,j) += JxW[qp]*n_L*dphi_j(2)*psi_i;
              }
        }
    } // end of the quadrature point qp-loop

  // Pin p = 0 if requested
  if ((_pressure_penalty != 0.) &&
      c.elem->contains_point(_pressure_point))
    {
      DenseSubMatrix<Number> &Kpp = *c.elem_subjacobians[p_var][p_var];
      DenseSubVector<Number> &Fp = *c.elem_subresiduals[p_var];
      const unsigned int n_p_dofs = c.dof_indices_var[p_var].size(); 

      Point zero(0.);
      Number p = c.point_value(p_var, zero);
      Number p_value = 0.;

      unsigned int dim = get_mesh().mesh_dimension();
      FEType fe_type = c.element_fe_var[p_var]->get_fe_type();
      Point p_master = FEInterface::inverse_map(dim, fe_type, c.elem, zero);

      std::vector<Real> point_phi(n_p_dofs);
      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          point_phi[i] = FEInterface::shape(dim, fe_type, c.elem, i, p_master);
        }

      for (unsigned int i=0; i != n_p_dofs; i++)
        {
          Fp(i) += _pressure_penalty * (p - p_value) * point_phi[i];
          if (compute_jacobians)
            for (unsigned int j=0; j != n_p_dofs; j++)
              {
		Kpp(i,j) += _pressure_penalty * point_phi[i] * point_phi[j];
              }
        }
    }

  return compute_jacobians;
}

      

bool BioheatSystem::side_time_derivative (bool request_jacobian,
                                          DiffContext &context)
{
  const bool compute_jacobians = request_jacobian && _analytic_jacobians;

  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // See if we're on a Neumann boundary
  boundary_id_type boundary_id =
    this->get_mesh().boundary_info->boundary_id(c.elem, c.side);
  libmesh_assert (boundary_id != BoundaryInfo::invalid_id);

  const bool neumann_u =
    neumann_values.count(boundary_id) &&
    (std::find(neumann_variables[boundary_id].begin(),
               neumann_variables[boundary_id].end(),
               u_var) != 
               neumann_variables[boundary_id].end());

  const bool neumann_v =
    neumann_values.count(boundary_id) &&
    (std::find(neumann_variables[boundary_id].begin(),
               neumann_variables[boundary_id].end(),
               v_var) != 
               neumann_variables[boundary_id].end());

  const bool neumann_w =
    neumann_values.count(boundary_id) &&
    (std::find(neumann_variables[boundary_id].begin(),
               neumann_variables[boundary_id].end(),
               w_var) != 
               neumann_variables[boundary_id].end());

  const bool neumann_T = _evaluate_T &&
    neumann_values.count(boundary_id) &&
    (std::find(neumann_variables[boundary_id].begin(),
               neumann_variables[boundary_id].end(),
               T_var) != 
               neumann_variables[boundary_id].end());

  // If we're not on a heterogenous Neumann boundary we have nothing
  // to do
  if (!(neumann_u || neumann_v || neumann_w || neumann_T))
    return compute_jacobians;

  // We're on a Neumann boundary so we need boundary values
  libmesh_assert(neumann_values.count(boundary_id));

  // Subdomain id (for liquid fraction field)
  const subdomain_id_type sbdid = c.elem->subdomain_id();

  // Element Jacobian * quadrature weight for side integration
  const std::vector<Real> &JxW_side = c.side_fe_var[u_var]->get_JxW();

  // Physical location of side quadrature points
  const std::vector<Point> &xyz_side = c.side_fe_var[u_var]->get_xyz();

  // The velocity shape functions at side quadrature points.
  const std::vector<std::vector<Real> >& phi_side =
    c.side_fe_var[u_var]->get_phi();

  // The temperature shape functions at side quadrature points.
  const std::vector<std::vector<Real> >& Tphi_side =
    c.side_fe_var[T_var]->get_phi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
  const unsigned int n_T_dofs = c.dof_indices_var[T_var].size();

  // The subvectors we need to fill:
  // (no submatrices - linear Neumann conditions look like constant
  // forcing functions)
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubVector<Number> &Fu = *c.elem_subresiduals[u_var];
  DenseSubVector<Number> &Fv = *c.elem_subresiduals[v_var];
  DenseSubVector<Number> &Fw = *c.elem_subresiduals[w_var];
  DenseSubVector<Number> &FT = *c.elem_subresiduals[T_var];

  // Get the liquid fraction on this subdomain
  libmesh_assert(n_L_values[sbdid]);
  FunctionBase<Number> &n_L_func = *n_L_values[sbdid];
      
  unsigned int n_sidepoints = c.side_qrule->n_points();
  for (unsigned int qp=0; qp != n_sidepoints; qp++)
    {
      const Number n_L = n_L_func(xyz_side[qp], time);
      const Number n_S = (1 - n_L);

      // One function value for each vector component, one for pressure, 
      DenseVector<Number> neumann(dim+1+_evaluate_T);
      (*neumann_values[boundary_id])(xyz_side[qp],time,neumann);

      // if (_evaluate_T)
      const Number K_eff = n_L*_k_L + n_S*_k_S;

      if (neumann_u || neumann_v || neumann_w)
        {
          const Number u_flux = neumann(0),
                       v_flux = dim>1?neumann(1):0,
                       w_flux = dim>2?neumann(2):0;

	  // We know that n_u_dofs == n_v_dofs etc. so we can compute
	  // contributions for all at the same time.
          for (unsigned int i=0; i != n_u_dofs; i++)
            {
              const Real phi_i = phi_side[i][qp];

              if (neumann_u)
                Fu(i) += JxW_side[qp] * u_flux * phi_i;
              if (neumann_v && dim > 1)
                Fv(i) += JxW_side[qp] * v_flux * phi_i;
              if (neumann_w && dim > 2)
                Fw(i) += JxW_side[qp] * w_flux * phi_i;

              // No jacobian terms for Neumann conditions
            }
        }
      if (neumann_T)
        {
          const Number T_grad = neumann(dim+1);

          for (unsigned int i=0; i != n_T_dofs; i++)
            {
              const Real Tphi_i = Tphi_side[i][qp];

              FT(i) += JxW_side[qp] * K_eff * T_grad * Tphi_i;
            }
        }
    }

  return compute_jacobians;
}

bool BioheatSystem::side_constraint (bool request_jacobian,
                                    DiffContext &context)
{
  const bool compute_jacobians = request_jacobian && _analytic_jacobians;

/*
  FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

  // See if we're on a Dirichlet boundary for anything
  boundary_id_type boundary_id =
    this->get_mesh().boundary_info->boundary_id(c.elem, c.side);
  libmesh_assert (boundary_id != BoundaryInfo::invalid_id);

  // If we're not on a Dirichlet boundary we have nothing to do
  if (!dirichlet_values.count(boundary_id))
    return compute_jacobians;

  const bool dirichlet_u =
    dirichlet_values.count(boundary_id) &&
    (std::find(dirichlet_variables[boundary_id].begin(),
               dirichlet_variables[boundary_id].end(),
               u_var) != 
               dirichlet_variables[boundary_id].end());

  const bool dirichlet_v =
    dirichlet_values.count(boundary_id) &&
    (std::find(dirichlet_variables[boundary_id].begin(),
               dirichlet_variables[boundary_id].end(),
               v_var) != 
               dirichlet_variables[boundary_id].end());

  const bool dirichlet_w =
    dirichlet_values.count(boundary_id) &&
    (std::find(dirichlet_variables[boundary_id].begin(),
               dirichlet_variables[boundary_id].end(),
               w_var) != 
               dirichlet_variables[boundary_id].end());

  const bool dirichlet_p =
    dirichlet_values.count(boundary_id) &&
    (std::find(dirichlet_variables[boundary_id].begin(),
               dirichlet_variables[boundary_id].end(),
               p_var) != 
               dirichlet_variables[boundary_id].end());

  const bool dirichlet_T =
    _evaluate_T &&
    dirichlet_values.count(boundary_id) &&
    (std::find(dirichlet_variables[boundary_id].begin(),
               dirichlet_variables[boundary_id].end(),
               w_var) != 
               dirichlet_variables[boundary_id].end());

  // We're on a Dirichlet boundary it's got to be Dirichlet for
  // *something*
  libmesh_assert (dirichlet_u || dirichlet_v || dirichlet_w || 
                  dirichlet_T || dirichlet_p);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.

  // Element Jacobian * quadrature weight for side integration
  const std::vector<Real> &JxW_side = c.side_fe_var[u_var]->get_JxW();

  // Physical location of side quadrature points
  const std::vector<Point> &xyz_side = c.side_fe_var[u_var]->get_xyz();

  // Current time
  const Real time = c.time;

  // The velocity, pressure, and temperature shape functions at side
  // quadrature points.
  const std::vector<std::vector<Real> >& phi_side =
    c.side_fe_var[u_var]->get_phi();

  const std::vector<std::vector<Real> >& pphi_side =
    c.side_fe_var[p_var]->get_phi();

  const std::vector<std::vector<Real> >& Tphi_side =
    c.side_fe_var[T_var]->get_phi();

  // The number of local degrees of freedom in each of u, v and w
  const unsigned int n_u_dofs = c.dof_indices_var[u_var].size(); 

  // The number of local degrees of freedom in T
  const unsigned int n_T_dofs = c.dof_indices_var[T_var].size(); 

  // The number of local degrees of freedom in p
  const unsigned int n_p_dofs = c.dof_indices_var[p_var].size(); 

  // The subvectors and submatrices we need to fill:
  const unsigned int dim = this->get_mesh().mesh_dimension();
  DenseSubMatrix<Number> &Kuu = *c.elem_subjacobians[u_var][u_var];
  DenseSubMatrix<Number> &Kvv = *c.elem_subjacobians[v_var][v_var];
  DenseSubMatrix<Number> &Kww = *c.elem_subjacobians[w_var][w_var];
  DenseSubMatrix<Number> &Kpp = *c.elem_subjacobians[p_var][p_var];
  DenseSubMatrix<Number> &KTT = *c.elem_subjacobians[T_var][T_var];

  DenseSubVector<Number> &Fu = *c.elem_subresiduals[u_var];
  DenseSubVector<Number> &Fv = *c.elem_subresiduals[v_var];
  DenseSubVector<Number> &Fw = *c.elem_subresiduals[w_var];
  DenseSubVector<Number> &Fp = *c.elem_subresiduals[p_var];
  DenseSubVector<Number> &FT = *c.elem_subresiduals[T_var];

  // Dirichlet velocity boundary conditions imposed at each timestep
  // via the penalty method.

// TODO - configurable penalty

  // The penalty value.  \f$ \frac{1}{\epsilon} \f$
  const Real penalty = 1.e10;

  unsigned int n_sidepoints = c.side_qrule->n_points();

  for (unsigned int qp=0; qp != n_sidepoints; qp++)
    {
      // One function value for each vector component, one for pressure, 
      DenseVector<Number> dirichlet(dim+1+_evaluate_T);
      (*dirichlet_values[boundary_id])(xyz_side[qp],time,dirichlet);

      const Real JxWxPen = JxW_side[qp] * penalty;

      if (dirichlet_u || dirichlet_v || dirichlet_w)
        {
          Number u = c.side_value(u_var, qp),
                 v = c.side_value(v_var, qp),
                 w = c.side_value(w_var, qp);

          const Number u_value = dirichlet(0),
                       v_value = dim>1?dirichlet(1):0,
                       w_value = dim>2?dirichlet(2):0;

          for (unsigned int i=0; i != n_u_dofs; i++)
            {
              const Real phi_i = phi_side[i][qp];

              if (dirichlet_u)
                Fu(i) += JxWxPen * (u - u_value) * phi_i;
              if (dirichlet_v && dim > 1)
                Fv(i) += JxWxPen * (v - v_value) * phi_i;
              if (dirichlet_w && dim > 2)
                Fw(i) += JxWxPen * (w - w_value) * phi_i;

              if (compute_jacobians)
                for (unsigned int j=0; j != n_u_dofs; j++)
                  {
                    const Real phi_j = phi_side[j][qp];

                    if (dirichlet_u)
                      Kuu(i,j) += JxWxPen * phi_i * phi_j;
                    if (dirichlet_v && dim > 1)
                      Kvv(i,j) += JxWxPen * phi_i * phi_j;
                    if (dirichlet_w && dim > 2)
                      Kww(i,j) += JxWxPen * phi_i * phi_j;
                  }
            }
        }

      if (dirichlet_p)
        {
          Number p = c.side_value(p_var, qp);

          const Number p_value = dirichlet(dim);

          for (unsigned int i=0; i != n_p_dofs; i++)
            {
              const Real pphi_i = pphi_side[i][qp];

              FT(i) += JxWxPen * (p - p_value) * pphi_i;

              if (compute_jacobians)
                for (unsigned int j=0; j != n_p_dofs; j++)
                  {
                    const Real pphi_j = pphi_side[j][qp];

                    Kpp(i,j) += JxWxPen * pphi_i * pphi_j;
                  }
            }
        }

      if (dirichlet_T)
        {
          Number T = c.side_value(T_var, qp);

          const Number T_value = dirichlet(dim+1);

          for (unsigned int i=0; i != n_T_dofs; i++)
            {
              const Real Tphi_i = Tphi_side[i][qp];

              FT(i) += JxWxPen * (T - T_value) * Tphi_i;

              if (compute_jacobians)
                for (unsigned int j=0; j != n_T_dofs; j++)
                  {
                    const Real Tphi_j = Tphi_side[j][qp];

                    KTT(i,j) += JxWxPen * Tphi_i * Tphi_j;
                  }
            }
        }
    }
*/

  return compute_jacobians;
}
