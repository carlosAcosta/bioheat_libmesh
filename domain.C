
#include "domain.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "mesh_modification.h"
#include "mesh_refinement.h"
#include "mesh_triangle_interface.h"

// Headers for OD2 construction
#include "boundary_info.h"
#include "elem.h"
#include "face_tri3.h"
#include "face_quad4.h"

// Headers for cylinder construction
#include "point_locator_tree.h"

// Useful generation function for wave problems
void build_od2_square (UnstructuredMesh &mesh,
                       unsigned int nx,
                       unsigned int ny,
                       Real xmin=0., Real xmax=1.,
                       Real ymin=0., Real ymax=1.,
                       const ElemType type=INVALID_ELEM);

// Mesh construction
void build_domain (UnstructuredMesh &mesh, FEMParameters &param)
{
  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the square [-1,1]^D.  We instruct the mesh generator
  // to build a mesh of 8x8 \p Quad4 elements in 2D, or \p Hex8
  // elements in 3D.  HERMITE elements only need first order
  // nodes

  ElemType elemtype;
  if (param.dimension == 1)
    {
      elemtype = EDGE2;
    }
  else if (param.dimension == 2)
    {
      if (param.elementtype == "tri" ||
          param.elementtype == "unstructured")
        elemtype = TRI3;
      else
        elemtype = QUAD4;
    }
  else
    {
      libmesh_assert(param.dimension == 3);
      if (param.elementtype == "tri")
        elemtype = TET4;
      else
        elemtype = HEX8;
    }
  mesh.set_mesh_dimension(param.dimension);
 
  if (param.domaintype == "sphere")
    {
      if (param.dimension == 2)
        {
          if (param.elementtype == "unstructured")
            {
/*
	      // Try to match the desired per-triangle area of the
              // structured circle
              const Real desired_area =
                M_PI / 40.0 * std::pow(0.5, param.coarsegridx);

              const Real desired_h_boundary = std::sqrt(desired_area);

              const Real h_boundary = 
                std::cos(libMesh::pi*(0.5- 1./static_cast<Real>(n_points)));
              std::cout << "h_boundary=" << h_boundary << std::endl;
            
              // Add the points which make up the circle
              libmesh_assert (n_points > 0);
              const Real
                dtheta = 2*libMesh::pi / static_cast<Real>(n_points);

              for (unsigned int i=0; i<n_points; ++i)
                {
                  mesh.add_point( Point(center(0)+ radius*cos(i*dtheta),
                                        center(1)+ radius*sin(i*dtheta)) );
                }
*/


            }
          else
            MeshTools::Generation::build_sphere
              (mesh, 0.5, param.coarsegridx, elemtype);

          if (param.domain_edge_width != 1.0 ||
              param.domain_edge_length != 1.0)
            MeshTools::Modification::scale
              (mesh, param.domain_edge_width,
                     param.domain_edge_length);
        }
      else
        libmesh_error();
    }

  else if (param.domaintype == "square" || param.domaintype == "tube")
    {
      if (param.dimension == 1)
        MeshTools::Generation::build_line
          (mesh, param.coarsegridx,
           param.domain_xmin, param.domain_xmin + param.domain_edge_width,
           EDGE2);
      else if (param.dimension == 2)
        {
          if (param.elementtype == "unstructured")
            {
              MeshTools::Generation::build_delaunay_square
                (mesh, param.coarsegridx, param.coarsegridy,
                 param.domain_xmin, param.domain_xmin + param.domain_edge_width,
                 param.domain_ymin, param.domain_ymin + param.domain_edge_length,
                 elemtype);
            }
          else
            MeshTools::Generation::build_square
              (mesh, param.coarsegridx, param.coarsegridy,
               param.domain_xmin, param.domain_xmin + param.domain_edge_width,
               param.domain_ymin, param.domain_ymin + param.domain_edge_length,
               elemtype);
        }
      else if (param.dimension == 3)
        {
          MeshTools::Generation::build_cube
            (mesh, param.coarsegridx, param.coarsegridy, param.coarsegridz,
             param.domain_xmin, param.domain_xmin + param.domain_edge_width,
             param.domain_ymin, param.domain_ymin + param.domain_edge_length,
             param.domain_zmin, param.domain_zmin + param.domain_edge_height,
             elemtype);
        }
      else
        libmesh_error();
    }
  else if (param.domaintype == "file")
    {
      mesh.read(param.domainfile);

      if (param.domain_xmin != 0.0 ||
          param.domain_ymin != 0.0 ||
          param.domain_zmin != 0.0 ||
          param.domain_edge_width != 1.0 ||
          param.domain_edge_length != 1.0 ||
          param.domain_edge_height != 1.0)
        {
          std::cerr << "Mesh file scaling is not yet supported"
                    << std::endl;
          libmesh_error();
        }

// FIXME - all_tri currently invalidates boundary ids on triangle/hybrid meshes
//      if (param.elementtype == "tri")
//        MeshTools::Modification::all_tri(mesh);
    }
  else if (param.domaintype == "od2")
    {
      libmesh_assert (param.dimension == 2);

      build_od2_square
        (mesh, param.coarsegridx, param.coarsegridy,
         param.domain_xmin, param.domain_xmin + param.domain_edge_width,
         param.domain_ymin, param.domain_ymin + param.domain_edge_length,
         elemtype);
    }
  else if (param.domaintype == "cylinder")
    {
      libmesh_assert (param.dimension == 3);

      Mesh circle;
      circle.set_mesh_dimension(2);

      MeshTools::Generation::build_sphere
        (circle, 0.5, param.coarsegridx, QUAD4);
      MeshTools::Generation::build_extrusion
        (mesh, circle, param.coarsegridz, RealVectorValue(0,0,1));

      if (param.domain_edge_width != 1.0 ||
          param.domain_edge_length != 1.0 ||
          param.domain_edge_height != 1.0)
        MeshTools::Modification::scale
          (mesh, param.domain_edge_width,
                 param.domain_edge_length,
                 param.domain_edge_height);
    }
  else if (param.domaintype == "tube")
    {
      if (
//          param.domain_xmin != 0.0 ||
          param.domain_ymin != 0.0 ||
//          param.domain_zmin != 0.0 || 
//          param.domain_edge_width != 1.0 ||
          param.domain_edge_length != 1.0)
//          param.domain_edge_height != 1.0)
        {
          std::cerr << "Theta must range from 0 to 1"
                    << std::endl;
          libmesh_error();
        }

      PointLocatorTree point_locator(mesh);

      MeshBase::element_iterator       el     = mesh.elements_begin();
      const MeshBase::element_iterator end_el = mesh.elements_end();

      for (; el != end_el; ++el)
        {
          Elem* elem = *el;

          // The boundary we sew up causes problems with BoundaryInfo
          // later, so for now remove all boundary ids
          mesh.boundary_info->remove(elem);

          for (unsigned int n = 0; n != elem->n_nodes(); ++n)
            {
              Node &node = *(elem->get_node(n));
              if (std::abs(node(1) - param.domain_edge_length) < TOLERANCE)
                {
                  Point wrapped_point = node;
                  wrapped_point(1) = 0;
                  Elem *wrapped_elem =
                    const_cast<Elem *>(point_locator(wrapped_point));
                  libmesh_assert(elem);

                  bool found_matching_point = false;
                  for (unsigned int wn = 0; wn != wrapped_elem->n_nodes();
                       ++wn)
                    {
                      Node *new_node = wrapped_elem->get_node(wn);
                      if (new_node->absolute_fuzzy_equals(wrapped_point))
                        {
                          elem->set_node(n) = new_node;
                          found_matching_point = true;
                          break;
                        }
                    }
                  libmesh_assert(found_matching_point);
                }
            }
        }

      MeshBase::node_iterator       nd     = mesh.nodes_begin();
      const MeshBase::node_iterator end_nd = mesh.nodes_end();
      if (param.dimension == 3 || param.dimension == 2)
        {
          for (; nd != end_nd; ++nd)
            {
              Node &node = **nd;
              Real r = node(0), theta = node(1); //, z = node(2);
              node(0) = r * cos(2.0*M_PI*theta);
              node(1) = r * sin(2.0*M_PI*theta);
              // node(2) = z;

	      // The boundary we sew up causes problems with
	      // BoundaryInfo later, so for now remove all boundary
	      // ids
              mesh.boundary_info->remove(&node);

            }
        }
      else
        {
          libmesh_error();
        }
    }
  else
    {
      libmesh_error();
    }

  if (param.elementorder == 2)
    mesh.all_second_order();

  MeshRefinement mesh_refinement(mesh);
  mesh_refinement.uniformly_refine(param.coarserefinements);
}


void build_od2_square (UnstructuredMesh &mesh,
                       unsigned int nx,
                       unsigned int ny,
                       Real xmin, Real xmax,
                       Real ymin, Real ymax,
                       const ElemType type)
{
  libmesh_assert (type == QUAD4 || type == TRI3);

  // Start with a quad mesh
  MeshTools::Generation::build_square (mesh,
                                       nx,
                                       ny,
                                       xmin, xmax,
                                       ymin, ymax,
                                       QUAD4);

  // Begin generating OD2 triangle mesh

  // Loop over the elements, adding the centroid of each to the mesh.
  // Store pointers
  std::vector<Node*> centroids;
  {
    MeshBase::element_iterator       el     = mesh.elements_begin();
    const MeshBase::element_iterator end_el = mesh.elements_end();

    for (; el != end_el; ++el)
      {
        Elem* elem = *el;

        centroids.push_back ( mesh.add_point(elem->centroid()) );
      }
  }

  // Loop over the elements again, this time create the triangles,
  // save them in a vector
  std::vector<Elem*> triangles;

  // BoundaryInfo Storage for element ids, sides, and BC ids
  std::vector<Elem*>              saved_boundary_elements;
  std::vector<short int>          saved_bc_ids;
  std::vector<unsigned short int> saved_bc_sides;
  {
    MeshBase::element_iterator       el     = mesh.elements_begin();
    const MeshBase::element_iterator end_el = mesh.elements_end();

    for (; el != end_el; ++el)
      {
        Elem* elem = *el;

        // Examine left and right neighbors
        const bool left_edge  = (elem->neighbor(3)==NULL);
        const bool right_edge = (elem->neighbor(1)==NULL);

        // Examine top and bottom neighbors
        const bool top_edge = (elem->neighbor(2)==NULL);
        const bool bot_edge = (elem->neighbor(0)==NULL);


        // Can't make OD2 mesh out of 1 element
        libmesh_assert (!(left_edge && right_edge));

        if (left_edge)
          {
            // Build 5 triangles
            Tri3* tri1 = new Tri3;
            tri1->set_node(0) = elem->get_node(0);
            tri1->set_node(1) = centroids[elem->id()];
            tri1->set_node(2) = elem->get_node(3);
            triangles.push_back(tri1);

            // Also save BC data
            saved_boundary_elements.push_back(tri1);
            saved_bc_sides.push_back(2);
            saved_bc_ids.push_back(3);

            Tri3* tri2 = new Tri3;
            tri2->set_node(0) = elem->get_node(0);
            tri2->set_node(1) = elem->get_node(1);
            tri2->set_node(2) = centroids[elem->id()];
            triangles.push_back(tri2);

            if (bot_edge)
              {
                saved_boundary_elements.push_back(tri2);
                saved_bc_sides.push_back(0);
                saved_bc_ids.push_back(0);
              }

            Tri3* tri3 = new Tri3;
            tri3->set_node(0) = elem->get_node(1);
            tri3->set_node(1) = centroids[elem->neighbor(1)->id()];
            tri3->set_node(2) = centroids[elem->id()];
            triangles.push_back(tri3);

            Tri3* tri4 = new Tri3;
            tri4->set_node(0) = centroids[elem->id()];
            tri4->set_node(1) = centroids[elem->neighbor(1)->id()];
            tri4->set_node(2) = elem->get_node(2);
            triangles.push_back(tri4);

            Tri3* tri5 = new Tri3;
            tri5->set_node(0) = centroids[elem->id()];
            tri5->set_node(1) = elem->get_node(2);
            tri5->set_node(2) = elem->get_node(3);
            triangles.push_back(tri5);

            if (top_edge)
              {
                saved_boundary_elements.push_back(tri5);
                saved_bc_sides.push_back(1);
                saved_bc_ids.push_back(2);
              }
          }

        else if (right_edge)
          {
            // Build 3 triangles
            Tri3* tri1 = new Tri3;
            tri1->set_node(0) = centroids[elem->id()];
            tri1->set_node(1) = elem->get_node(2);
            tri1->set_node(2) = elem->get_node(3);
            triangles.push_back(tri1);

            if (top_edge)
              {
                saved_boundary_elements.push_back(tri1);
                saved_bc_sides.push_back(1);
                saved_bc_ids.push_back(2);
              }

            Tri3* tri2 = new Tri3;
            tri2->set_node(0) = elem->get_node(0);
            tri2->set_node(1) = elem->get_node(1);
            tri2->set_node(2) = centroids[elem->id()];
            triangles.push_back(tri2);

            if (bot_edge)
              {
                saved_boundary_elements.push_back(tri2);
                saved_bc_sides.push_back(0);
                saved_bc_ids.push_back(0);
              }

            Tri3* tri3 = new Tri3;
            tri3->set_node(0) = centroids[elem->id()];
            tri3->set_node(1) = elem->get_node(1);
            tri3->set_node(2) = elem->get_node(2);
            triangles.push_back(tri3);


            saved_boundary_elements.push_back(tri3);
            saved_bc_sides.push_back(1);
            saved_bc_ids.push_back(1);
          }

        else if (!left_edge && !right_edge)
          {
            // Build 4 triangles
            Tri3* tri1 = new Tri3;
            tri1->set_node(0) = centroids[elem->id()];
            tri1->set_node(1) = elem->get_node(2);
            tri1->set_node(2) = elem->get_node(3);
            triangles.push_back(tri1);

            if (top_edge)
              {
                saved_boundary_elements.push_back(tri1);
                saved_bc_sides.push_back(1);
                saved_bc_ids.push_back(2);
              }

            Tri3* tri2 = new Tri3;
            tri2->set_node(0) = elem->get_node(0);
            tri2->set_node(1) = elem->get_node(1);
            tri2->set_node(2) = centroids[elem->id()];
            triangles.push_back(tri2);

            if (bot_edge)
              {
                saved_boundary_elements.push_back(tri2);
                saved_bc_sides.push_back(0);
                saved_bc_ids.push_back(0);
              }

            Tri3* tri3 = new Tri3;
            tri3->set_node(0) = elem->get_node(1);
            tri3->set_node(1) = centroids[elem->neighbor(1)->id()];
            tri3->set_node(2) = centroids[elem->id()];
            triangles.push_back(tri3);

            Tri3* tri4 = new Tri3;
            tri4->set_node(0) = centroids[elem->id()];
            tri4->set_node(1) = centroids[elem->neighbor(1)->id()];
            tri4->set_node(2) = elem->get_node(2);
            triangles.push_back(tri4);
          }

        else
          {
            libmesh_error();
          }
      }
  }

  // Loop again, delete the Quads
  {
    MeshBase::element_iterator       el     = mesh.elements_begin();
    const MeshBase::element_iterator end_el = mesh.elements_end();

    for (; el != end_el; ++el)
      {
        Elem* elem = *el;

        libmesh_assert (elem->type() == QUAD4);

        mesh.delete_elem(elem);
      }
  }
  // Now add the triangles
  {
    std::vector<Elem*>::iterator el        = triangles.begin();
    const std::vector<Elem*>::iterator end = triangles.end();
    for (; el != end; ++el)
      {
        mesh.add_elem(*el);
      }
  }


  // And add boundary IDs
  for (unsigned int e=0; e<saved_boundary_elements.size(); ++e)
    mesh.boundary_info->add_side(saved_boundary_elements[e],
                                 saved_bc_sides[e],
                                 saved_bc_ids[e]);

  // Prepare for use (find_neighbors, etc.)
  mesh.prepare_for_use();

  // If user requested quads, find pairs of triangles to replace with
  // a quad.
  if (type == QUAD4)
    {
      std::vector<Elem*> new_quads;
      std::set<unsigned int> flagged_tris;

      {
        MeshBase::element_iterator       el     = mesh.elements_begin();
        const MeshBase::element_iterator end_el = mesh.elements_end();

        for (; el != end_el; ++el)
          {
            Elem* elem = *el;

            // Skip elements on the boundary
            bool skip = false;
            for (unsigned int s=0; s<elem->n_sides(); s++)
              if (elem->neighbor(s) == NULL)
                {
                  skip = true;
                  break;
                }

            if (skip)
              continue;

            // Has this element already been flagged as part of a pair?
            if (flagged_tris.find(elem->id()) != flagged_tris.end())
              continue;

            // We could hit either the top or bottom element of the pair
            // first.  We already know the neighbor numbering of the triangle pairs
            // because of how we created the triangles.

            Elem* top_partner = NULL;
            Elem* bot_partner = NULL;

            // Try to find a top partner.
            if (elem->neighbor(1)->neighbor(0)->id() == elem->id())
              {
                top_partner = elem->neighbor(1);
                bot_partner = elem;
              }

            // Try to find a bottom partner
            else if (elem->neighbor(0)->neighbor(1)->id() == elem->id())
              {
                bot_partner = elem->neighbor(0);
                top_partner = elem;
              }

            // If a partner was found...
            if (top_partner && bot_partner)
              {
                //              std::cout << "Elements " << elem->id()
                //                        << " and " << partner->id()
                //                        << " form a pair!" << std::endl;

                // Insert these element IDs into the set
                flagged_tris.insert(top_partner->id());
                flagged_tris.insert(bot_partner->id());

                // Form a QUAD, save it in the new_quads vector.
                Quad4* quad4 = new Quad4;
                quad4->set_node(0) = bot_partner->get_node(0);
                quad4->set_node(1) = bot_partner->get_node(1);
                quad4->set_node(2) = top_partner->get_node(2);
                quad4->set_node(3) = top_partner->get_node(0);
                new_quads.push_back(quad4);
              }
          } // End first loop over elements
      }

      // Delete all the flagged tris
      {
        for (std::set<unsigned int>::iterator it=flagged_tris.begin();
             it != flagged_tris.end(); ++it)
          {
            mesh.delete_elem( mesh.elem(*it) );
          }
      }

      // And add in all the new Quads
      {
        std::vector<Elem*>::iterator el        = new_quads.begin();
        const std::vector<Elem*>::iterator end = new_quads.end();
        for (; el != end; ++el)
          {
            mesh.add_elem(*el);
          }
      }

      // Prepare for use (find_neighbors, etc.)
      mesh.prepare_for_use();
    }
}
