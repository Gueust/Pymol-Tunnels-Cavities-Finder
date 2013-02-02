// Copyright (c) 2009-2010  INRIA Sophia-Antipolis (France) (Initial part)
// Jean-Baptiste Lespiau (Ecole Polytechnique) and Julie Bernauer (Inria)
// All rights reserved.
#ifndef CGAL_CGO_WRITER_H
#define CGAL_CGO_WRITER_H

#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
#include <CGAL/Skin_surface_refinement_policy_3.h>

//#include <ESBTL/Fcc_lattice_generation.h>
#include <iostream>
#include <list>
#include <algorithm>
#include <iterator>

void print_normal(double x1, double y1, double z1, 
		  double x2, double y2, double z2,
		  double x3, double y3, double z3, std::ostream &out) {

  out << "NORMAL,"
      << (y2-y1)*(z3-z2) - (z2-z1)*(y3-y2) << ", "
      << (z2-z1)*(x3-x2) - (x2-x1)*(z3-z2) << ", "
      << (x2-x1)*(y3-y2) - (y2-y1)*(x3-x2) << "," 
      << std::endl;
}


template <class SkinSurface, class Polyhedron>
void write_cgo(SkinSurface &skin,
				   Polyhedron &p,
				   std::ostream &out)
{
  typedef typename Polyhedron::Vertex_iterator                  Vertex_iterator;
  typedef typename Polyhedron::Facet_iterator                   Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator HFC;
  typedef typename Polyhedron::Vertex_handle                    Vertex_handle;
  typedef typename Polyhedron::Traits::Vector_3                 Vector;

  CGAL::Skin_surface_refinement_policy_3<SkinSurface, Polyhedron> policy(skin);

  // Write header
  out << "from pymol import cmd" <<  std::endl
      << "from pymol.cgo import * # get constants" 
      << std::endl;

  out << "obj = [" << std::endl
      << "  BEGIN, TRIANGLES,"  << std::endl
      << "  COLOR, 0.8, 0.8, 0.8," << std::endl;

  /** It is not possible with Compiled Graphics Objects (CGOs)
    * which are Pymol abstraction of OpenGL objects to build
    * directly a Polygon. So we are forced to do it with triangles
    */

  for(Facet_iterator fi = p.facets_begin(); fi != p.facets_end(); ++fi) {
    // First halfedge
    HFC hc = fi->facet_begin();
    HFC hc_end = hc;
    CGAL_assertion( CGAL::circulator_size(j) >= 3);

    // We go to all the points of the polygon
    double x,y,z;
    bool first = true;
    do {
      Vertex_handle vh = (*hc).vertex();

      if (first) {
	Vector n = policy.normal(vh);
	out << "NORMAL,"
	    << n.x() << ", "
	    << n.y() << ", "
	    << n.z() << "," 
	    << std::endl;
	x = (*hc).vertex()->point().x();
	y = (*hc).vertex()->point().y();
	z = (*hc).vertex()->point().z();
	first = false;
	continue;
      }
      /*
      print_normal( (*hc_end).vertex()->point().x(),
		    (*hc_end).vertex()->point().x(),
		    (*hc_end).vertex()->point().x(),
		    x, y, z, 
		    (*hc).vertex()->point().x(),
		    (*hc).vertex()->point().y(),
		    (*hc).vertex()->point().z(), out); */
      // Source of the triangle
      out << "  VERTEX, " 
	  << (*hc_end).vertex()->point().x() << ", "
	  << (*hc_end).vertex()->point().y() << ", " 
	  << (*hc_end).vertex()->point().z() << ", "
	  << std::endl;
      // Previous point
      out << "  VERTEX, " 
	  << x  << ", "
	  << y  << ", "
	  << z  << ", "
	  <<  std::endl;
      // Current Point
      out << "  VERTEX, " 
	  << (*hc).vertex()->point().x() << ", "
	  << (*hc).vertex()->point().y() << ", "
	  << (*hc).vertex()->point().z() << ", " 
	  << std::endl <<  std::endl;
      // We update the previous point coordinates
      x = (*hc).vertex()->point().x();
      y = (*hc).vertex()->point().y();
      z = (*hc).vertex()->point().z();
    } while (++hc != hc_end);
    out << std::endl;
  }

  out << "  END" << std::endl
      << "]"   << std::endl
      << std::endl
      << "cmd.load_cgo(obj,'cgo01')" << std::endl;
}

#endif // CGAL_CGO_WRITER_H
