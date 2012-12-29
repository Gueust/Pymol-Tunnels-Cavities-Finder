// Adapted from ESBTL file skin_surface_pdb_reader.cpp using CGAL
//
// Copyright (c) 2009-2010  INRIA Sophia-Antipolis (France) (Initial part)
// Jean-Baptiste Lespiau (Ecole Polytechnique) and Julie Bernauer (Inria)
// All rights reserved.
//
//
// ESBTL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ESBTL is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ESBTL.  If not, see <http://www.gnu.org/licenses/>.
//
//
// Additional permission under GNU GPL version 3 section 7
//
// If you modify this Library, or any covered work, by linking or
// combining it with CGAL (or a modified version of that library), the
// licensors of this Library grant you additional permission to convey
// the resulting work. Corresponding Source for a non-source form of
// such a combination shall include the source code for the parts of CGAL
// used as well as that of the covered work.
//
//
// Initial author(s)     :  Nico Kruithof (?)
// Additional author : Jean-Baptiste Lespiau

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/Skin_surface_polyhedral_items_3.h>
#include <list>

#include <ESBTL/CGAL/EPIC_kernel_with_atom.h>
#include <ESBTL/default.h>

typedef ESBTL::CGAL::EPIC_kernel_with_atom                  K;
typedef ESBTL::CGAL::Default_system                         System;
typedef CGAL::Skin_surface_traits_3<K>                      Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::FT                                  FT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Weighted_point::Point                               Bare_point;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;

#include <list>
#include <fstream>
#include "skin_surface_writer.h"
#include "extract_balls_from_pdb.h"


int main(int argc, char *argv[]) {
  const char *filename;
  if (argc == 2) {
    filename = argv[1];
  } else {
    filename = "../data/1t7i.pdb";
  }
  
 
  std::list<Weighted_point> l;
  double                    shrinkfactor = 0.5;
  //Container for molecular system
  std::vector<System> systems;
  
  // Retrieve input balls:
  extract_balls_from_pdb<K>(filename,systems,std::back_inserter(l));
  
  // Construct skin surface:
  std::cout << "Constructing skin surface..." <<std::endl;
  Skin_surface_3 skin_surface(l.begin(), l.end(), shrinkfactor);

  // Extract mesh from the skin surface:
  std::cout << "Meshing skin surface..." <<std::endl;
  Polyhedron p;
  CGAL::mesh_skin_surface_3(skin_surface, p);

  std::ofstream out("mesh.off");
  out << p;

  return 0;
}