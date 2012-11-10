// Adapted from ESBTL file triangulation.cpp using CGAL
//
// Copyright (c) 2009-2010  INRIA Sophia-Antipolis (France) (Initial part)
// Jean-Baptiste Lespiau (Ecole Polytechnique) and Julie Bernauer (Inria)
// All rights reserved.
//
//
//ESBTL is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//ESBTL is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with ESBTL.  If not, see <http://www.gnu.org/licenses/>.
//
//
//Additional permission under GNU GPL version 3 section 7
//
//If you modify this Library, or any covered work, by linking or
//combining it with CGAL (or a modified version of that library), the
//licensors of this Library grant you additional permission to convey
//the resulting work. Corresponding Source for a non-source form of
//such a combination shall include the source code for the parts of CGAL
//used as well as that of the covered work.
//
//
//
// Initial author(s)     :  Sébastien Loriot
// Additional author : Jean-Baptiste Lespiau


#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <list>
#include <vector>

#include <ESBTL/constants.h>
#include <ESBTL/molecular_system.h>
#include <ESBTL/PDB.h>
#include <ESBTL/line_selectors.h>
#include <ESBTL/builder.h>
#include <ESBTL/line_reader.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/weighted_atom_iterator.h>
#include <ESBTL/occupancy_handlers.h>
#include <ESBTL/CGAL/EPIC_kernel_with_atom.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Alpha_shape_3.h>

typedef ESBTL::CGAL::EPIC_kernel_with_atom Kernel;
typedef ESBTL::CGAL::Default_system My_system;
typedef CGAL::Delaunay_triangulation_3<Kernel>                Delaunay;


typedef ESBTL::Generic_classifier<ESBTL::Radius_of_atom<double,My_system::Atom> >  T_Atom_classifier;
typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> >             Accept_none_occupancy_policy;

//Delaunay
typedef CGAL::Delaunay_triangulation_3<Kernel>                       Delaunay;
//Regular
typedef CGAL::Regular_triangulation_euclidean_traits_3<Kernel>       Regular_traits;
typedef CGAL::Regular_triangulation_3<Regular_traits>                   Regular;
//Alpha_shape
typedef Regular_traits                                                  Alpha_gt;
typedef CGAL::Alpha_shape_vertex_base_3<Alpha_gt>                       Alpha_Vb;
typedef CGAL::Alpha_shape_cell_base_3<Alpha_gt>                         Alpha_Fb;
typedef CGAL::Triangulation_data_structure_3<Alpha_Vb,Alpha_Fb>         Alpha_Tds;
typedef CGAL::Regular_triangulation_3<Alpha_gt,Alpha_Tds>               Alpha_Triangulation_3;
typedef CGAL::Alpha_shape_3<Alpha_Triangulation_3>                      Alpha_shape_3;
typedef Alpha_shape_3::Alpha_iterator               Alpha_iterator;

//Iterator for regular triangulation
typedef ESBTL::Weighted_atom_iterator<My_system::Model,
                                      CGAL::Weighted_point<Kernel::Point_3,double>,
                                      ESBTL::Weight_of_atoms<T_Atom_classifier> > Weighted_atom_iterator;

//Skin surface
typedef CGAL::Skin_surface_traits_3<Kernel>                          Skin_traits;
typedef CGAL::Skin_surface_3<Skin_traits>                               Skin_surface_3;

int main(int argc, char** argv){

  if (argc != 2 ){
    std::cerr << "Please provide a filename"  << std::endl;
    return EXIT_FAILURE;
  }

  ESBTL::PDB_line_selector_two_systems sel;

  std::vector<My_system> systems;
  ESBTL::All_atom_system_builder<My_system> builder(systems,sel.max_nb_systems());



  T_Atom_classifier atom_classifier;

  if (ESBTL::read_a_pdb_file(argv[1],sel,builder,Accept_none_occupancy_policy())){
    if ( systems.empty() || systems[0].has_no_model() ){
      std::cerr << "No atoms found" << std::endl;
      return EXIT_FAILURE;
    }
    //Consider only the first model of the first system
    const My_system::Model& model=*systems[0].models_begin();
    unsigned nb_atm=0;
    unsigned nb_hetatm=0;

    for (My_system::Model::Atoms_const_iterator it_atm=model.atoms_begin();it_atm!=model.atoms_end();++it_atm){
      if (it_atm->is_hetatm())
        ++nb_hetatm;
      else
        ++nb_atm;
    }
    std::cout << "(nb_atm,nb_hetatm) = (" <<nb_atm<<","<<nb_hetatm<<")" << std::endl;

    /*
    Delaunay dt(model.atoms_begin(),model.atoms_end());
    std::cout << dt.finite_vertices_begin()->point().atom_name() << std::endl;
    std::cout << "Delaunay " << dt.number_of_vertices()<< std::endl;


    Regular rt(Weighted_atom_iterator(model.atoms_begin(),&atom_classifier),Weighted_atom_iterator(model.atoms_end(),&atom_classifier));
    std::cout << "Regular  " << rt.number_of_vertices()<< std::endl;
    std::cout << rt.finite_vertices_begin()->point().weight() << std::endl;
    std::cout << atom_classifier.get_properties(rt.finite_vertices_begin()->point().point()).value() << std::endl;
    std::cout << rt.finite_vertices_begin()->point().atom_name() << std::endl;
 */

    /*
     *   compute alpha shape
  Alpha_shape_3 as(lp.begin(),lp.end());
  std::cout << "Alpha shape computed in REGULARIZED mode by default"
	    << std::endl;

  // find optimal alpha value
  Alpha_iterator opt = as.find_optimal_alpha(1);
  std::cout << "Optimal alpha value to get one connected component is "
	    <<  *opt    << std::endl;
  as.set_alpha(*opt);
  assert(as.number_of_solid_components() == 1);
     */
    Alpha_shape_3 as(Weighted_atom_iterator(model.atoms_begin(),&atom_classifier),Weighted_atom_iterator(model.atoms_end(),&atom_classifier));
    std::cout << "Alpha    " << as.number_of_vertices()<< std::endl;
    std::cout << as.finite_vertices_begin()->point().atom_name() << std::endl;

    // find optimal alpha value
    Alpha_iterator opt = as.find_optimal_alpha(1);
    std::cout << "Optimal alpha value to get one connected component is " <<  *opt    << std::endl;
    as.set_alpha(*opt);

    std::cout << "Fin de l'alpha" << std::endl;
   // Skin_surface_3 skin(Weighted_atom_iterator(model.atoms_begin(),&atom_classifier),Weighted_atom_iterator(model.atoms_end(),&atom_classifier),0.5);
   // std::cout << "Skin    " << std::endl;

  }
  else
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}


/*
 *
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

#include <fstream>
#include <list>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Gt;

typedef CGAL::Alpha_shape_vertex_base_3<Gt>          Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>            Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb>  Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>         Alpha_shape_3;

typedef Gt::Point_3                                  Point;
typedef Alpha_shape_3::Alpha_iterator               Alpha_iterator;

int main()
{
  std::list<Point> lp;

  //read input
  std::ifstream is("../data/bunny_1000");
  int n;
  is >> n;
  std::cout << "Reading " << n << " points " << std::endl;
  Point p;
  for( ; n>0 ; n--)    {
    is >> p;
    lp.push_back(p);
  }

  // compute alpha shape
  Alpha_shape_3 as(lp.begin(),lp.end());
  std::cout << "Alpha shape computed in REGULARIZED mode by default"
	    << std::endl;

  // find optimal alpha value
  Alpha_iterator opt = as.find_optimal_alpha(1);
  std::cout << "Optimal alpha value to get one connected component is "
	    <<  *opt    << std::endl;
  as.set_alpha(*opt);
  assert(as.number_of_solid_components() == 1);
  return 0;
}7
*/
