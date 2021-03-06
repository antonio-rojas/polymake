/* Copyright (c) 1997-2018
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
   http://www.polymake.org

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

#include "polymake/client.h"
#include "polymake/graph/lattice_builder.h"

namespace polymake { namespace matroid {

using namespace graph;
using namespace graph::lattice;

template <typename IMatrix>
perl::Object lattice_of_flats(const GenericIncidenceMatrix<IMatrix>& mat_hyperplanes, int total_rank)
{
   const bool is_dual = mat_hyperplanes.rows() < mat_hyperplanes.cols();
   if (is_dual) {
      const int total = mat_hyperplanes.rows();
      BasicClosureOperator<> cop = BasicClosureOperator<>(total, T(mat_hyperplanes));
      TrivialCut<BasicDecoration> cut;
      BasicDecorator<> dec = BasicDecorator<>(mat_hyperplanes.cols(), total_rank, Set<int>());
      Lattice<BasicDecoration, Sequential> init_lattice;
      return static_cast<perl::Object>(lattice_builder::compute_lattice_from_closure<BasicDecoration>(cop, cut, dec,0, lattice_builder::Dual(), init_lattice));
   } else {
      const int total = mat_hyperplanes.cols();
      BasicClosureOperator<> cop = BasicClosureOperator<>(total, mat_hyperplanes);
      TrivialCut<BasicDecoration> cut;
      BasicDecorator<> dec = BasicDecorator<>(0, Set<int>());
      Lattice<BasicDecoration, Sequential> init_lattice;
      return static_cast<perl::Object>(lattice_builder::compute_lattice_from_closure<BasicDecoration>(cop, cut, dec,0, lattice_builder::Primal(), init_lattice));
   }
}

FunctionTemplate4perl("lattice_of_flats(IncidenceMatrix, $)");

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
