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
#include "polymake/polytope/ppl_interface.h"
#include "polymake/polytope/generic_lp_client.h"

namespace polymake { namespace polytope { namespace ppl_interface {

template <typename Scalar=Rational>
auto create_LP_solver()
{
   return cached_LP_solver<Rational>(new LP_Solver<Rational>(), true);
}

}

void ppl_lp_client(perl::Object p, perl::Object lp, bool maximize)
{
   generic_lp_client<Rational>(p, lp, maximize, ppl_interface::LP_Solver<Rational>());
}

Function4perl(&ppl_lp_client, "ppl_lp_client(Polytope<Rational>, LinearProgram<Rational>, $)");

InsertEmbeddedRule("function ppl.simplex: create_LP_solver<Scalar> [Scalar==Rational] () : c++ (name => 'ppl_interface::create_LP_solver') : returns(cached);\n");

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
