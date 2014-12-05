/* Copyright (c) 1997-2014
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
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/matroid/revlex_bases.h"
#include <algorithm>

namespace polymake { namespace matroid {


std::string bases_to_revlex_encoding(const Array<Set<int> >& bases,
                                     int rank,
                                     int n_elements)
{
   return bases_to_revlex_encoding_impl(bases, rank, n_elements);
}

Array<Set<int> > bases_from_revlex_encoding(const std::string& revlex_encoding,
                                            int rank,
                                            int n_elements,
                                            perl::OptionSet options)
{
   const bool check_basis_exchange_axiom = options["check_basis_exchange_axiom"];
   return bases_from_revlex_encoding_impl(revlex_encoding, rank, n_elements, check_basis_exchange_axiom);
}

UserFunction4perl("# Encode the bases of a given matroid as a string."
                  "# All possible binom(n,r) tuples of indices are traversed in revlex order."
                  "# If the current tuple is a basis, a '*' is recorded, else a '0'."
                  "# @param Array<Set> bases the list of bases of the matroid"
                  "# @param Int r the rank of the matroid"
                  "# @param Int n the number of elements of the matroid"
                  "# @return String",
                  &bases_to_revlex_encoding,
                  "bases_to_revlex_encoding(Array<Set> $$)");

UserFunction4perl("# Decode the bases of a given matroid from a string containing"
                  "# all possible binom(n,r) tuples of indices in revlex order."
                  "# If the current tuple is a basis, a '*' is recorded, else a '0'."
                  "# @param String encoding the revlex encoding of the list of bases of the matroid"
                  "# @param Int r the rank of the matroid"
                  "# @param Int n the number of elements of the matroid"
                  "# @option Bool check_basis_exchange_axiom whether to perform the check of the axiom after construction"
                  "# @return Array<Set>",
                  &bases_from_revlex_encoding,
                  "bases_from_revlex_encoding(String $$ { check_basis_exchange_axiom => 0 })");


} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End: