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

#ifndef _POLYMAKE_SIMPLE_ROOTS_H
#define _POLYMAKE_SIMPLE_ROOTS_H

#include "polymake/SparseMatrix.h"
#include "polymake/Rational.h"
#include "polymake/QuadraticExtension.h"
#include "polymake/hash_map"
#include "polymake/hash_set"
#include "polymake/linalg.h"

namespace polymake { namespace polytope {

template<typename E>
perl::Object
root_system_impl(const SparseMatrix<E>& simple_roots)
{
   hash_set<SparseVector<E>> R_old, R_new;
   hash_map<SparseVector<E>,int> index_of;
   int index(0);
   for (; index < simple_roots.rows(); ++index) {
      index_of[simple_roots.row(index)] = index;
      R_new += simple_roots.row(index);
   }

   while (R_new != R_old) {
      R_old = R_new;
      for (const auto& r: R_old)
         for (const auto& s: R_old)
            R_new += reflect(s, r);
   }
   const int n = R_new.size();

   SparseMatrix<E> V(n, simple_roots.cols());
   for (const auto& r: Set<SparseVector<E>>(entire(R_new))) { // ensure canonical ordering of roots
      if (!index_of.exists(r))
         index_of[r] = index++;
      V.row(index_of[r]) = r;
   }
   
   Array<Array<int>> gens(simple_roots.rows());
   for (int i=0; i<gens.size(); ++i) {
      Array<int> gen(n);
      const SparseVector<E> h(simple_roots.row(i));
      for (int j=0; j<n; ++j)
         gen[j] = index_of[reflect(V.row(j), h)];
      gens[i] = gen;
   }

   perl::Object a("group::PermutationAction");
   a.take("GENERATORS") << gens;
   perl::Object g("group::Group");
   g.take("VECTOR_ACTION") << a;
   perl::Object R("VectorConfiguration", mlist<E>());
   R.take("VECTORS") << V;
   R.take("GROUP") << g;
   return R;
}
      
SparseMatrix<Rational> simple_roots_type_A (const int n);
SparseMatrix<Rational> simple_roots_type_B (const int n);
SparseMatrix<Rational> simple_roots_type_C (const int n);
SparseMatrix<Rational> simple_roots_type_D (const int n);
SparseMatrix< QuadraticExtension<Rational> > simple_roots_type_E6 ();
SparseMatrix< QuadraticExtension<Rational> > simple_roots_type_E7 ();
SparseMatrix<Rational> simple_roots_type_E8 ();
SparseMatrix<Rational> simple_roots_type_F4 ();
SparseMatrix<Rational> simple_roots_type_G2 ();
SparseMatrix<QuadraticExtension<Rational> > simple_roots_type_H3 ();
SparseMatrix<QuadraticExtension<Rational> > simple_roots_type_H4 ();

} }

#endif // _POLYMAKE_SIMPLE_ROOTS_H

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
