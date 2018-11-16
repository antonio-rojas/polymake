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

#ifndef POLYMAKE_POLYTOPE_GENERIC_CONVEX_HULL_CLIENT_H
#define POLYMAKE_POLYTOPE_GENERIC_CONVEX_HULL_CLIENT_H

#include "polymake/client.h"
#include "polymake/polytope/convex_hull.h"

namespace polymake { namespace polytope {

//! check column dimensions of given matrices, extend empty one if needed
//! @param A first matrix
//! @param B second matrix
//! @param prepend_0_column prepend a column of zero elements to non-empty matrices
//! @return success indicator
template <typename E>
bool align_matrix_column_dim(Matrix<E>& A, Matrix<E>& B, bool prepend_0_column)
{
   const int d = std::max(A.cols(), B.cols());
   for (auto m : { &A, &B }) {
      if (m->cols() != d) {
         if (m->rows() != 0 || m->cols() != 0)
            return false;
         m->resize(0, d);
      }
      if (prepend_0_column && d)
         *m = zero_vector<E>() | *m;
   }
   return true;
}

template <typename Scalar, typename Solver>
void generic_convex_hull_primal(perl::Object& p, bool isCone, const Solver& solver)
{
   Matrix<Scalar> Points = p.give("RAYS | INPUT_RAYS"),
               Lineality = p.lookup("LINEALITY_SPACE | INPUT_LINEALITY");

   if (!align_matrix_column_dim(Points, Lineality, isCone))
      throw std::runtime_error("convex_hull_primal - dimension mismatch between RAYS|INPUT_RAYS and LINEALITY_SPACE|INPUT_LINEALITY");

   const auto sol = solver.enumerate_facets(Points, Lineality, isCone);
   if (isCone) {
      p.take("FACETS") << sol.first.minor(All, range_from(1));
      p.take("LINEAR_SPAN") << sol.second.minor(All, range_from(1));
   } else {
      p.take("FACETS") << sol.first;
      p.take("AFFINE_HULL") << sol.second;
   }
}

template <typename Scalar, typename Solver>
void generic_convex_hull_dual(perl::Object& p, bool isCone, const Solver& solver)
{
   Matrix<Scalar> H = p.give("FACETS | INEQUALITIES"),
                 EQ = p.lookup("LINEAR_SPAN | EQUATIONS");

   if (!align_matrix_column_dim(H, EQ, isCone))
      throw std::runtime_error("convex_hull_dual - dimension mismatch between FACETS|INEQUALITIES and LINEAR_SPAN|EQUATIONS");

   // * we handle the case of polytopes with empty exterior description somewhat special:
   //    empty facet matrix implies that the polytope must be empty!
   //    empty inequalities cannot occur due to the far face initial rule
   // * this also covers the case when the ambient dimension is empty for polytopes
   // * for cones an empty exterior description describes the whole space and is handled
   //   correctly in the interfaces
   if (isCone || H.rows() > 0 || EQ.rows() > 0) {
      try {
         const auto sol = solver.enumerate_vertices(H, EQ, isCone);
         if (isCone) {
            p.take("RAYS") << sol.first.minor(All, range_from(1));
            p.take("LINEALITY_SPACE") << sol.second.minor(All, range_from(1));
         } else {
            p.take("RAYS") << sol.first;
            p.take("LINEALITY_SPACE") << sol.second;
         }
         p.take("POINTED") << (sol.second.rows()==0);
         p.take("LINEALITY_DIM") << sol.second.rows();
         return;
      }
      catch (const infeasible&) { }
   }
   const int d = H.cols();
   p.take("RAYS") << Matrix<Scalar>(0, d);
   p.take("LINEALITY_SPACE") << Matrix<Scalar>(0, d);
   p.take("LINEALITY_DIM") << 0;
   p.take("POINTED") << true;
}

} }

#endif // POLYMAKE_POLYTOPE_GENERIC_CONVEX_HULL_CLIENT_H

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
