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

#ifndef POLYMAKE_POLYTOPE_CONVEX_HULL_H
#define POLYMAKE_POLYTOPE_CONVEX_HULL_H

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/Bitset.h"
#include "polymake/internal/linalg_exceptions.h"

namespace polymake { namespace polytope {

// first = facets or vertices (cone: rays)
// second = affine hull (cone: linear span) or lineality space
template <typename Scalar>
using convex_hull_result = std::pair<Matrix<Scalar>, Matrix<Scalar>>;

// not every solver offers direct support of redundant point/inequality elimination
enum class CanEliminateRedundancies { no, yes };

template <typename Scalar, CanEliminateRedundancies can_eliminate = CanEliminateRedundancies::no>
class ConvexHullSolver {
public:
   virtual ~ConvexHullSolver() {}

   virtual convex_hull_result<Scalar> enumerate_facets(const Matrix<Scalar>& points, const Matrix<Scalar>& linealities, bool isCone) const = 0;

   virtual convex_hull_result<Scalar> enumerate_vertices(const Matrix<Scalar>& inequalities, const Matrix<Scalar>& equations, bool isCone) const = 0;
};

template <typename Scalar>
class ConvexHullSolver<Scalar, CanEliminateRedundancies::yes>
   : public ConvexHullSolver<Scalar, CanEliminateRedundancies::no> {
public:
   // first = non-redundant points (vertices) as row indexes into points
   // second = non-redundant linealities as row indexes into points (when i < points.rows()) or linealities (i - points.rows())
   virtual std::pair<Bitset, Set<int>> get_non_redundant_points(const Matrix<Scalar>& points, const Matrix<Scalar>& linealities, bool isCone) const = 0;

   // first = non-redundant inequalities (facets) as row indexes into inequalities
   // second = non-redundant equations as row indexes into inequalities (when i < inequalities.rows()) or equations (i - inequalities.rows())
   virtual std::pair<Bitset, Set<int>> get_non_redundant_inequalities(const Matrix<Scalar>& inequalities, const Matrix<Scalar>& equations, bool isCone) const = 0;
};

template <typename Scalar, CanEliminateRedundancies can_eliminate = CanEliminateRedundancies::no>
using cached_convex_hull_solver = CachedObjectPointer<ConvexHullSolver<Scalar, can_eliminate>, Scalar>;

template <typename Scalar, CanEliminateRedundancies can_eliminate = CanEliminateRedundancies::no>
const ConvexHullSolver<Scalar, can_eliminate>& get_convex_hull_solver()
{
   static cached_convex_hull_solver<Scalar, can_eliminate> solver_ptr("polytope::create_convex_hull_solver");
   return can_eliminate == CanEliminateRedundancies::yes ? solver_ptr.get(can_eliminate) : solver_ptr.get();
}

// convenience wrappers
template <typename Scalar, typename Matrix1, typename Matrix2>
convex_hull_result<Scalar> enumerate_vertices(const GenericMatrix<Matrix1, Scalar>& inequalities, const GenericMatrix<Matrix2, Scalar>& equations,
                                              bool isCone)
{
   const auto& solver = get_convex_hull_solver<Scalar>();
   return solver.enumerate_vertices(convert_to_persistent_dense(inequalities.top()), convert_to_persistent_dense(equations.top()), isCone);
}

template <typename Scalar, typename Matrix1, typename Matrix2>
convex_hull_result<Scalar> enumerate_facets(const GenericMatrix<Matrix1, Scalar>& points, const GenericMatrix<Matrix2, Scalar>& linealities,
                                            bool isCone)
{
   const auto& solver = get_convex_hull_solver<Scalar>();
   return solver.enumerate_facets(convert_to_persistent_dense(points.top()), convert_to_persistent_dense(linealities.top()), isCone);
}

// convenience wrappers for empty equations/lineality
template <typename Scalar, typename Matrix1>
convex_hull_result<Scalar> enumerate_vertices(const GenericMatrix<Matrix1, Scalar>& inequalities, bool isCone)
{
   return enumerate_vertices(inequalities, Matrix<Scalar>(0, inequalities.cols()), isCone);
}

template <typename Scalar, typename Matrix1>
convex_hull_result<Scalar> enumerate_facets(const GenericMatrix<Matrix1, Scalar>& points, bool isCone)
{
   return enumerate_vertices(points.top(), Matrix<Scalar>(0, points.cols()), isCone);
}

// try to compute vertices, return empty matrices in the infeasible case
template <typename Scalar, typename Matrix1, typename Matrix2>
convex_hull_result<Scalar> try_enumerate_vertices(const GenericMatrix<Matrix1, Scalar>& inequalities, const GenericMatrix<Matrix2, Scalar>& equations,
                                                  bool isCone)
{
   try {
      return enumerate_vertices(inequalities, equations, isCone);
   }
   catch (const infeasible&) {
      const int d = std::max(inequalities.cols(), equations.cols());
      return { Matrix<Scalar>(0, d), Matrix<Scalar>(0, d) };
   }
}

template <typename Scalar, typename Matrix1, typename Matrix2>
std::pair<Bitset, Set<int>> get_non_redundant_points(const GenericMatrix<Matrix1, Scalar>& points, const GenericMatrix<Matrix2, Scalar>& linealities,
                                                     bool isCone)
{
   const auto& solver = get_convex_hull_solver<Scalar, CanEliminateRedundancies::yes>();
   return solver.get_non_redundant_points(convert_to_persistent_dense(points.top()), convert_to_persistent_dense(linealities.top()), isCone);
}

template <typename Scalar, typename Matrix1, typename Matrix2>
std::pair<Bitset, Set<int>> get_non_redundant_inequalities(const GenericMatrix<Matrix1, Scalar>& inequalities, const GenericMatrix<Matrix2, Scalar>& equations,
                                                           bool isCone)
{
   const auto& solver = get_convex_hull_solver<Scalar, CanEliminateRedundancies::yes>();
   return solver.get_non_redundant_inequalities(convert_to_persistent_dense(inequalities.top()), convert_to_persistent_dense(equations.top()), isCone);
}

} }

#endif // POLYMAKE_POLYTOPE_CONVEX_HULL_H

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
