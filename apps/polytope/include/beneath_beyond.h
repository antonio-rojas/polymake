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

#ifndef POLYMAKE_POLYTOPE_BENEATH_BEYOND_H
#define POLYMAKE_POLYTOPE_BENEATH_BEYOND_H

#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/polytope/convex_hull.h"

namespace polymake { namespace polytope {

template <typename Scalar>
class BeneathBeyondConvexHullSolver : public ConvexHullSolver<Scalar, CanEliminateRedundancies::yes> {
public:
   ~BeneathBeyondConvexHullSolver() override;

   convex_hull_result<Scalar>
   enumerate_facets(const Matrix<Scalar>& Points, const Matrix<Scalar>& Lineality, const bool isCone) const override;

   convex_hull_result<Scalar>
   enumerate_vertices(const Matrix<Scalar>& Inequalities, const Matrix<Scalar>& Equations, const bool isCone) const override;

   std::pair<Bitset, Set<int>> get_non_redundant_points(const Matrix<Scalar>& points, const Matrix<Scalar>& linealities, bool isCone) const override;

   std::pair<Bitset, Set<int>> get_non_redundant_inequalities(const Matrix<Scalar>& inequalities, const Matrix<Scalar>& equations, bool isCone) const override;

   Array<Set<int>>
   placing_triangulation(const Matrix<Scalar>& Points, const Matrix<Scalar>& Lineality = Matrix<Scalar>()) const;
};

} } // end namespace polymake


#endif // POLYMAKE_POLYTOPE_BENEATH_BEYOND_H

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
