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

#ifndef POLYMAKE_INTERNAL_SPARSE_LINALG_H
#define POLYMAKE_INTERNAL_SPARSE_LINALG_H

#include "polymake/vector"

namespace pm {

template <typename E>
std::enable_if_t<is_field<E>::value, E>
det(SparseMatrix<E> M)
{
   const int dim = M.rows();
   if (!dim) return one_value<E>();

   std::vector<int> column_permutation(dim);
   copy_range(entire(sequence(0, dim)), column_permutation.begin());
   E result = one_value<E>();

   for (auto pivotrow=entire(rows(M)); !pivotrow.at_end(); ++pivotrow) {
      if (pivotrow->empty()) return zero_value<E>();

      auto pivot=pivotrow->begin();
      const int pr=pivotrow.index(), pc=pivot.index();   // row and column index
      result *= *pivot;
      if (column_permutation[pr] != pc) {
         std::swap(column_permutation[pr], column_permutation[pc]);
         negate(result);
      }

      auto beneath=cross_direction(pivot);
      ++beneath;
      while (!beneath.at_end()) {
         // delete all elements below pivot
         int r=beneath.index();
         const E factor=(*beneath)/(*pivot);
         ++beneath;
         M[r] -= factor * M[pr];
      }
   }
   return result;
}

template <typename E>
std::enable_if_t<is_field<E>::value, SparseVector<E>>
reduce(SparseMatrix<E> M, SparseVector<E> V)
{
   const int n_cols=M.cols();
   int col=0;
   for (auto pivotrow=entire(rows(M)); !pivotrow.at_end() && col < n_cols; ++pivotrow) {
      if (pivotrow->empty()) continue;

      auto pivot=pivotrow->begin();
      const E pivotelem=*pivot;

      (*pivotrow) /= pivotelem;

      auto in_col = cross_direction(pivotrow->begin());
      for (++in_col; !in_col.at_end(); ) {
         const E factor=*in_col;
         const int r2=in_col.index();
         ++in_col;
         M.row(r2) -= (*pivotrow) * factor;
      }
      const E factor = V[pivot.index()];
      V -= (*pivotrow) * factor;
      ++col;
   }
   return V;
}

template <typename E>
std::enable_if_t<is_field<E>::value, SparseMatrix<E>>
inv(SparseMatrix<E> M)
{
   const int dim=M.rows();
   SparseMatrix<E> L=unit_matrix<E>(dim), R=unit_matrix<E>(dim);

   for (auto c=entire(cols(M)); !c.at_end(); ++c) {
      if (c->empty()) throw degenerate_matrix();

      auto in_col=c->begin();
      auto in_row=cross_direction(in_col);
      int pr=in_col.index(), pc=c.index();
      const E pivotelem=*in_col;
      M.row(pr) /= pivotelem;  L.row(pr) /= pivotelem;  ++in_col;
      while (! in_col.at_end()) {
         const E factor=*in_col;
         int r=in_col.index();  ++in_col;
         M.row(r) -= factor * M.row(pr);  L.row(r) -= factor * L.row(pr);
      }
      ++in_row;
      while (! in_row.at_end()) {
         R.col(in_row.index()) -= (*in_row) * R.col(pc);
         M.row(pr).erase(in_row++);
      }
   }
   R.permute_cols(attach_operation(rows(M), BuildUnary<operations::front_index>()));
   return R*L;
}

template <typename E, bool ensure_nondegenerate=true>
std::enable_if_t<is_field<E>::value, Vector<E>>
lin_solve(SparseMatrix<E> A, Vector<E> B)
{
   const int m=A.rows(), n=A.cols();
   int non_empty_rows=m-n;
   if (ensure_nondegenerate && non_empty_rows<0) throw underdetermined();

   for (auto r=entire(rows(A)); !r.at_end(); ++r) {
      const int pr=r.index();
      if (r->empty()) {
         if (ensure_nondegenerate && --non_empty_rows<0) throw degenerate_matrix();
         if (!is_zero(B[pr])) throw infeasible();
         continue;
      }
      auto in_row=r->begin();
      auto in_col=cross_direction(in_row);
      const E pivotelem=*in_row;
      if (!is_one(pivotelem)) {
         (*r) /= pivotelem;
         B[pr] /= pivotelem;
      }

      for (++in_col; !in_col.at_end(); ) {
         const E factor=*in_col;
         const int r2=in_col.index();
         ++in_col;
         A.row(r2) -= (*r) * factor;
         B[r2] -= B[pr] * factor;
      }
   }

   Vector<E> result(A.cols());
   for (auto r=entire<reversed>(rows(A)); !r.at_end(); ++r) {
      if (r->empty()) continue;
      auto in_row=r->begin();
      auto in_col=cross_direction(in_row);
      const E& elem=result[in_row.index()]=B[r.index()];
      while (!(--in_col).at_end())
         B[in_col.index()] -= elem * (*in_col);
   }
   return result;
}

} // end namespace pm

#endif // POLYMAKE_INTERNAL_SPARSE_LINALG_H

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
