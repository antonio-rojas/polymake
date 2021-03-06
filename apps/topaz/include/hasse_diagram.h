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

#ifndef POLYMAKE_TOPAZ_HASSE_DIAGRAM_H
#define POLYMAKE_TOPAZ_HASSE_DIAGRAM_H

#include "polymake/PowerSet.h"
#include "polymake/graph/BasicLatticeTypes.h"

namespace polymake { namespace topaz {

using graph::lattice::BasicDecoration;
using graph::lattice::FaceIndexingData;

// Iterate either over all maximal cells (for the artificial top node)
// or over all subsets less one element.
class simplicial_closure_iterator {
public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = Set<int>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = ptrdiff_t;

  template <typename Iterable>
  explicit simplicial_closure_iterator(const Iterable& f)
  {
    for (auto f_it = entire(f); !f_it.at_end(); ++f_it)
      data.push_back(*f_it);
    it = entire(data);
  }

  reference operator* () const { return it.operator*(); }
  pointer operator->() const { return it.operator->(); }
  simplicial_closure_iterator& operator++ () { ++it; return *this; }
  const simplicial_closure_iterator operator++ (int) { simplicial_closure_iterator copy = *this; operator++(); return copy; }
  bool at_end() const { return it.at_end(); }

protected:
  std::list<Set<int>> data;
  pm::iterator_range<std::list<Set<int>>::const_iterator> it;
};


template <typename Decoration = BasicDecoration>
class SimplicialClosure {
public:
  using ClosureData = Set<int>;

  SimplicialClosure(const IncidenceMatrix<>& facets_)
    : facets(facets_)
    , total_size(facets.cols()) {}

  ClosureData closure_of_empty_set() const
  {
    return sequence(0,facets.cols()+1);
  }

  ClosureData compute_closure_data(const Decoration& face) const
  {
    return face.face;
  }

  FaceIndexingData get_indexing_data(const ClosureData& data)
  {
    int &fi = face_index_map[data];
    return FaceIndexingData(fi, fi == -1, fi == -2);
  }

  simplicial_closure_iterator get_closure_iterator(const ClosureData& face) const
  {
    return face.size() > total_size
           ? simplicial_closure_iterator(rows(facets))
           : simplicial_closure_iterator(all_subsets_less_1(face));
  }

protected:
  const IncidenceMatrix<> facets;
  const int total_size;
  FaceMap<> face_index_map;
};

class SimplicialDecorator {
protected:
  const Set<int> artificial_set;
  int top_rank;
public:
  SimplicialDecorator(int top_rank_arg, const Set<int>& artificial_set_arg)
    : artificial_set(artificial_set_arg)
    , top_rank(top_rank_arg) {}

  template <typename TSet>
  BasicDecoration compute_initial_decoration(const GenericSet<TSet, int>& face) const
  {
    return BasicDecoration(artificial_set, top_rank);
  }

  template <typename TSet>
  BasicDecoration compute_decoration(const GenericSet<TSet, int>& face,
                                     const BasicDecoration& predecessor_data) const
  {
    return BasicDecoration(face, face.top().size());   
  }

  BasicDecoration compute_artificial_decoration(const NodeMap<Directed, BasicDecoration> &decor,
               const std::list<int>& max_nodes) const
  {
    return BasicDecoration(Set<int>(), 0);
  }
};

graph::Lattice<graph::lattice::BasicDecoration>
hasse_diagram_from_facets(const Array<Set<int>>& facets, const graph::lattice::RankRestriction& rr = graph::lattice::RankRestriction());

perl::Object upper_hasse_diagram(perl::Object complex, int boundary_rank);

} }

#endif
