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
#include "polymake/vector"
#include "polymake/list"
#include "polymake/linalg.h"
#include "polymake/polytope/lrs_interface.h"
#include "polymake/Map.h"
#include "polymake/Graph.h"
#include "polymake/Set.h"
#include "polymake/Series.h"
#include "polymake/IncidenceMatrix.h"

namespace polymake { namespace polytope {
namespace {

template <typename E, typename Matrix, typename Vector1, typename Vector2> inline
void assign_facet_through_points(const GenericMatrix<Matrix,E>& M,
                                 const GenericVector<Vector1,E>& V_cut,
                                 GenericVector<Vector2,E>& f)
{
   f=null_space(M)[0];
   if (f*V_cut > 0) f.negate();
}
}

template <typename SetTop>
perl::Object truncation(perl::Object p_in, const GenericSet<SetTop>& trunc_vertices, perl::OptionSet options)
{
   if (options.exists("cutoff") && options.exists("noc")) 
      throw std::runtime_error("truncation: cannot specify cutoff and noc options simultaneously");
   
   const bool pointed = p_in.give("POINTED");
   if (!pointed)
      throw std::runtime_error("truncation: input should be pointed");
      
   const bool noc = options["noc"],
      relabel = options["relabel"];

   const IncidenceMatrix<> VIF=p_in.give("VERTICES_IN_FACETS");
   const Graph<> G=p_in.give("GRAPH.ADJACENCY");

   bool inequalities;

   const int n_vertices=VIF.cols(), n_facets=VIF.rows();
   typedef Map<int,int> vertex_map_type;
   vertex_map_type vertex_map;          // truncated vertex => the first of the new vertices

   int n_vertices_out, n_trunc_vertices;
   if (trunc_vertices.top().empty()) 
      throw std::runtime_error("truncation: no vertices to truncate specified");
   if (trunc_vertices.top().front() < 0 || trunc_vertices.top().back() >= n_vertices)
      throw std::runtime_error("vertex numbers out of range");

   perl::Object p_out("Polytope<Rational>");
   if (pm::identical<SetTop, Set<int> >::value)
      p_out.set_description() << p_in.name() << " with vertices " << trunc_vertices << " truncated" << endl;

   n_trunc_vertices=trunc_vertices.top().size();
   n_vertices_out=n_vertices-n_trunc_vertices;

   for (typename Entire<SetTop>::const_iterator v=entire(trunc_vertices.top()); !v.at_end(); ++v) {
      vertex_map[*v]=n_vertices_out;
      n_vertices_out+=G.degree(*v);
   }

   int n_facets_out=n_facets+n_trunc_vertices;
   IncidenceMatrix<> VIF_out(n_facets_out, n_vertices_out);

   // first inherit the original facets along with untouched vertices in them
   if (n_trunc_vertices < n_vertices)
      copy(entire(rows(VIF.minor(All,~keys(vertex_map)))), rows(VIF_out).begin());

   int new_facet=n_facets;
   for (vertex_map_type::iterator tv=vertex_map.begin();  !tv.at_end();  ++tv, ++new_facet) {
      int new_vertex=tv->second;
      for (Entire< AdjacencyMatrix< Graph<> >::row_type >::const_iterator
              nb=entire(G.adjacent_nodes(tv->first));  !nb.at_end();  ++nb, ++new_vertex) {
         // the new vertex inherits the ridge from the truncated vertex,
         // and it belongs to the new facet
         (VIF_out.col(new_vertex) = VIF.col(tv->first) * VIF.col(*nb)) += new_facet;
      }
   }

   std::vector<std::string> labels_out;
   if (relabel) {
      std::vector<std::string> labels(n_vertices);
      read_labels(p_in, "VERTEX_LABELS", labels);
      labels_out.resize(n_vertices_out);
      copy(entire(select(labels, ~keys(vertex_map))), labels_out.begin());

      for (vertex_map_type::iterator tv=vertex_map.begin();  !tv.at_end();  ++tv) {
         int new_vertex=tv->second;
         for (Entire< AdjacencyMatrix< Graph<> >::row_type >::const_iterator
                 nb=entire(G.adjacent_nodes(tv->first));  !nb.at_end();  ++nb, ++new_vertex) {
            labels_out[new_vertex]=labels[tv->first] + '-' + labels[*nb];
         }
      }
   }

   if (noc) {
      if (p_in.exists("COMBINATORIAL_DIM")) {
         const int dim=p_in.give("COMBINATORIAL_DIM");
         p_out.take("COMBINATORIAL_DIM") << dim;
         inequalities = (dim == 1);
      } else {
         inequalities = true;
      }
   } else {
      Rational cutoff_factor(1,2);
      if (options["cutoff"] >> cutoff_factor && (cutoff_factor<=0 || cutoff_factor>1))
         throw std::runtime_error("cutoff factor must be within (0,1]");

      std::vector<int> renumber_vertices;
      if (cutoff_factor==1) {
         renumber_vertices.resize(n_vertices);
         copy(sequence(0).begin(), select(renumber_vertices, ~keys(vertex_map)).begin());
      }
      const int dim = p_in.give("CONE_DIM");
      inequalities = (cutoff_factor == 1 || dim == 2);
      
      const Matrix<Rational> V=p_in.give("VERTICES"),
         F=p_in.give("FACETS"),
         AH=p_in.give("AFFINE_HULL");

      Matrix<Rational> F_out=F / zero_matrix<Rational>(n_trunc_vertices, F.cols());
      lrs_interface::solver S;
      Matrix<Rational> orth(AH);
      if (orth.cols()) orth.col(0).fill(0);

      Rows< Matrix<Rational> >::iterator new_facet = rows(F_out).begin()+n_facets;
      for (vertex_map_type::iterator tv=vertex_map.begin();  !tv.at_end();  ++tv, ++new_facet) {
         const int v_cut_off=tv->first;
         Matrix<Rational> basis(G.out_degree(v_cut_off), V.cols());
         const bool simple_vertex=basis.rows()+AH.rows()==V.cols()-1;

         Rows< Matrix<Rational> >::iterator b=rows(basis).begin();
         for (Entire< Graph<>::adjacent_node_list >::const_iterator nb_v=entire(G.adjacent_nodes(v_cut_off)); !nb_v.at_end(); ++nb_v, ++b) {
            if (vertex_map.exists(*nb_v))
               *b = (1-cutoff_factor/2) * V[v_cut_off] + cutoff_factor/2 * V[*nb_v];
            else
               *b = (1-cutoff_factor) * V[v_cut_off] + cutoff_factor * V[*nb_v];
         }
         if (simple_vertex) {
            // calculate a hyperplane thru the basis points
            assign_facet_through_points(basis/orth, V[v_cut_off], new_facet->top());
         } else {
            // look for a valid separating hyperplane furthest from the vertex being cut off
            *new_facet=S.solve_lp(basis, orth, V[v_cut_off], false).second;
         }

         if (cutoff_factor==1) {
            // we must take care of coinciding vertices
            int new_vertex=tv->second;
            b=rows(basis).begin();
            for (Entire< Graph<>::adjacent_node_list >::const_iterator nb_v=entire(G.adjacent_nodes(v_cut_off));
                 !nb_v.at_end();  ++nb_v, ++new_vertex) {

               if (!simple_vertex && !is_zero((*new_facet)*(*b))) continue;       // doesn't touch this vertex

               int other_vertex;
               vertex_map_type::iterator otv=vertex_map.find(*nb_v);
               if (!otv.at_end()) {
                  // pairs of coinciding new vertices should not be handled twice
                  if (v_cut_off < *nb_v) continue;
                  other_vertex=otv->second;
                  // number of the opposite new vertex can be found only by enumeration of adjacent_nodes...
                  for (Entire< Graph<>::adjacent_node_list >::const_iterator other_nb_v=G.adjacent_nodes(*nb_v).begin();
                       *other_nb_v != v_cut_off;  ++other_nb_v, ++other_vertex) ;
               } else {
                  other_vertex=renumber_vertices[*nb_v];
               }
               // the new vertex coincides with the neighbor and disappears
               VIF_out.col(other_vertex) += VIF_out.col(new_vertex).back(); // the new facet
               VIF_out.col(new_vertex).clear();
               if (relabel) labels_out[new_vertex].clear();
            }
         }
      }

      // if cutoff is 1 some old facets might not be facets anymore
      // if dim is 1 then vertices are facets. So we truncate facets
      if (inequalities){
        F_out /= unit_vector<Rational>(F_out.cols(), 0);
        p_out.take("INEQUALITIES")  << F_out;
        p_out.take("EQUATIONS") << AH;
      } else {
        p_out.take("FACETS")  << F_out;
        p_out.take("AFFINE_HULL") << AH;
      }

      if (cutoff_factor==1) {
         VIF_out.squeeze();
         if (relabel) {
            labels_out.resize(std::remove(labels_out.begin()+n_vertices-n_trunc_vertices,
                                          labels_out.end(), std::string())
                              - labels_out.begin());
         }
      }
   }

   p_out.take("N_VERTICES") << VIF_out.cols();
   
   // if cutoff is 1 some old facets might not be facets anymore
   if (!inequalities){
     p_out.take("VERTICES_IN_FACETS") << VIF_out;
   }

   if (relabel)
      p_out.take("VERTEX_LABELS") << labels_out;
   
   return p_out;
}

perl::Object truncation(perl::Object p_in, const pm::all_selector&, perl::OptionSet options)
{
   const int n_verts = p_in.give("N_VERTICES");
   perl::Object p_out=truncation(p_in,sequence(0,n_verts),options);
   p_out.set_description() << p_in.name() << " with all vertices truncated" << endl;
   return p_out;
}

perl::Object truncation(perl::Object p_in, int vertex, perl::OptionSet options)
{
   perl::Object p_out=truncation(p_in,scalar2set(vertex),options);
   p_out.set_description() << p_in.name() << " with vertex " << vertex << " truncated" << endl;
   return p_out;
}

perl::Object truncation(perl::Object p_in, const Array<int>& verts, perl::OptionSet options)
{
   Set<int> trunc_vertices;
   for (Entire< Array<int> >::const_iterator vi = entire(verts); !vi.at_end(); ++vi)
      trunc_vertices += *vi;

   if (verts.size() != trunc_vertices.size())
      throw std::runtime_error("truncation: repeating vertex numbers in the list");
   
   return truncation(p_in,trunc_vertices,options);
}

UserFunctionTemplate4perl("# @category Producing a polytope from polytopes"
                          "# "
                          "# Cut off one or more vertices of a polyhedron."
                          "# "
                          "# The exact location of the cutting hyperplane(s) can be controlled by the"
                          "# option //cutoff//, a rational number between 0 and 1."
                          "# When //cutoff//=0, the hyperplane would go through the chosen vertex, thus cutting off nothing."
                          "# When //cutoff//=1, the hyperplane touches the nearest neighbor vertex of a polyhedron."
                          "# "
                          "# Alternatively, the option //noc// (no coordinates) can be specified to produce a"
                          "# pure combinatorial description of the resulting polytope, which corresponds to"
                          "# the cutoff factor 1/2."
                          "# @param Polytope P"
                          "# @param Set<Int> trunc_vertices the vertex/vertices to be cut off;"
                          "#   A single vertex to be cut off is specified by its number."
                          "#   Several vertices can be passed in a Set or in an anonymous array of indices: [n1,n2,...]"
                          "#   Special keyword __All__ means that all vertices are to be cut off."
                          "# @option Rational cutoff controls the exact location of the cutting hyperplane(s);"
                          "#   rational number between 0 and 1; default value: 1/2"
                          "# @option Bool noc produces a pure combinatorial description (in contrast to //cutoff//)"
                          "# @option Bool relabel creates an additional section [[VERTEX_LABELS]];"
                          "#   New vertices get labels of the form 'LABEL1-LABEL2', where LABEL1 is the original label"
                          "#   of the truncated vertex, and LABEL2 is the original label of its neighbor."
                          "# @return Polytope"
                          "# @author Kerstin Fritzsche (initial version)",
                          "truncation(Polytope * {cutoff=>undef, noc=>undef, relabel=>undef})");
} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End: