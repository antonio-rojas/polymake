/*
 * Normaliz
 * Copyright (C) 2007-2013  Winfried Bruns, Bogdan Ichim, Christof Soeger
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

//---------------------------------------------------------------------------

// The matrix Map is assumed to be such that the rank of Map equals
// the number of columns of Map and no test is performed in the constructor

//---------------------------------------------------------------------------
#ifndef SIMPLEX_H
#define SIMPLEX_H

//---------------------------------------------------------------------------

#include <vector>
#include <list>

#include "libnormaliz.h"
#include "cone.h"
#include "HilbertSeries.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using std::list;
using std::vector;

template<typename Integer> class Full_Cone;

template<typename Integer>
class Simplex {
    size_t dim;
    Integer volume;
    vector<key_t> key;
    Matrix<Integer> Generators;
    vector< Integer > diagonal;
    vector< Integer > multiplicators;
    Matrix<Integer> Support_Hyperplanes;

//---------------------------------------------------------------------------
public:
//                      Construction and destruction
    Simplex(const Matrix<Integer>& Map);  //contructor of the first in lexicographic
    //order simplex inside Map, the rank of Map is assumed to equal the number of
    //columns of Map
    Simplex(const vector<key_t>& k, const Matrix<Integer>& Map); //main constuctor
    //the rank of M is assumed to equal the number of columns of M

//                          Data acces
    size_t read_dimension() const;              // returns dim
    void write_volume(const Integer& vol);  // writes volume
    Integer read_volume() const;            // returns volume
    vector<key_t> read_key() const;          // returns key
    Matrix<Integer> read_generators() const;        // returns generators
    vector<Integer> read_diagonal() const;    // returns diagonal
    vector<Integer> read_multiplicators() const;    // returns multiplicators
    Matrix<Integer> read_support_hyperplanes() const;  // returns the support hyperplanes

};
//class Simplex end 


template<typename Integer>
class SimplexEvaluator {
    Full_Cone<Integer> * C_ptr;
    size_t dim;
    //Integer volume;
    Integer det_sum; // sum of the determinants of all evaluated simplices
    mpq_class mult_sum; // sum of the multiplicities of all evaluated simplices
    size_t collected_elements_size;
    Matrix<Integer> Generators;
    Matrix<Integer> TGenerators;
    Matrix<Integer> GenCopy;
    Matrix<Integer> InvGenSelRows; // selected rows of inverse of Gen
    Matrix<Integer> InvGenSelCols; // selected columns of inverse of Gen
    Matrix<Integer> Sol;
    Matrix<Integer> InvSol;
    vector< Integer > GDiag; // diagonal of generator matrix after trigonalization
    vector< Integer > TDiag; // diagonal of transpose of generaor matrix after trigonalization
    vector< bool > Excluded;
    vector< Integer > Indicator; 
    vector< long > gen_degrees;
    vector< num_t > hvector;  //h-vector of the current evaluation
    HilbertSeries Hilbert_Series; //this is the summed Hilbert Series
    list< vector<Integer> > Candidates;
    list< vector<Integer> > Hilbert_Basis;
    list< vector<Integer> > Collected_Elements;
    //temporary objects are kept to prevent repeated alloc and dealloc
    Matrix<Integer> RS; // right hand side to hold order vector
    // Matrix<Integer> RSmult; // for multiple right hand sides


    //checks whether a new element is reducible by the local Hilbert basis
    bool is_reducible_interior(const vector< Integer >& new_element);

    bool isDuplicate(const vector<Integer>& cand) const;

	void addMult(const Integer& volume);

//---------------------------------------------------------------------------

public:

    SimplexEvaluator(Full_Cone<Integer>& fc);

    // full evaluation of the simplex, writes data back to the cone,
    // returns volume
    Integer evaluate(SHORTSIMPLEX<Integer>& s);

    // moves the union of Hilbert basis / deg1 elements to the cone
    // for partial triangulation it merges the sorted list
    void transfer_candidates();
    // returns sum of the determinants of all evaluated simplices
    Integer getDetSum() const;
    // returns sum of the multiplicities of all evaluated simplices
    mpq_class getMultiplicitySum() const;
    // returns sum of the Hilbert Series of all evaluated simplices
    const HilbertSeries& getHilbertSeriesSum() const;
};
//class SimplexEvaluater end

}

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
