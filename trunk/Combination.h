#ifndef COMBINATION_H
#define COMBINATION_H

#include <vector>

namespace Combination {
  inline bool add(std::vector<int> & cnt, int elem, int vmax);
  inline bool add(std::vector<int> & cnt, int elem, const std::vector<int> & vmax);
  inline bool next_vector(std::vector<int> & vec, int vmax);
  inline bool next_vector(std::vector<int> & vec, const std::vector<int> & vmax);
  //
  inline int  makeIndexVector( int imax, std::vector< std::vector<int> > & nvec);
  inline int  makeIndexVector( const std::vector<int> & imax, std::vector< std::vector<int> > & nvec);
  inline int  getIndexVector(std::vector<int> & nvec, std::vector< std::vector<int> > & allVectors);

}; // end of Utils

#include "Combination.icc"

#endif
