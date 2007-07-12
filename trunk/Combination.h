#ifndef COMBINATION_H
#define COMBINATION_H

#include <vector>

namespace Combination {
  template<typename T> inline bool add(std::vector<T> & cnt, int elem, T vmax);
  template<typename T> inline bool add(std::vector<T> & cnt, int elem, const std::vector<T> & vmax);
  template<typename T> inline bool next_vector(std::vector<T> & vec, T vmax);
  template<typename T> inline bool next_vector(std::vector<T> & vec, const std::vector<T> & vmax);
  //
  template<typename T> inline size_t  makeIndexVector( T imax, std::vector< std::vector<T> > & nvec);
  template<typename T> inline size_t  makeIndexVector( const std::vector<T> & imax, std::vector< std::vector<T> > & nvec);
  template<typename T> inline int     getIndexVector(std::vector<T> & nvec, std::vector< std::vector<T> > & allVectors);

}; // end of Utils

#include "Combination.icc"

#endif
