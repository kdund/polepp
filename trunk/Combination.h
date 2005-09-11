#ifndef COMBINATION_H
#define COMBINATION_H

#include <vector>

namespace Combination {
  bool add(std::vector<int> & cnt, int elem, int vmax);
  bool add(std::vector<int> & cnt, int elem, const std::vector<int> & vmax);
  bool next_vector(std::vector<int> & vec, int vmax);
  bool next_vector(std::vector<int> & vec, const std::vector<int> & vmax);
}; // end of Utils

#endif
