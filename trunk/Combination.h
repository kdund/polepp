#ifndef COMBINATION_H
#define COMBINATION_H

#include <vector>

namespace Combination {
  bool add(std::vector<int> & cnt, int elem, int max);
  bool add(std::vector<int> & cnt, int elem, std::vector<int> max);
  bool next_vector(std::vector<int> & vec, int max);
  bool next_vector(std::vector<int> & vec, std::vector<int> max);
}; // end of Utils

#endif
