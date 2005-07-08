#include "Combination.h"
namespace Combination {
  bool add(std::vector<int> & cnt, int elem, int max) {
    if (elem<0) return false;
    if (elem>int(cnt.size())-1) return false;
    cnt[elem]++;
    if (cnt[elem]>max) {
      cnt[elem] = 0;
      --elem;
      return add(cnt,elem,max);
    }
    return true;
  }

  bool next_vector(std::vector<int> & vec, int max) {
    int elem = int(vec.size())-1;
    return add(vec,elem,max);
  }
};
