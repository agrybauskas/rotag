#ifndef _COMBINATORICS_H_
#define _COMBINATORICS_H_

#include <vector>

template <class T>
std::vector<std::vector<T>>
  permutation(int size,
              std::vector<std::vector<T>> base,
              std::vector<T> list,
              std::vector<std::vector<T>> permuted_list) {
  const int base_size = base.size();
  if(base_size == size) {
    /* permuted_list.push_back(base); */
  } else {
    for(size_t i = 0; i < list.size(); i++) {
      // const std::vector<std::vector<double>> updated_base =
      //   base.push_back(list[base_size][i]);
      // permutation(size, updated_base, list, permuted_list);
    }
  }
  return permuted_list;
}

#endif
