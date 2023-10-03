#ifndef _COMBINATORICS_H_
#define _COMBINATORICS_H_

#include <vector>

template <class T>
std::vector<std::vector<T>>
  permutation(int size,
              std::vector<std::vector<T>> list,
              std::vector<T> base = {},
              std::vector<std::vector<T>> permuted_list = {{}}) {
  const int base_size = base.size();
  if(base_size == size) {
    permuted_list.push_back(base);
  } else {
    for(size_t i = 0; i < list[base_size].size(); i++) {
      std::vector<std::vector<T>> updated_base = base;
      updated_base.push_back(list[base_size]);
      permutation(size, list, updated_base, permuted_list);
    }
  }
  return permuted_list;
}

#endif
