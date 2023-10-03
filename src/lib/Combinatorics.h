#ifndef _COMBINATORICS_H_
#define _COMBINATORICS_H_

#include <vector>
#include <iostream>

template <class T>
std::vector<std::vector<T>>
  permutation(size_t size,
              std::vector<std::vector<T>> non_permuted_list,
              std::vector<std::vector<T>> permuted_list = {{}},
              std::vector<T> base = {}) {
  const size_t base_size = base.size();
  if(base_size == size) {
    permuted_list.push_back(base);
  } else {
    for(size_t i = 0; i < non_permuted_list[base_size].size(); i++) {
      std::vector<T> updated_base = base;
      updated_base.push_back(non_permuted_list[base_size][i]);
      permutation(size, non_permuted_list, permuted_list , updated_base);
    }
  }
  return permuted_list;
}

#endif
