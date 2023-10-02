#include "Combinatorics.h"

#include <iostream>
#include <vector>

std::vector<std::vector<double>>
  permutation(int size,
              std::vector<std::vector<double>> base,
              std::vector<double> list,
              std::vector<std::vector<double>> permuted_list) {
  const int base_size = base.size();
  if(base_size == size) {
    permuted_list.push_back(base);
  } else {
    for(size_t i = 0; i < list.size(); i++) {
      // const std::vector<std::vector<double>> updated_base =
      //   base.push_back(list[base_size][i]);
      // permutation(size, updated_base, list, permuted_list);
    }
  }
  return permuted_list;
}
