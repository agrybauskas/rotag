#ifndef SRC_LIB_COMBINATORICS_H_
#define SRC_LIB_COMBINATORICS_H_

#include <vector>
#include <iostream>

template <class T>
void permutation(size_t size,
                 std::vector<std::vector<T>> item_list,
                 std::vector<std::vector<T>> *permuted_list,
                 std::vector<T> base = {}) {
    const size_t base_size = base.size();
    if (base_size == size) {
        (*permuted_list).push_back(base);
    } else {
        const size_t last_idx = item_list.size() > 0 ? item_list.size() - 1 : 0;
        for (size_t i = 0; i < item_list[last_idx].size(); i++) {
            std::vector<T> updated_base = base;
            updated_base.push_back(item_list[last_idx][i]);
            permutation(size, item_list, permuted_list, updated_base);
        }
    }
}

#endif  // SRC_LIB_COMBINATORICS_H_
