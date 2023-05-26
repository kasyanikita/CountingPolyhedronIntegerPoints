#ifndef _GLOBAL_DEFS_H_
#define _GLOBAL_DEFS_H_

#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <valarray>
#include <cassert>

#include "boost_tee_logging.h"
#include "timer.h"

namespace GroupIP
{
    using namespace std;

    using int_t = int64_t;
    using uint_t = uint64_t;

    void print_matrix(std::vector<std::vector<int_t>> &M, std::string name)
    {
        std::cout << name << ":\n";
        for (int i = 0; i < M.size(); ++i)
        {
            for (int j = 0; j < M[i].size(); ++j)
            {
                std::cout << M[i][j] << " ";
            }
            std::cout << '\n';
        }
        std::cout << "\n";
    }

    template <class T>
    void print_vector(std::vector<T> v, std::string name)
    {
        std::cout << name << ": ";
        for (int i = 0; i < v.size(); ++i)
        {
            std::cout << v[i] << " ";
        }
        std::cout << '\n';
    }
} // namespace GroupIP

#endif
