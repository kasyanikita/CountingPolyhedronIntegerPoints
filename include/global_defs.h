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

// #include "boost_tee_logging.h"
// #include "timer.h"

namespace GroupIP
{   
    int return66();
    using namespace std;

    using int_t = int64_t;
    using uint_t = uint64_t;

    using Matrix = std::vector<std::vector<int_t>>;
    using Vector = std::vector<int_t>;
    using uVector = std::vector<uint_t>;
} // namespace GroupIP

#endif
