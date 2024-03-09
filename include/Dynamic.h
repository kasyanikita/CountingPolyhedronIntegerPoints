#ifndef COUNTINGINTEGERPOINTS_DYNAMIC_H_
#define COUNTINGINTEGERPOINTS_DYNAMIC_H_

#include <stdexcept>
#include "global_defs.h"
#include "group_element.h"
#include "ExpPoly.h"
#include "snf_class.h"
#include "auxiliary_dynamic.h"
#include "Tools.h"
#include <flint/fmpz_mat.h>

namespace GroupIP
{
    class Dynamic
    {
        SNFClass snf_;
        AuxiliaryDynamic aux_;
        int_t n;
        std::vector<int_t> c_;
        std::vector<int_t> b_;
        std::vector<std::vector<ExpPoly>> dp;
        std::vector<std::vector<bool>> isComputed;
        std::vector<std::vector<ExpPoly>> init_dp();

    public:
     Dynamic(const std::vector<std::vector<int_t>> &, const std::vector<int_t> &,
             const std::vector<int_t>&);
     void init();
     ExpPoly &d(uint_t k, GroupElement g);
     ExpPoly operator()(uint_t k, const GroupElement &ge);
    };
}

#endif  // COUNTINGINTEGERPOINTS_DYNAMIC_H_