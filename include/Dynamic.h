#ifndef COUNTINGINTEGERPOINTS_DYNAMIC_H_
#define COUNTINGINTEGERPOINTS_DYNAMIC_H_

#include <stdexcept>
#include "global_defs.h"
#include "group_element.h"
#include "ExpPoly.h"
#include "snf_class.h"
#include "Tools.h"
#include <flint/fmpz_mat.h>

namespace GroupIP
{
    class Dynamic
    {
        int_t n_;
        Vector c_;
        Matrix h_;
        std::vector<GroupIP::GroupElement> g_;
        std::vector<std::vector<ExpPoly>> dp;
        std::vector<std::vector<bool>> isComputed;
        std::vector<std::vector<ExpPoly>> init_dp(const Vector&) const;

    public:
     Dynamic(const Vector &, const std::vector<GroupIP::GroupElement> &,
             const Matrix &);
     void Init(const Vector &);
     ExpPoly &d(uint_t k, GroupElement g);
     ExpPoly operator()(uint_t k, const GroupElement &ge);
    };
}

#endif  // COUNTINGINTEGERPOINTS_DYNAMIC_H_