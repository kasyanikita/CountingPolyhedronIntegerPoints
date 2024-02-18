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
        std::vector<int_t> den;

        std::vector<std::vector<ExpPoly>> init_dp();
        std::vector<int_t> get_denominator();
        void normalize();

    public:
     Dynamic(const std::vector<std::vector<int_t>> &, const std::vector<int_t> &,
             const std::vector<int_t>&);
     void init();
     ExpPoly &d(uint_t k, GroupElement g);
     std::vector<std::vector<ExpPoly>> get_table();
     void init_first_layer();
     std::vector<int_t> get_den();
     ExpPoly get_final_poly();
    //  void start();
     void new_start();
     ExpPoly operator()(uint_t k, const GroupElement &ge);
    };
}

#endif  // COUNTINGINTEGERPOINTS_DYNAMIC_H_