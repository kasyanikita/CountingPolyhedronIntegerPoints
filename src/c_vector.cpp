#include "c_vector.h"

namespace GroupIP {

    Vector CVector::get_c_vector(int_t alpha) {
        bool flag = true;
        Vector c;
        std::vector<std::vector<std::vector<int_t>>> A_subs;
        for (int i = 0; i < A_.size(); ++i) {
            auto A_sub = get_sub_matrix(A_, i);
            A_subs.push_back(A_sub);
        }
        while (flag) {
            c = gen_rand_vector(A_subs[0].size(), -alpha, alpha);
            for (int i = 0; i < c.size(); ++i) {
            if (c[i] == 0) c[i] = 1;
            }
            flag = check_subs(A_subs, c);
        }
        return c;
    }

}  // namespace GroupIP