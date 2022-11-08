#pragma once

#include "global_defs.h"

namespace GroupIP
{
    class ExpPoly
    {
        std::unordered_map<double, double> poly;

    public:
        friend ExpPoly operator*(double c, ExpPoly &exp_poly);
        friend ExpPoly operator*(ExpPoly &exp_poly, double c);
        ExpPoly(const std::vector<double> &exp, const std::vector<double> &base)
        {
            assert(exp.size() == base.size());
            for (size_t i = 0; i < exp.size(); ++i)
            {
                auto it = poly.find(exp[i]);
                if (it != poly.end())
                {
                    poly[exp[i]] += base[i];
                }
                else
                {
                    poly.insert(std::make_pair(exp[i], base[i]));
                }
            }
        }
        ExpPoly operator+(ExpPoly &rhs) const
        {
            ExpPoly res(*this);
            for (auto [k, v] : rhs.poly) {
                auto it = res.poly.find(k);
                if (it != res.poly.end()) {
                    res.poly[k] += v;
                } else {
                    res.poly[k] = v;
                }
            }
            return res;
        }

        void print_poly()
        {
            for (auto [k, v] : poly)
            {
                if (v == 1)
                {
                    std::cout << "e^" << k << "x"
                              << " ";
                }
                else
                {
                    std::cout << v << "*e^" << k << "x"
                              << " ";
                }
            }
            std::cout << '\n';
        }
    };

    ExpPoly operator*(double c, ExpPoly &exp_poly)
    {
        ExpPoly res(exp_poly);
        for (auto [k, v] : res.poly)
        {
            res.poly[k] *= c;
        }
        return res;
    }

    ExpPoly operator*(ExpPoly &exp_poly, double c)
    {
        return c * exp_poly;
    }

} // namespace GroupIP