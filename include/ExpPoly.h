#pragma once

#include "global_defs.h"

namespace GroupIP
{
    class ExpPoly
    {
    public:
        using exp_t = int_t;
        using coeff_t = uint_t;
    private:

        std::unordered_map<exp_t, coeff_t> poly;

    public:
        friend ExpPoly operator*(coeff_t c, ExpPoly &exp_poly);
        friend ExpPoly operator*(ExpPoly &exp_poly, coeff_t c);
        friend std::ostream &operator<<(std::ostream &out, ExpPoly &exp_poly);

        ExpPoly(const std::vector<exp_t> &exp, const std::vector<coeff_t> &coeff)
        {
            assert(exp.size() == coeff.size());
            for (size_t i = 0; i < exp.size(); ++i)
            {
                auto it = poly.find(exp[i]);
                if (it != poly.end())
                {
                    poly[exp[i]] += coeff[i];
                }
                else
                {
                    poly.insert(std::make_pair(exp[i], coeff[i]));
                }
            }
        }

        ExpPoly operator+(ExpPoly &rhs) const
        {
            ExpPoly res(*this);
            for (auto [k, v] : rhs.poly)
            {
                auto it = res.poly.find(k);
                if (it != res.poly.end())
                {
                    res.poly[k] += v;
                }
                else
                {
                    res.poly[k] = v;
                }
            }
            return res;
        }

        ExpPoly monomial_multiply(exp_t exp, coeff_t coeff)
        {
            std::vector<coeff_t> coeffs;
            std::vector<exp_t> exps;
            for (const auto & [e, c] : poly)
            {
                coeffs.push_back(coeff * c);
                exps.push_back(exp + e);
            }
            ExpPoly res(exps, coeffs);
            return res;
        }

        std::vector<coeff_t> get_coeffs() 
        {
            std::vector<coeff_t> coeffs;
            for (const auto& [exp, coeff]: poly)
            {
                coeffs.push_back(coeff);
            }
            return coeffs;
        }

        std::vector<exp_t> get_exps()
        {
            std::vector<exp_t> exps;
            for (const auto &[exp, coeff] : poly)
            {
                exps.push_back(exp);
            }
            return exps;
        }

        std::unordered_map<exp_t, coeff_t> get_poly()
        {
            return poly;
        }
    };

    ExpPoly operator*(ExpPoly::coeff_t c, ExpPoly &exp_poly)
    {
        ExpPoly res(exp_poly);
        for (auto [k, v] : res.poly)
        {
            res.poly[k] *= c;
        }
        return res;
    }

    ExpPoly operator*(ExpPoly &exp_poly, ExpPoly::coeff_t c)
    {
        return c * exp_poly;
    }

    std::ostream &operator<<(std::ostream &out, ExpPoly &exp_poly)
    {
        int i = 0;
        for (auto [exp, coeff] : exp_poly.poly)
        {
            ++i;
            if ((coeff == 1) && (exp == 1))
            {
                out << "e^x";
            }
            else if ((coeff == 1) && (exp != 1))
            {
                out << "e^" << exp << "x";
            }
            else if (((coeff != 1) && (exp == 1)))
            {
                out << coeff << "*e^x";
            }
            else
            {
                out << coeff << "*e^" << exp << "x";
            }
            if (i != exp_poly.poly.size()) 
            {
                std::cout << "+";
            }
        }
        return out;
    }

} // namespace GroupIP
