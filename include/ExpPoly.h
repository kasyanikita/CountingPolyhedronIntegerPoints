#pragma once

#include "global_defs.h"

namespace GroupIP
{
    class ExpPoly
    {
        std::vector<double> _base;
        std::vector<double> _exponent;

    public:
        friend ExpPoly operator*(double c, ExpPoly &x);
        friend ExpPoly operator*(ExpPoly &x, double c);
        ExpPoly(const std::vector<double> &exp) : _exponent(exp), _base(exp.size(), 1)
        {
        }
        ExpPoly operator+(ExpPoly &rhs) const
        {
            assert(_exponent.size() == rhs._exponent.size());
        }

        const std::vector<double>& get_base() const {
            return _base;
        }

        const std::vector<double>& get_exponent() const {
            return _exponent;
        }
    };

    ExpPoly operator*(double c, ExpPoly &poly)
    {
        ExpPoly res(poly._exponent);
        for (auto &x : res._base)
        {
            x *= c;
        }
        return res;
    }

    ExpPoly operator*(ExpPoly &poly, double c)
    {
        return c * poly;
    }

} // namespace GroupIP