#pragma once

#include "global_defs.h"

namespace GroupIP
{
    inline int_t modulo(int_t x, int_t mod)
    {
        auto res(x % mod);
        if (res < 0)
            res += mod;
        return res;
    }

    class GroupElement
    {
        friend GroupElement operator*(int_t c, const GroupElement &x);
        friend GroupElement operator*(const GroupElement &x, int_t c);

        std::vector<int_t> _components;
        std::vector<int_t> _mod;

        void normalize_components()
        {
            assert(_components.size() == _mod.size());
            for (size_t i = 0; i < _components.size(); ++i)
            {
                _components[i] = modulo(_components[i], _mod[i]);
            }
        }

        bool compare(const GroupElement& ge)
        {
            if (_components == ge._components && _mod == ge._mod)
                return true;
            return false;
        }

    public:
        GroupElement(std::vector<int_t> mod) : _mod(std::move(mod)) {}

        void assign(std::vector<int_t> comp)
        {
            assert(comp.size() == _mod.size());
            _components = std::move(comp);
            normalize_components();
        }

        GroupElement invert() const
        {
            GroupElement res(*this);

            for (size_t i = 0; i < res._components.size(); ++i)
            {
                res._components[i] = -res._components[i];
            }
            res.normalize_components();

            return res;
        }

        bool operator==(const GroupElement& ge)
        {
            return compare(ge);
        }

        bool operator!=(const GroupElement &ge)
        {
            return !compare(ge);
        }

        GroupElement operator-(const GroupElement &rhs) const
        {
            return *this + rhs.invert();
        }

        GroupElement &operator+=(const GroupElement &rhs)
        {
            assert(_components.size() == rhs._components.size());
            assert(_mod == rhs._mod);

            for (size_t i = 0; i < _components.size(); ++i)
            {
                _components[i] += rhs._components[i];
            }
            normalize_components();

            return *this;
        }

        GroupElement operator+(const GroupElement &rhs) const
        {
            GroupElement res(*this);
            res += rhs;
            return res;
        }

        bool operator<(const GroupElement &rhs) const
        {
            assert(_components.size() == rhs._components.size());
            assert(_mod == rhs._mod);
            return std::lexicographical_compare(_components.begin(), _components.end(),
                                                rhs._components.begin(), rhs._components.end());
        }

        size_t get_idx() const
        {
            size_t res = 0;
            int_t mult = 1;
            for (size_t i = 0; i < _mod.size(); ++i)
            {
                res += static_cast<size_t>(_components[i] * mult);
                mult *= _mod[i];
            }
            return res;
        }

        const std::vector<int_t> &get_components() const
        {
            return _components;
        }

        const std::vector<int_t> &get_mod() const
        {
            return _mod;
        }
    };

    GroupElement operator*(int_t c, const GroupElement &x)
    {
        GroupElement res(x);

        for (size_t i = 0; i < res._components.size(); ++i)
        {
            res._components[i] *= c;
        }
        res.normalize_components();

        return res;
    }

    GroupElement operator*(const GroupElement &x, int_t c)
    {
        return c * x;
    }

} // namespace GroupIP
