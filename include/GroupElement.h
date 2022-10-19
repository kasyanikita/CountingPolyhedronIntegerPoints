#ifndef INCLUDE_GROUPELEMENT_H_
#define INCLUDE_GROUPELEMENT_H_

#include <vector>
#include <algorithm>
#include <cassert>

namespace GroupNS {

    using intT = int64_t;
    using uintT = uint64_t;

    inline intT modulo(intT a, intT b)
    {
        if (b < 0) return modulo(a, -b);
        const int result = a % b;
        return result >= 0 ? result : result + b;
    }

    class GroupElement {
        std::vector<intT> components;
        std::vector<intT> mod;
        void reduce_components() {
            for(size_t i = 0; i < components.size(); ++i) {
                components[i] = modulo(components[i], mod[i]);
            }
        }
    public:
        GroupElement(const std::vector<intT>& m): mod(m) { }

        void assign(const std::vector<intT>& comp) {
            assert(comp.size() == mod.size());
            components = comp;
            reduce_components();
        }

        GroupElement operator-(const GroupElement& rhs) const {
            assert(components.size() == rhs.components.size());
            assert(mod == rhs.mod);
            std::vector<intT> res_comp;
            for (size_t i = 0; i < components.size(); ++i)
            {
                res_comp.push_back(components[i] - rhs.components[i]);
            }
            GroupElement result(mod);
            result.assign(res_comp);
            return result;
        }

        GroupElement operator+(const GroupElement& rhs) const {
            assert(components.size() == rhs.components.size());
            assert(mod == rhs.mod);
            std::vector<intT> res_comp;
            for (size_t i = 0; i < components.size(); ++i) {
                res_comp.push_back(components[i] + rhs.components[i]);
            }
            GroupElement result(mod);
            result.assign(res_comp);
            return result;
        }

        GroupElement& operator+=(const GroupElement& rhs) {
            assert(components.size() == rhs.components.size());
            assert(mod == rhs.mod);
            for(size_t i = 0; i < components.size(); ++i) {
                components[i] += rhs.components[i];
            }
            reduce_components();
            return *this;
        }

        bool operator<(const GroupElement& rhs) const {
            assert(components.size() == rhs.components.size());
            assert(mod == rhs.mod);
            return std::lexicographical_compare(components.begin(), components.end(),
                                                rhs.components.begin(), rhs.components.end());
        }

        size_t get_id() const {
            size_t res = 0;
            intT mult = 1;
            for (size_t i = 0; i < mod.size(); ++i)
            {
                res += static_cast<size_t>(components[i] * mult);
                mult *= mod[i];
            }
            return res;
        }

        const std::vector<intT>& get_components() const {
            return components;
        }

        const std::vector<intT>& get_mod() const {
            return mod;
        }
    };

    GroupElement operator*(intT x, const GroupElement &elem) {
        GroupElement res(elem.get_mod());
        auto comp = elem.get_components();
        std::vector<intT> res_comp;
        for (size_t i = 0; i < comp.size(); ++i) {
            res_comp.push_back(x * comp[i]);
        }
        res.assign(res_comp);
        return res;
    }

    GroupElement operator*(const GroupElement &elem, intT x) {
        return x * elem;
    }
}

#endif  // INCLUDE_GROUPELEMENT_H_
