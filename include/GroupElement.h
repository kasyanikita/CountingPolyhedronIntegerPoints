#ifndef INCLUDE_GROUPELEMENT_H_
#define INCLUDE_GROUPELEMENT_H_

#include <vector>
#include <cassert>

class GroupElement {
    std::vector<int> components;
    std::vector<int> mod;
    void reduce_components() {
        for(size_t i = 0; i < components.size(); ++i) {
            components[i] %= mod[i];
        }
    }
public:
    GroupElement(const std::vector<int>& comp, const std::vector<int>& m): components(comp), mod(m) {
        assert(components.size() == mod.size());
        reduce_components();
    }

    GroupElement operator+(const GroupElement& rhs) const {
        assert(components.size() == rhs.components.size());
        assert(mod == rhs.mod);
        std::vector<int> res_comp;
        for (size_t i = 0; i < components.size(); ++i) {
            res_comp.push_back(components[i] + rhs.components[i]);
        }
        GroupElement result(res_comp, mod);
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

    const std::vector<int>& get_components() const {
        return components;
    }

    const std::vector<int>& get_mod() const {
        return mod;
    }
};

#endif  // INCLUDE_GROUPELEMENT_H_
