#ifndef COUNTINGINTEGERPOINTS_GROUPELEMENT_H_
#define COUNTINGINTEGERPOINTS_GROUPELEMENT_H_

#include "global_defs.h"

namespace GroupIP
{
    inline int_t modulo(int_t x, int_t mod);

    class GroupElement
    {
        friend GroupElement operator*(int_t c, const GroupElement &x);
        friend GroupElement operator*(const GroupElement &x, int_t c);

        std::vector<int_t> _components;
        std::vector<int_t> _mod;
        int_t order = -1;

        void normalize_components();
        bool compare(const GroupElement& ge);

    public:
        GroupElement(std::vector<int_t> mod) : _mod(std::move(mod)) {}
        void assign(std::vector<int_t> comp);
        GroupElement invert() const;
        bool operator==(const GroupElement& ge);
        bool operator!=(const GroupElement &ge);
        GroupElement operator-(const GroupElement &rhs) const;
        GroupElement &operator+=(const GroupElement &rhs);
        GroupElement operator+(const GroupElement &rhs) const;
        bool operator<(const GroupElement &rhs) const;
        int_t getOrder();
        size_t get_idx() const;
        const std::vector<int_t> &get_components() const;
        const std::vector<int_t> &get_mod() const;
    };

    GroupElement operator*(int_t c, const GroupElement &x);
    GroupElement operator*(const GroupElement &x, int_t c);

} // namespace GroupIP

#endif  // COUNTINGINTEGERPOINTS_GROUPELEMENT_H_
