#include "group_element.h"

namespace GroupIP {

int_t modulo(int_t x, int_t mod) {
  auto res(x % mod);
  if (res < 0) res += mod;
  return res;
}

void GroupElement::normalize_components() {
  assert(_components.size() == _mod.size());
  for (size_t i = 0; i < _components.size(); ++i) {
    _components[i] = modulo(_components[i], _mod[i]);
  }
}

bool GroupElement::compare(const GroupElement& ge) {
  if (_components == ge._components && _mod == ge._mod) return true;
  return false;
}

void GroupElement::assign(std::vector<int_t> comp) {
  assert(comp.size() == _mod.size());
  _components = std::move(comp);
  normalize_components();
}

GroupElement GroupElement::invert() const {
  GroupElement res(*this);

  for (size_t i = 0; i < res._components.size(); ++i) {
    res._components[i] = -res._components[i];
  }
  res.normalize_components();

  return res;
}

bool GroupElement::operator==(const GroupElement& ge) { return compare(ge); }

bool GroupElement::operator!=(const GroupElement& ge) { return !compare(ge); }

GroupElement GroupElement::operator-(const GroupElement& rhs) const {
  return *this + rhs.invert();
}

GroupElement& GroupElement::operator+=(const GroupElement& rhs) {
  assert(_components.size() == rhs._components.size());
  assert(_mod == rhs._mod);

  for (size_t i = 0; i < _components.size(); ++i) {
    _components[i] += rhs._components[i];
  }
  normalize_components();

  return *this;
}

GroupElement GroupElement::operator+(const GroupElement& rhs) const {
  GroupElement res(*this);
  res += rhs;
  return res;
}

bool GroupElement::operator<(const GroupElement& rhs) const {
  assert(_components.size() == rhs._components.size());
  assert(_mod == rhs._mod);
  return std::lexicographical_compare(_components.begin(), _components.end(),
                                      rhs._components.begin(),
                                      rhs._components.end());
}

size_t GroupElement::get_idx() const {
  size_t res = 0;
  int_t mult = 1;
  for (size_t i = 0; i < _mod.size(); ++i) {
    res += static_cast<size_t>(_components[i] * mult);
    mult *= _mod[i];
  }
  return res;
}

int_t GroupElement::getOrder() {
  if (order != -1) {
    return order;
  }
  
  std::vector<int_t> e(_components.size(), 0);
  int r = 1;
  auto sum = *this;
  while (sum.get_components() != e) {
    sum += *this;
    ++r;
  }
  order = r;

  return order;
}

const std::vector<int_t>& GroupElement::get_components() const { return _components; }

const std::vector<int_t>& GroupElement::get_mod() const { return _mod; }

GroupElement operator*(int_t c, const GroupElement& x) {
  GroupElement res(x);

  for (size_t i = 0; i < res._components.size(); ++i) {
    res._components[i] *= c;
  }
  res.normalize_components();

  return res;
}

GroupElement operator*(const GroupElement& x, int_t c) { return c * x; }
}