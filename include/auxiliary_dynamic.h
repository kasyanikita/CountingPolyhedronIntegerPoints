#ifndef AUXILIARY_DYNAMICS_H_
#define AUXILIARY_DYNAMICS_H_

#include "snf_class.h"
#include "group_element.h"
#include "Tools.h"

class AuxiliaryDynamic {
 public:
  AuxiliaryDynamic(const std::vector<std::vector<int_t>>&,
                   const std::vector<int_t>&, SNFClass&);
  void Init();
  const std::vector<GroupIP::GroupElement>& GetG() const;
  const std::vector<std::vector<int_t>>& GetH() const;

 private:
  std::vector<std::vector<int_t>> A_;
  std::vector<int_t> b_;
  SNFClass& snf_;
  std::vector<GroupIP::GroupElement> g_;
  std::vector<std::vector<int_t>> h_;

  void CalcG();
  void CalcH();
};

#endif  // AUXILIARY_DYNAMICS_H_
