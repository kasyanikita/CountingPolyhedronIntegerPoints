#ifndef COUNTINGINTEGERPOINTS_CVECTOR_H_
#define COUNTINGINTEGERPOINTS_CVECTOR_H_

#include "global_defs.h"
#include "Tools.h"


namespace GroupIP {

    class HyperplaneAvoidSolver {
      private:
        Vector c_;
        Matrix A_;

      public:
        HyperplaneAvoidSolver(const Matrix& A) : A_(A) {}
        Vector get_vector(int_t);
    };
}

#endif  // COUNTINGINTEGERPOINTS_CVECTOR_H_