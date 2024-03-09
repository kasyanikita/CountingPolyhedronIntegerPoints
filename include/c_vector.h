#ifndef COUNTINGINTEGERPOINTS_CVECTOR_H_
#define COUNTINGINTEGERPOINTS_CVECTOR_H_

#include "global_defs.h"
#include "Tools.h"


namespace GroupIP {

    class CVector {
      private:
        Vector c_;
        Matrix A_;

      public:
        CVector(const Matrix& A): A_(A) {}
        Vector get_c_vector(int_t);
    };
}

#endif  // COUNTINGINTEGERPOINTS_CVECTOR_H_