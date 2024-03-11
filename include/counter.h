#ifndef COUNTINGINTEGERPOINTS_COUNTER_H_
#define COUNTINGINTEGERPOINTS_COUNTER_H_

#include <gmpxx.h>

#include <chrono>
#include <iostream>

#include "Dynamic.h"
#include "ExpPoly.h"
#include "ToddPoly_FFT.h"
#include "Tools.h"
#include "flint/fmpz_mat.h"
#include "hyperplane_avoid_solver.h"

class Counter {
 private:
  Matrix A_;
  Vector b_;
  Vector c_;
  mpf_class count_simple_cone(Matrix &A, Vector &b, Vector &c);

 public:
  Counter(Matrix &A, Vector &b);
  mpf_class count_integer_points(Matrix &A, Vector &b);
};

#endif  // COUNTINGINTEGERPOINTS_COUNTER_H_