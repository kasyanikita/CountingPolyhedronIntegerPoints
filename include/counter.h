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

mpf_class PolyEvaluation(Matrix &A, Vector &b, Vector &c);

#endif  // COUNTINGINTEGERPOINTS_COUNTER_H_