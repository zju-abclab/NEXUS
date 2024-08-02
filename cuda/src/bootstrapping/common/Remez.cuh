#pragma once

#include <NTL/mat_RR.h>

#include <future>
#include <string>
#include <vector>

#include "Point.cuh"
#include "Polynomial.cuh"
#include "RemezParam.h"
#include "func.cuh"
#define _USE_MATH_DEFINES

using namespace std;
using namespace NTL;

namespace boot {
class Remez {
 private:
  RR width, sc, approx, max_err, min_err, current_err;
  Point *sample_point, *extreme_point;
  long extreme_count;
  RR *coeff;

 public:
  RemezParam params;
  long boundary_K;
  double log_width;
  long deg;

  Remez(RemezParam _params, long _boundary_K, double _log_width, long _deg);

  virtual RR function_value(RR x) = 0;
  void better_initialize();
  void initialize();
  void getcoeffwitherr();
  void getextreme_local(Point *local_extreme_point, long &local_extreme_count, long k);
  void getextreme();
  void choosemaxs();

  void generate_optimal_poly(Polynomial &poly);
  void showcoeff();
};
}  // namespace boot
