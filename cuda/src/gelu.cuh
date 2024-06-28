#pragma once
#include "ckks_evaluator.cuh"
#include "phantom.h"

namespace nexus {
using namespace std;
using namespace phantom;

class GELUEvaluator {
 private:
  int d_g = 2;
  int d_f = 2;

  CKKSEvaluator *ckks = nullptr;

 public:
  GELUEvaluator(CKKSEvaluator &ckks) : ckks(&ckks) {}

  void gelu(PhantomCiphertext &x, PhantomCiphertext &res);
};
}  // namespace nexus
