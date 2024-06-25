#pragma once
#include "ckks_evaluator.cuh"
#include "phantom.h"

using namespace std;
using namespace phantom;
using namespace nexus;

namespace nexus {
class GELUEvaluator {
 private:
  CKKSEvaluator *ckks = nullptr;

 public:
  GELUEvaluator(CKKSEvaluator &ckks) : ckks(&ckks) {}
  void gelu(PhantomCiphertext &x, PhantomCiphertext &res);

  vector<double> gelu_plain(vector<double> &input);
};
}  // namespace nexus
